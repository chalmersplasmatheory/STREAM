#include <gsl/gsl_interp.h>
#include <string>
#include "DREAM/MultiInterpolator1D.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


using namespace STREAM;
using namespace DREAM;
using namespace std;

/**
 * Define options for a temporal 'data' section in the
 * specified module.
 */
void STREAM::SimulationGenerator::DefineDataT_3D(
    const string& modname, Settings *s, const string& name
) {
    len_t dims[3] = {0};
    s->DefineSetting(modname + "/" + name + "/t", "Time grid on which the prescribed data is defined.", 0, (real_t*)nullptr);
    s->DefineSetting(modname + "/" + name + "/tinterp", "Interpolation method to use for time grid interpolation.", (int_t)DREAM::OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR);
    s->DefineSetting(modname + "/" + name + "/x", "Table where each element is the prescribed data.", 3, dims, (real_t*)nullptr);
}

/**
 * Load data from the 'data' seciton of the specified module.
 * The data is expected to depend on time only.
 *
 * modname: Name of module to load data from.
 * s:       Settings object to load data from.
 * name:    Name of group containing data structure (default: "data").
 */
FVM::Interpolator1D ***STREAM::SimulationGenerator::LoadDataT_3D(
    const string& modname, Settings *s, const string& name
) {
    len_t nx[3], nt;

    const real_t *t = s->GetRealArray(modname + "/" + name + "/t", 1, &nt);
    const real_t *x = s->GetRealArray(modname + "/" + name + "/x", 3, nx);

    if (nt != nx[2])
        throw SettingsException(
            "%s: Inconsistent dimensions of data. Data has "
            LEN_T_PRINTF_FMT " elements while the time vector has "
            LEN_T_PRINTF_FMT " elements.",
            (modname+"/"+name).c_str(), nx, nt
        );

    enum DREAM::OptionConstants::prescribed_data_interp tinterp =
        (enum DREAM::OptionConstants::prescribed_data_interp)s->GetInteger(modname + "/" + name + "/tinterp");

    // Select Interpolator1D interpolation method
    enum FVM::Interpolator1D::interp_method interp1_meth;
    switch (tinterp) {
        case DREAM::OptionConstants::PRESCRIBED_DATA_INTERP_NEAREST:
            interp1_meth = FVM::Interpolator1D::INTERP_NEAREST; break;
        case DREAM::OptionConstants::PRESCRIBED_DATA_INTERP_LINEAR:
            interp1_meth = FVM::Interpolator1D::INTERP_LINEAR; break;

        default:
            throw SettingsException(
                "%s: Unrecognized interpolation method on time grid: %d.",
                modname.c_str(), tinterp
            );
    }
    printf("nx0=%lu\n", nx[0]);
    FVM::Interpolator1D ***interp = new FVM::Interpolator1D**[nx[0]];
    interp[0] = new FVM::Interpolator1D*[nx[0]*nx[1]];
    
    bool onlyZeroElements;
    for (len_t i = 0; i < nx[0]; i++) {
        if (i > 0)
            interp[i] = interp[i-1] + nx[1];
    	for (len_t j = 0; j < nx[1]; j++) {
    	    onlyZeroElements = true;
    	    
    	    for (len_t k = 0; k < nt; k++) {
		if (x[(i*nx[1]+j)*nx[2] + k] != 0) {
		    onlyZeroElements = false;
		    break;
		}
		
	    }
	    if (onlyZeroElements)
	        interp[i][j] = nullptr;
	    else {
	        real_t *new_x = new real_t[nt];
       	        real_t *new_t = new real_t[nt];
	    
   	        // Copy data
  	        for (len_t k = 0; k < nt; k++) {
		    new_t[k] = t[k];
  		    new_x[k] = x[(i*nx[1]+j)*nx[2] + k];
	        }
	    
    	        interp[i][j] = new FVM::Interpolator1D(nt, 1, new_t, new_x, interp1_meth);
    	    }
    	}
    }
    
    return interp;
}
