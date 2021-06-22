/**
 * Construct equation for the ion densities.
 */

#include <vector>
#include "DREAM/Equations/Fluid/IonTransientTerm.hpp"
#include "DREAM/Equations/Fluid/IonPrescribedParameter.hpp"
#include "DREAM/Equations/Fluid/LyOpaqueDIonRateEquation.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesIdentityTerm.hpp"
#include "STREAM/Equations/IonTransport.hpp"


using namespace STREAM;
using namespace std;


#define MODULENAME "eqsys/n_i"


/**
 * Define options for the ion module.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_Ions(DREAM::Settings *s) {
    // Include all settings from DREAM...
    DREAM::SimulationGenerator::DefineOptions_Ions(s);

    // TODO transport settings
}

/**
 * Construct the equation governing the evolution of the
 * ion densities for each charge state.
 */
void SimulationGenerator::ConstructEquation_Ions(
    EquationSystem *eqsys, DREAM::Settings *s,
    DREAM::ADAS *adas, DREAM::AMJUEL *amjuel
) {
    const real_t t0 = 0;
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    len_t nZ, ntypes;
    const int_t *_Z  = s->GetIntegerArray(MODULENAME "/Z", 1, &nZ);
    const int_t *itypes = s->GetIntegerArray(MODULENAME "/types", 1, &ntypes);
    const int_t *iopacity_modes = s->GetIntegerArray(MODULENAME "/opacity_modes", 1, &ntypes);

    // Parse list of ion names (stored as one contiguous string,
    // each substring separated by ';')
    vector<string> ionNames = s->GetStringList(MODULENAME "/names");

    // Automatically name any unnamed ions
    if (ionNames.size() < nZ) {
        for (len_t i = ionNames.size(); i < nZ; i++)
            ionNames.push_back("Ion " + to_string(i));
    } else if (ionNames.size() > nZ) {
        throw DREAM::SettingsException(
            "ions: Too many ion names given: %zu. Expected " LEN_T_PRINTF_FMT ".",
            ionNames.size(), nZ
        );
    }

    // Get list of tritium species
    vector<string> tritiumNames = s->GetStringList(MODULENAME "/tritiumnames");

    // Verify that exactly one type per ion species is given
    if (nZ != ntypes)
        throw DREAM::SettingsException(
            "ions: Expected the lengths of 'Z' and 'types' to match."
        );

    // Data type conversion
    len_t *Z = new len_t[nZ];
    for (len_t i = 0; i < nZ; i++)
        Z[i] = (len_t)_Z[i];

    enum DREAM::OptionConstants::ion_data_type *types =
        new enum DREAM::OptionConstants::ion_data_type[ntypes];
    for (len_t i = 0; i < ntypes; i++)
        types[i] = (enum DREAM::OptionConstants::ion_data_type)itypes[i];

    // Verify that all non-prescribed elements are in ADAS
    for (len_t i = 0; i < nZ; i++) {
        if (!adas->HasElement(Z[i]) && types[i] != DREAM::OptionConstants::ION_DATA_PRESCRIBED)
            throw DREAM::SettingsException(
                "ions: The DREAM ADAS database does not contain '%s' (Z = " LEN_T_PRINTF_FMT ")",
                ionNames[i].c_str(), Z[i]
            );
    }
    
    enum DREAM::OptionConstants::ion_opacity_mode *opacity_mode = new enum DREAM::OptionConstants::ion_opacity_mode[ntypes];
    for (len_t i = 0; i < ntypes; i++)
        opacity_mode[i] = (enum DREAM::OptionConstants::ion_opacity_mode)iopacity_modes[i];
        
    // Verify that all non-prescribed elements have ground state opaque coefficients available
    for (len_t i = 0; i < nZ; i++) {
        if (Z[i]!=1 && opacity_mode[i]==DREAM::OptionConstants::OPACITY_MODE_GROUND_STATE_OPAQUE){
            throw DREAM::SettingsException(
            	"ions: There are no rate coefficients implemented for plasmas opaque to radiative transitions to the ground state for other species than hydrogen isotopes"
            );
        }
    }

    /////////////////////
    /// LOAD ION DATA ///
    /////////////////////
    // Count number of prescribed/dynamic charge states
    len_t nZ0_prescribed=0, nZ_prescribed=0, nZ_dynamic=0, nZ0_dynamic=0;
    len_t *prescribed_indices = new len_t[nZ];
    len_t *dynamic_indices = new len_t[nZ];
    for (len_t i = 0; i < nZ; i++) {
        switch (types[i]) {
            case DREAM::OptionConstants::ION_DATA_PRESCRIBED:
                nZ0_prescribed += Z[i] + 1;
                prescribed_indices[nZ_prescribed++] = i;
                break;

            case DREAM::OptionConstants::ION_DATA_TYPE_DYNAMIC:
            case DREAM::OptionConstants::ION_DATA_EQUILIBRIUM:
                nZ0_dynamic += Z[i] + 1;
                dynamic_indices[nZ_dynamic++] = i;
                break;

            default:
                throw DREAM::SettingsException(
                    "ions: Unrecognized ion model type specified: %d.",
                    types[i]
                );
        }
    }

    // Load ion data
    real_t *dynamic_densities = DREAM::SimulationGenerator::LoadDataIonR(
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ0_dynamic, "initial"
    );
    DREAM::MultiInterpolator1D *prescribed_densities =
        DREAM::SimulationGenerator::LoadDataIonRT(
            MODULENAME, fluidGrid->GetRadialGrid(), s, nZ0_prescribed, "prescribed"
        );

    // Construct ion handler
    DREAM::IonHandler *ih = new DREAM::IonHandler(fluidGrid->GetRadialGrid(), eqsys->GetUnknownHandler(), Z, nZ, ionNames, tritiumNames);
    eqsys->SetIonHandler(ih);

    // Initialize ion equations
    DREAM::FVM::Operator *eqn = new DREAM::FVM::Operator(fluidGrid);

    DREAM::IonPrescribedParameter *ipp = nullptr;
    if (nZ0_prescribed > 0)
        ipp = new DREAM::IonPrescribedParameter(fluidGrid, ih, nZ_prescribed, prescribed_indices, prescribed_densities);

    const len_t id_ni = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES);
    // Construct dynamic equations
    len_t nDynamic = 0, nEquil = 0;
    for (len_t iZ = 0; iZ < nZ; iZ++) {
        switch (types[iZ]) {
            case DREAM::OptionConstants::ION_DATA_PRESCRIBED: 
                break;
            // 'Dynamic' and 'Equilibrium' differ by a transient term
            case DREAM::OptionConstants::ION_DATA_TYPE_DYNAMIC:
                nDynamic++;
                eqn->AddTerm(
                    new DREAM::IonTransientTerm(fluidGrid, ih, iZ, id_ni)
                );
                [[fallthrough]];
            case DREAM::OptionConstants::ION_DATA_EQUILIBRIUM:
                nEquil++;
                // TODO Update these
                if(ih->GetZ(iZ)==1 && opacity_mode[iZ]==DREAM::OptionConstants::OPACITY_MODE_GROUND_STATE_OPAQUE){
		            eqn->AddTerm(new DREAM::LyOpaqueDIonRateEquation(
		                fluidGrid, ih, iZ, eqsys->GetUnknownHandler(),
		                true, false, false, amjuel
		            ));		            
                }else{
		            eqn->AddTerm(new DREAM::IonRateEquation(
		                fluidGrid, ih, iZ, adas, eqsys->GetUnknownHandler(),
		                true, false, false
		            ));
                }
                break;

            default:
                throw DREAM::SettingsException(
                    "ions: Unrecognized ion model type specified: %d.",
                    types[iZ]
                );
        }
    }

    // Set equation description
    string desc;
    if (ipp != nullptr && nEquil > 0) {
        if (nEquil == nDynamic)
            desc = "Prescribed + dynamic";
        else
            desc = "Prescribed + dynamic + equilibrium";
    } else if (ipp != nullptr)
        desc = "Fully prescribed";
    else {
        if (nEquil == nDynamic)
            desc = "Fully dynamic";
        else if (nDynamic == 0)
            desc = "Fully equilibrium";
        else
            desc = "Dynamic + equilibrium";
    }

    if (ipp != nullptr)
        eqn->AddTerm(ipp);
        
    // Add ion equation to system
	eqsys->SetOperator(id_ni, id_ni, eqn, desc);

    // Initialize dynamic ions
    const len_t Nr = fluidGrid->GetNr();
    real_t *ni = new real_t[ih->GetNzs() * Nr];

    for (len_t i = 0; i < ih->GetNzs() * Nr; i++)
        ni[i] = 0;

    // Begin by evaluating prescribed densities
    if (ipp != nullptr) {
        ipp->Rebuild(t0, 1, nullptr);
        ipp->Evaluate(ni);
    }

    // ...and then fill in with the initial dynamic ion values
    for (len_t i = 0, ionOffset = 0; i < nZ_dynamic; i++) {
        len_t Z   = ih->GetZ(dynamic_indices[i]);
        len_t idx = ih->GetIndex(dynamic_indices[i], 0);

        for (len_t Z0 = 0; Z0 <= Z; Z0++) {
            for (len_t ir = 0; ir < Nr; ir++)
                ni[(idx+Z0)*Nr+ir] = dynamic_densities[ionOffset+ir];
            ionOffset += Nr;
        }
    }

    eqsys->SetInitialValue(id_ni, ni, t0);
    ih->Rebuild();

    delete [] types;
}

void SimulationGenerator::ConstructEquation_ion_transport(EquationSystem *eqsys, DREAM::Settings *s) {
    DREAM::FVM::Operator *op_ion_transport = new DREAM::FVM::Operator(eqsys->GetFluidGrid());
    DREAM::FVM::Operator *op_n_i =  new DREAM::FVM::Operator(eqsys->GetFluidGrid());
    
    DREAM::IonHandler *ions = eqsys->GetIonHandler();
    
    for (len-t iz=1; iz<ions->GetNZ(); iz++) {// Ska vara iz=1 och inte 0?
        op_ion_transport->AddTerm(new DREAM::IonSpeciesIdentityTerm(eqsys->GetFluidGrid(), iz, -1.0));
        op_n_i->AddTerm(new IonTransport(eqsys->GetFluidGrid(), ions, iz, eqsys->GetConfinementTime(), eqsys->GetUnknownHandler()));
    }
    
    eqsys->SetOperator(OptionConstants::UQTY_ION_TRANSPORT, OptionConstants::UQTY_ION_TRANSPORT, op_ion_transport, "-n_i^(j)/tau_D");
    eqsys->SetOperator(OptionConstants::UQTY_ION_TRANSPORT, DREAM::OptionConstants::UQTY_ION_SPECIES, op_n_i);
}
