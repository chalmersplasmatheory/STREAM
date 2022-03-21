#ifndef _STREAM_ELECTRON_CYCLOTRON_HEATING_HPP
#define _STREAM_ELECTRON_CYCLOTRON_HEATING_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "STREAM/Equations/OpticalThickness.hpp"

namespace STREAM {
    class ElectronCyclotronHeating : public DREAM::FVM::EquationTerm {
    private:
        EllipticalRadialGridGenerator *radials;
        OpticalThickness *OT;
	
	real_t parentheses_ECH, dpECH_dTe, dpECH_dne;

        len_t id_Tcold, id_ncold;
        
        // Input parameters
        real_t P_inj, f_o, f_x, theta; 

    public:
        ElectronCyclotronHeating(
            DREAM::FVM::Grid*,
            EllipticalRadialGridGenerator*, DREAM::FVM::UnknownQuantityHandler*, OpticalThickness*, real_t, real_t, real_t, real_t);
        ~ElectronCyclotronHeating();

        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; } 
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 2; } 
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*); 

        virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*);
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*); 
        virtual void SetVectorElements(real_t*, const real_t*);
    };
}

#endif/*_STREAM_ELECTRON_CYCLOTRON_HEATING_HPP*/
