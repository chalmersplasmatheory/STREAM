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
	
	real_t parentheses_ECH, eta_polo, eta_polx;
        real_t dT_cold, dn_cold;

        len_t id_Tcold, id_ncold;
        
        // Input parameters
        real_t P_inj; 

    public:
        ElectronCyclotronHeating(
            DREAM::FVM::Grid*,
            EllipticalRadialGridGenerator*, DREAM::FVM::UnknownQuantityHandler*);
        ~ElectronCyclotronHeating();

        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; } // What should return?
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 6; } // What should return?
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*); // beräkna parantes och optical thickness

        virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*); // Derivator
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*); // Använda från SetJacobianBlock?
        virtual void SetVectorElements(real_t*, const real_t*); // Använda från Rebuilt
    };
}

#endif/*_STREAM_ELECTRON_CYCLOTRON_HEATING_HPP*/
