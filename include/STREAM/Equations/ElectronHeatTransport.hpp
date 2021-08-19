#ifndef _STREAM_ELECTRON_HEAT_TRANSPORT_HPP
#define _STREAM_ELECTRON_HEAT_TRANSPORT_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"

namespace STREAM {
    class ElectronHeatTransport : public DREAM::FVM::EquationTerm {
    private:
        ConfinementTime *confinementTime;
        
        real_t invtau;

        EllipticalRadialGridGenerator *radials;

        // Precomputed coefficient used for calculating
        // derivatives of the diffusion coefficient Drr...
        real_t dI_p, dI_wall, dT_cold, dW_i, dN_i;

        // IDs of unknown quantities used by the operator...
        len_t id_Ip, id_Iwall=0, id_Tcold, id_Wi, id_Ni, id_Wcold;

    public:
        ElectronHeatTransport(
            DREAM::FVM::Grid*, ConfinementTime*,
            EllipticalRadialGridGenerator*, DREAM::FVM::UnknownQuantityHandler*);
        ~ElectronHeatTransport();

        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 6; }
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*);

        virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*);
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*);
        virtual void SetVectorElements(real_t*, const real_t*);
    };
}

#endif/*_STREAM_ELECTRON_HEAT_TRANSPORT_HPP*/
