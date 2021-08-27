#ifndef _STREAM_RUNAWAY_ELECTRON_TRANSPORT_HPP
#define _STREAM_RUNAWAY_ELECTRON_TRANSPORT_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/RunawayElectronConfinementTime.hpp"

namespace STREAM {
    class RunawayElectronTransport : public DREAM::FVM::EquationTerm {
    private:
        EllipticalRadialGridGenerator *radials;
        RunawayElectronConfinementTime *reConfinementTime;

        real_t invtau;
        real_t dI_p, dI_wall, dE_field;

        len_t id_Ip, id_Iwall, id_Efield;

    public:
        RunawayElectronTransport(
            DREAM::FVM::Grid*, EllipticalRadialGridGenerator*,
            RunawayElectronConfinementTime*, DREAM::FVM::UnknownQuantityHandler*
        );
        ~RunawayElectronTransport();

        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 4; }
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*);

        virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*);
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*);
        virtual void SetVectorElements(real_t*, const real_t*);
    };
}

#endif/*_STREAM_RUNAWAY_ELECTRON_TRANSPORT_HPP*/
