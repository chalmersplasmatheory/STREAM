#ifndef _STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP
#define _STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class ElectronTransportDiffusion : public FVM::DiffusionTerm {
    private:
        enum OptionConstants::momentumgrid_type mgtype;
        FVM::Interpolator1D *coeffD;

        FVM::UnknownQuantityHandler *unknowns;
        
        EllipticalRadialGridGenerator *radials;

        // Precomputed coefficient used for calculating
        // derivatives of the diffusion coefficient Drr...
        real_t *dD=nullptr;

        // IDs of unknown quantities used by the operator...
        len_t id_n_cold;
        len_t id_T_cold;
        
        real_t a;

        void AllocateDiffCoeff();
        virtual void SetPartialDiffusionTerm(len_t, len_t) override;

    public:
        ElectronTransportDiffusion(FVM::Grid*, enum OptionConstants::momentumgrid_type, EllipticalRadialGridGenerator *aB, ConfinementTime *tau, FVM::UnknownQuantityHandler*);
        ~ElectronTransportDiffusion();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_DIFFUSION_HPP*/
