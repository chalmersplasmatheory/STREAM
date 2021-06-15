#ifndef _STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP
#define _STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace STREAM {
    class ElectronTransportDiffusion : public DREAM::FVM::DiffusionTerm {
    private:
        enum DREAM::OptionConstants::momentumgrid_type mgtype;
        ConfinementTime *coefftauinv; 
        DREAM::FVM::UnknownQuantityHandler *unknowns;
        
        EllipticalRadialGridGenerator *radials;

        // Precomputed coefficient used for calculating
        // derivatives of the diffusion coefficient Drr...
        real_t *dI_p=nullptr;
        real_t *dI_wall=nullptr;
        real_t *dT_cold=nullptr;
        real_t *dW_i=nullptr;
        real_t *dn_i=nullptr;

        // IDs of unknown quantities used by the operator...
        len_t id_Ip, id_Iwall, id_Tcold, id_Wi, id_ni;

        void AllocateDiffCoeff(); 
        virtual void SetPartialDiffusionTerm(len_t, len_t) override;

    public:
        ElectronTransportDiffusion(DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, EllipticalRadialGridGenerator*, ConfinementTime, DREAM::FVM::UnknownQuantityHandler*);
        ~ElectronTransportDiffusion();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP*/
