#ifndef _STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP
#define _STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace STREAM {
    class ElectronTransportDiffusion : public FVM::DiffusionTerm { //Ska vara subklass?
    private:
        enum DREAM::OptionConstants::momentumgrid_type mgtype;
        DREAM::FVM::Interpolator1D *coefftauinv; // Rätt implementation av 1/tau?

        DREAM::FVM::UnknownQuantityHandler *unknowns;
        
        EllipticalRadialGridGenerator *radials;

        // Precomputed coefficient used for calculating
        // derivatives of the diffusion coefficient Drr...
        real_t *dtauinv=nullptr;

        // IDs of unknown quantities used by the operator...
        len_t id_n_cold;

        /* Dessa behövs va? */
        void AllocateDiffCoeff(); 
        virtual void SetPartialDiffusionTerm(len_t, len_t) override;

    public:
        ElectronTransportDiffusion(DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, EllipticalRadialGridGenerator*, DREAM::FVM::Interpolator1D*, DREAM::FVM::UnknownQuantityHandler*);
        ~ElectronTransportDiffusion();

        /* Dessa behövs va? */
        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_STREAM_EQUATIONS_ELECTRON_TRANSPORT_DIFFUSION_HPP*/
