#ifndef _STREAM_EQUATION_RADIATED_POWER_TERM_HPP
#define _STREAM_EQUATION_RADIATED_POWER_TERM_HPP

#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/AMJUEL.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

namespace STREAM{
    class RadiatedPowerTerm : public DREAM::RadiatedPowerTerm{ 
    private:
        DREAM::FVM::UnknownQuantityHandler *unknowns;
        DREAM::IonHandler *ions;
        PlasmaVolume *volumes;
        len_t id_T_cold, id_n_cold, id_lambda_i;
    protected:
        void SetWeights() override;
        void SetDiffWeights(len_t derivId, len_t nMultiples) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
    public:
        RadiatedPowerTerm(DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*, DREAM::IonHandler*, DREAM::ADAS*, DREAM::NIST*, DREAM::AMJUEL*, DREAM::OptionConstants::ion_opacity_mode*, bool, PlasmaVolume*);
    };
}

#endif /*_STREAM_EQUATION_RADIATED_POWER_TERM_HPP*/
