#ifndef _STREAM_EQUATIONS_NEUTRAL_INFLUX_HPP
#define _STREAM_EQUATIONS_NEUTRAL_INFLUX_HPP

#include "DREAM/IonHandler.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/SputteredRecycledCoefficient.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp" 
#include <gsl/gsl_math.h>

namespace STREAM{ 
    class NeutralInflux {
        private: 
            DREAM::IonHandler *ions;
            SputteredRecycledCoefficient *SRC;
            ConfinementTime *coefftauinv; 
            PlasmaVolume *PV; 

        public:
            NeutralInflux(DREAM::IonHandler *ihdl, SputteredRecycledCoefficient *SRC, ConfinementTime *coefftauinv, PlasmaVolume *PV);
            
            real_t EvaluateNeutralInflux_dnkj(const len_t iIon, const len_t kIon);
                        
            real_t EvaluateNeutralInflux(const len_t iIon);

            real_t EvaluateNeutralInflux_dIp(const len_t iIon);

            real_t EvaluateNeutralInflux_dIwall(const len_t iIon);

            real_t EvaluateNeutralInflux_dTcold(const len_t iIon);

            real_t EvaluateNeutralInflux_dWi(const len_t iIon);
            
            real_t EvaluateNeutralInflux_dNi(const len_t iIon);
    };
}

#endif /*_STREAM_EQUATIONS_NEUTRAL_INFLUX_HPP*/
