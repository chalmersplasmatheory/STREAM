#ifndef _STREAM_EQUATIONS_NEUTRAL_INFLUX_HPP
#define _STREAM_EQUATIONS_NEUTRAL_INFLUX_HPP

#include "STREAM/Equations/ConfinementTime.hpp"
#include <gsl/gsl_math.h>

namespace STREAM{
    class NeutralInflux {
        private: 
            DREAM::IonHandler *ions;
            SputteredRecycledCoefficient *SRC;
            ConfinementTime *coefftauinv; 
            /*PlasmaVolume *PV; ***Add when have PlasmaVolume class */

            real_t c1, c2, c3;
            
            real_t tauinv/*, V_p ***Add when have PlasmaVolume class */;

            real_t DeuteriumRecyclingCoefficient(real_t t);
        
        public:
            NeutralInflux(DREAM::IonHandler *ihdl, SputteredRecycledCoefficient *SRC, ConfinementTime *coefftauinv, /*PlasmaVolume *PV, ***Add when have PlasmaVolume class */ real_t c1, real_t c2, real_t c3);
            
            real_t EvaluateNeutralInflux_dni(real_t t, const len_t iIon);
                        
            real_t EvaluateNeutralInflux(real_t t, const len_t iIon);

            real_t EvaluateNeutralInflux_dIp(real_t t, const len_t iIon);

            real_t EvaluateNeutralInflux_dIwall(real_t t, const len_t iIon);

            real_t EvaluateNeutralInflux_dTe(real_t t, const len_t iIon);

            real_t EvaluateNeutralInflux_dWi(real_t t, const len_t iIon);
            
            real_t EvaluateNeutralInflux_dNi(real_t t, const len_t iIon);
    };
}
#endif /*_STREAM_EQUATIONS_NEUTRAL_INFLUX_HPP*/
