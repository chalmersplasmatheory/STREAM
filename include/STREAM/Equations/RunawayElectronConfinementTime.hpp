#ifndef _STREAM_EQUATIONS_RUNAWAY_ELECTRON_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATIONS_RUNAWAY_ELECTRON_CONFINEMENT_TIME_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <gsl/gsl_math.h>
#include "DREAM/Constants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"



/** 
 * Implementation of the runaway electron confinement time
 */

namespace STREAM{
        class RunawayElectronConfinementTime {

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; 
                        
                        EllipticalRadialGridGenerator *radials;
                        
                        len_t id_Ip, id_Iwall; 

                public:
                        real_t l_MK2;
                        real_t I_ref = 100000;
                        real_t B_v   = 1e-3;
                
                        ConfinementTime(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, real_t l_MK2);
                        
                        real_t EvaluateConfinementTime(len_t ir);

                        real_t EvaluateConfinementTime_dIp(len_t ir);

                        real_t EvaluateConfinementTime_dIwall(len_t ir);
    };
}

#endif/*_STREAM_EQUATIONS_CONFINEMENT_TIME_HPP */
