#ifndef _STREAM_EQUATION_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATION_CONFINEMENT_TIME_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <gsl/gsl_math.h>
#include "DREAM/Constants.hpp"



/** 
 * Implementation of the particle confinement time
 */

namespace STREAM{
        class ConfinementTime {

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; 
                        
                        EllipticalRadialGridGenerator *radials;
                        
                        len_t id_Ip, id_Imk2, id_Tcold, id_Wi, id_ni; 

                public:
                        real_t l_MK2;
                        real_t I_ref = 100000;
                        real_t B_v   = 1e-3;
                
                        ConfinementTime(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *aB, real_t l_MK2);
                        
                        real_t EvaluateConfinementTime(len_t ir, real_t t);

                        real_t EvaluateConfinementTime_dIp(len_t ir, real_t t);

                        real_t EvaluateConfinementTime_dIMK2(len_t ir, real_t t);

                        real_t EvaluateConfinementTime_dTe(len_t ir, real_t t);

                        real_t EvaluateConfinementTime_dWi(len_t ir, real_t t);
                        
                        real_t EvaluateConfinementTime_dni(len_t ir, real_t t);
}

#endif/*_STREAM_EQUATION_CONFINEMENT_TIME_HPP */
