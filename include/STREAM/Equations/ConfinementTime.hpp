#ifndef _STREAM_EQUATIONS_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATIONS_CONFINEMENT_TIME_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <gsl/gsl_math.h>
#include "DREAM/Constants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/ConnectionLength.hpp"



/** 
 * Implementation of the particle confinement time
 */

namespace STREAM{
        class ConfinementTime {

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; 
                        
                        EllipticalRadialGridGenerator *radials;
						DREAM::IonHandler *ions;
                        
                        len_t id_Tcold, id_Wi, id_Ni; 
                        
                        ConnectionLength *CL;
                        
                        len_t D_index;

                public:
                        ConfinementTime(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
					DREAM::IonHandler*, ConnectionLength*, len_t D_index=0);
                        
                        real_t EvaluateConfinementTime(len_t ir);

                        real_t EvaluateParallelConfinementTime(len_t ir);

                        real_t EvaluatePerpendicularConfinementTime(len_t ir);

                        real_t EvaluateConfinementTime_dIp(len_t ir);

                        real_t EvaluateConfinementTime_dIwall(len_t ir);

                        real_t EvaluateConfinementTime_dTcold(len_t ir);

                        real_t EvaluateConfinementTime_dWi(len_t ir);
                        
                        real_t EvaluateConfinementTime_dNi(len_t ir);
    };
}

#endif/*_STREAM_EQUATIONS_CONFINEMENT_TIME_HPP */
