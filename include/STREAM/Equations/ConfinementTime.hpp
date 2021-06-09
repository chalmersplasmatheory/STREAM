#ifndef _STREAM_EQUATION_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATION_CONFINEMENT_TIME_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
/* Behövs båda dessa? Är adresserna rätt?  */

/** 
 * Implementation of the particle confinement time
 */

namespace STREAM{
        class ConefinementTime : public FVM::EvaluableEquationTerm {

                protected:
                        FVM::UnknownQuantityHandler *unknowns; /* Ska denna vara med? */
                        
                        len_t id_Ip, id_Imk2, id_Tcold, id_Ti; 

                public:
                        real_t a, B, l_MK2;
                        real_t I_ref = 100000;
                        real_t B_v   = 1e-3; 
                
                        ConefinementTime(FVM::UnknownQuantityHandler *u);

                        /* Är real_t rätt data-typ?  */
                        real_t EvaluateConfinementTime();

                        real_t EvaluateConfinementTime_dIp();

                        real_t EvaluateConfinementTime_dIMK2();

                        real_t EvaluateConfinementTime_dTe();

                        real_t EvaluateConfinementTime_dTi();
}

#endif/*_STREAM_EQUATION_CONFINEMENT_TIME_HPP */
