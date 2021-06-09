#ifndef _STREAM_EQUATION_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATION_CONFINEMENT_TIME_HPP

/* #include "FVM/Equation/EvaluableEquationTerm.hpp" */
#include "FVM/UnknownQuantityHandler.hpp"
/* Behövs dessa? Är adresserna rätt?  */
#include <gsl/gsl_math-h>
/* Är detta rätt matte-bibliotek? */

/** 
 * Implementation of the particle confinement time
 */
/* Är det rätt namespace? */
namespace STREAM{
        class ConefinementTime /*: public FVM::EvaluableEquationTerm /* Ska denna vara med? */ */{

                protected:
                        FVM::UnknownQuantityHandler *unknowns; /* Ska denna vara med? */
                        
                        len_t id_Ip, id_Imk2, id_Tcold, id_Ti; 

                public:
                        /* Är detta korrekt eller ska det göras på något annat sät? */
                        real_t a, B, l_MK2;
                        real_t I_ref = 100000;
                        real_t B_v   = 1e-3;
                        
                        double k_B = 0.00008617342;
                
                        ConefinementTime(FVM::UnknownQuantityHandler *u, real_t a, real_t B, real_t l_MK2);

                        /* Är real_t rätt data-typ?  */
                        real_t EvaluateConfinementTime(len_t ir);

                        real_t EvaluateConfinementTime_dIp(len_t ir);

                        real_t EvaluateConfinementTime_dIMK2(len_t ir);

                        real_t EvaluateConfinementTime_dTe(len_t ir);

                        real_t EvaluateConfinementTime_dTi(len_t ir);
}

#endif/*_STREAM_EQUATION_CONFINEMENT_TIME_HPP */
