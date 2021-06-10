#ifndef _STREAM_EQUATION_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATION_CONFINEMENT_TIME_HPP

/* #include "FVM/Equation/EvaluableEquationTerm.hpp" */
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
/* Behövs dessa? Är adresserna rätt?  */
#include <gsl/gsl_math.h>
/* Är detta rätt matte-bibliotek? */
#include "DREAM/Constants.hpp"
/* Behövs denna? */


/** 
 * Implementation of the particle confinement time
 */
 
/* Är det rätt namespace? */
namespace STREAM{
        class ConfinementTime /*: public FVM::EvaluableEquationTerm /* Ska denna vara med? */ */{

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; /* Ska denna vara med? */
                        
                        EllipticalRadialGridGenerator *radials;
                        
                        len_t id_Ip, id_Imk2, id_Tcold, id_Wi, id_Ni; 

                public:
                        /* Är detta korrekt eller ska det göras på något annat sätt? */
                        real_t l_MK2;
                        real_t I_ref = 100000;
                        real_t B_v   = 1e-3;
                
                        ConfinementTime(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *aB, real_t l_MK2);

                        /* Är real_t rätt data-typ?  */
                        real_t EvaluateConfinementTime(len_t ir);

                        real_t EvaluateConfinementTime_dIp(len_t ir);

                        real_t EvaluateConfinementTime_dIMK2(len_t ir);

                        real_t EvaluateConfinementTime_dTe(len_t ir);

                        real_t EvaluateConfinementTime_dTi(len_t ir);
}

#endif/*_STREAM_EQUATION_CONFINEMENT_TIME_HPP */
