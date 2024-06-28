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
#include <array>
#include <functional>
#include "STREAM/Equations/Coefficients.hpp"

namespace STREAM
{
    class EquationSystem; // To avoid to include the file -> circular definitions
    
    class ConfinementTime 
    {

        protected:
                DREAM::FVM::UnknownQuantityHandler *unknowns; 
                
                EllipticalRadialGridGenerator *radials;
                DREAM::IonHandler *ions;
                
                len_t id_Tcold, id_Wi, id_Ni, id_Ip, id_Ncold, id_Efield; 
                
                ConnectionLength *CL;
                
                len_t D_index;
                
                static constexpr real_t ConversionFactor = 1e-6;
                static constexpr real_t nullset = 0.0;
                static constexpr real_t Conv_n20 = 1e-20;
                
                bool mixedConfLaw;
                
                EquationSystem *eqsys;
                
                OptionConstants::Conf_Time_type type;
                
                real_t GetOhmicPower(len_t ir) const; // Computes the ohmic power
                
                real_t Bohm_ConfinementTime(len_t ir) const; // Computes the confinement time from Bohm
                                        
                real_t INTOR_ConfinementTime(len_t ir) const; // Computes INTOR scaling law
                
                real_t ITER89_OL_ConfinementTime(len_t ir) const; // Computes ITER89 offset linear scaling law
                
                real_t OS_OL_ConfinementTime(len_t ir) const; // Computes Odajima-Shimomura offset linear scaling law
                
                real_t RL_OL_ConfinementTime(len_t ir) const; // Computes Rebut-Lallia offset linear scaling law

                real_t RLW_ConfinementTime(len_t ir) const; // Computes Rebut-Lallia-Watkins scaling law
                
                real_t Goldstone_scaling(PowerList const&, len_t ir, real_t const& ConvUnit = 1e-19) const;
                
                const PowerList GetCoeff() const;
                const PowerList GetCoeff(OptionConstants::Conf_Time_type _type) const;
                
                const real_t EvaluatePerpendicularConfinementTimeType(len_t ir) const;
                
        public:
                ConfinementTime(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
            DREAM::IonHandler*, ConnectionLength*, EquationSystem *eqsys, len_t D_index=0);
                
                real_t KappaOut();
                
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
