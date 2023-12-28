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
#include "STREAM/Equations/Coefficients.hpp"
//#include "STREAM/Settings/OptionConstants.hpp"
//#include "STREAM/EquationSystem.hpp" // MINE
//#include "STREAM/Equations/PlasmaVolume.hpp"
//#include "DREAM/EquationSystem.hpp"

/** 
 * Implementation of the particle confinement time
 */
/*namespace DREAM{
	class EquationSystem;
}*/
//enum STREAM::OptionConstants::Conft_Time_type;

namespace STREAM{
		class EquationSystem;
		//class OptionConstants::Conf
        class ConfinementTime {

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
                        
                        // BEGINNING OF MY LINES
                        EquationSystem *eqsys;
                        
                        OptionConstants::Conf_Time_type type;
                        
                        real_t GetOhmicPower(len_t ir);
                        
                        real_t Bohm_ConfinementTime(len_t ir);
                                                
						real_t INTOR_ConfinementTime(len_t ir);
						
						real_t ITER89_OL_ConfinementTime(len_t ir);
						
						real_t OS_OL_ConfinementTime(len_t ir);
						
						real_t RL_OL_ConfinementTime(len_t ir);
						
						real_t Goldstone_scaling(PowerList const&, len_t ir, real_t const& ConvUnit = 1e-19);
						
						const PowerList GetCoeff();
                        
                public:
                        ConfinementTime(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
					DREAM::IonHandler*, ConnectionLength*, EquationSystem *eqsys, len_t D_index=0); //MODIFIED
                        
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
