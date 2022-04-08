#ifndef _STREAM_EQUATIONS_CONFINEMENT_TIME_HPP
#define _STREAM_EQUATIONS_CONFINEMENT_TIME_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <gsl/gsl_math.h>
#include "DREAM/Constants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"



/** 
 * Implementation of the particle confinement time
 */

namespace STREAM{
        class ConfinementTime {

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; 
                        
                        EllipticalRadialGridGenerator *radials;
						DREAM::IonHandler *ions;
                        
                        len_t id_Ip, id_Iwall, id_Tcold, id_Wi, id_Ni; 
                        
                        len_t D_index;

                public:
                        real_t l_MK2;
                        real_t I_ref;
                        real_t B_v;

			// The following factor appears in the connection length. In "older" studies
			// (pre-Mineev 2014 (http://www-naweb.iaea.org/napc/physics/FEC/FEC2014/fec2014-preprints/255_PPCP320.pdf))
			// this factor was set to 1. But Mineev found, using 3D simulations of
			// the poloidal stray field, that it should rather be 3. To allow us to
			// easily switch between simulation modes, we put this constant here and
			// use it in the *.cpp file.
			real_t connectionLengthFactor = 1.0;
                
                        ConfinementTime(
							DREAM::FVM::UnknownQuantityHandler *u,
							EllipticalRadialGridGenerator *r,
							DREAM::IonHandler*, real_t l_MK2,
							real_t B_v, real_t I_ref, real_t connectionLengthFactor,
							len_t D_index=0
						);
                        
                        real_t EvaluateConfinementTime(len_t ir);

                        real_t EvaluateParallelConfinementTime(len_t ir);

                        real_t EvaluatePerpendicularConfinementTime(len_t ir);

                        real_t EvaluateConnectionLength(len_t ir);

                        real_t EvaluateConfinementTime_dIp(len_t ir);

                        real_t EvaluateConfinementTime_dIwall(len_t ir);

                        real_t EvaluateConfinementTime_dTcold(len_t ir);

                        real_t EvaluateConfinementTime_dWi(len_t ir);
                        
                        real_t EvaluateConfinementTime_dNi(len_t ir);
                        
                        void Initialize();
    };
}

#endif/*_STREAM_EQUATIONS_CONFINEMENT_TIME_HPP */
