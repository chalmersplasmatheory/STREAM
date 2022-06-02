#ifndef _STREAM_EQUATIONS_CONNECTION_LENGTH_HPP
#define _STREAM_EQUATIONS_CONNECTION_LENGTH_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <gsl/gsl_math.h>
#include "DREAM/Constants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"



/** 
 * Implementation of the connection length
 */

namespace STREAM{
        class ConnectionLength {

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; 
                        
                        EllipticalRadialGridGenerator *radials;
                        
                        len_t id_Ip, id_Iwall; 

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
			real_t connectionLengthFactor;
                
                        ConnectionLength(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
					real_t l_MK2, real_t B_v, real_t I_ref, real_t connectionLengthFactor);
                        
                        real_t EvaluateInverseConnectionLength(len_t ir);

                        real_t EvaluateInverseConnectionLength_dIp(len_t ir);

                        real_t EvaluateInverseConnectionLength_dIwall(len_t ir);

                        real_t EvaluateConnectionLength(len_t ir);

                        void Initialize();
    };
}

#endif/*_STREAM_EQUATIONS_CONNECTION_LENGTH_HPP */
