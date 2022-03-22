#ifndef _STREAM_EQUATIONS_OPTICAL_THICKNESS_HPP
#define _STREAM_EQUATIONS_OPTICAL_THICKNESS_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <gsl/gsl_math.h>
#include "DREAM/Constants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"


namespace STREAM{
        class OpticalThickness {

                protected:
                        DREAM::FVM::UnknownQuantityHandler *unknowns; 
                        
                        EllipticalRadialGridGenerator *radials;
                        
                        len_t id_Tcold, id_ncold;
                        
                        // Input parameters
		        len_t N;
		        real_t theta, phi; 

                public:
                        OpticalThickness(DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, len_t, real_t, real_t);
                        
                        real_t EvaluateOpticalThickness_o(len_t ir);
                        real_t EvaluateOpticalThickness_x(len_t ir);
                        
                        real_t EvaluateOpticalThickness_o_dTe(len_t ir);
                        real_t EvaluateOpticalThickness_x_dTe(len_t ir);
                        
                        real_t EvaluateOpticalThickness_o_dne(len_t ir);
                        real_t EvaluateOpticalThickness_x_dne(len_t ir);
                        
                        //void Initialize();
    };
}

#endif/*_STREAM_EQUATIONS_OPTICAL_THICKNESS_HPP */
