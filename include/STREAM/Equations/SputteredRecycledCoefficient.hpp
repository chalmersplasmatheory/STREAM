#ifndef _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP
#define _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP

#include "STREAM/stream.h"
#include "DREAM/config.h"
#include "FVM/Interpolator1D.hpp"

namespace STREAM{ 
    class SputteredRecycledCoefficient {
        private: 
            real_t **coefficientTable;
      
        public:
            SputteredRecycledCoefficient(real_t**);
            ~SputteredRecycledCoefficient();

            real_t GetSRCoefficient(len_t upper, len_t lower);
    };
}

/*

namespace STREAM{ 
    class SputteredRecycledCoefficient {
        private: 
            DREAM::FVM::Interpolator1D **coefficientTable;
            real_t **currCoefficientTable;
      
        public:
            SputteredRecycledCoefficient(DREAM::FVM::Interpolator1D**);
            ~SputteredRecycledCoefficient();
            
            AllocateForCurrTable();

            real_t GetSRCoefficient(len_t upper, len_t lower);
            
            NeedsRebuild(const real_t);
            Rebuild(const real_t);
    };
}


*/

#endif/*_STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP*/
