#ifndef _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP
#define _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP

#include "STREAM/stream.h"
#include "DREAM/config.h"

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
#endif/*_STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP*/
