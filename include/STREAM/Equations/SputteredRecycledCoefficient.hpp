#ifndef _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP
#define _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP

#include "STREAM/stream.h"
#include "DREAM/config.h"
#include <unordered_map>

namespace STREAM{
    class SputteredRecycledCoefficient {
        private: 
            std::unordered_map<len_t /*iIon1*/, std::unordered_map<len_t /*iIon2*/, real_t /*Y*/> > *coefficientTable;
      
        public:
            SputteredRecycledCoefficient(std::unordered_map<len_t /*iIon1*/, std::unordered_map<len_t /*iIon2*/, real_t /*Y*/> > *coefficientTable);
            void AddSRCoefficient(len_t upper, len_t lower, real_t coefficient);
            real_t GetSRCoefficient(len_t upper, len_t lower);
    };
}
#endif/*_STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP*/
