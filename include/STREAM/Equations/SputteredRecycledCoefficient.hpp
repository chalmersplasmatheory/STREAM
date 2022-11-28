#ifndef _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP
#define _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP

#include "DREAM/IonHandler.hpp"
#include "STREAM/stream.h"
#include "DREAM/config.h"
#include "FVM/Interpolator1D.hpp"

namespace STREAM{ 
    class SputteredRecycledCoefficient {
        private: 
            DREAM::FVM::Interpolator1D ***coefficientTable;
            real_t **currCoefficientTable;
            
            DREAM::IonHandler* ions;
      
        public:
            SputteredRecycledCoefficient(DREAM::FVM::Interpolator1D***, DREAM::IonHandler*);
            ~SputteredRecycledCoefficient();
            
            void AllocateForCurrTable();

            real_t GetSRCoefficient(len_t upper, len_t lower);
            
            bool NeedsRebuild(const real_t);
            bool Rebuild(const real_t);
    };
}

#endif/*_STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP*/
