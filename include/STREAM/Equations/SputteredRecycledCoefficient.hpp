#ifndef _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP
#define _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP

namespace STREAM{
    class SputteredRecycledCoefficient {
        private: 
            std::unordered_map< len_t /*iIon1*/, std::unordered_map<len_t /*iIon2*/, real_t /*Y*/> > 
                coefficientTable = {
                    {/* D */ , {
                                {/* C */,0},
                                {/* O */,0}
                               }
                    },
                    {/* C */ , {
                                {/* D */,0.03},
                                {/* C */,0},
                                {/* O */,1}
                               }
                    },
                    {/* O */ , {
                                {/* D */,0},
                                {/* C */,0},
                                {/* O */,1}
                               }
                    } 
                };
        /*
            coefficientTable[ D ] = coefficient_D_row;
                coefficient_D_row[ C ] = 0;
                coefficient_D_row[ O ] = 0;
            coefficientTable[ C ] = coefficient_C_row;
                coefficient_D_row[ D ] = 0.03;
                coefficient_D_row[ C ] = 0;
                coefficient_D_row[ O ] = 1;
            coefficientTable[ O ] = coefficient_O_row;
                coefficient_D_row[ D ] = 0;
                coefficient_D_row[ C ] = 0;
                coefficient_D_row[ O ] = 1;
        */
        public:
            SputteredRecycledCoefficient();
            void add(len_t upper, len_t lower, real_t coefficient);
            real_t get(len_t upper, len_t lower);
    };
}
#endif/*_STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP*/
