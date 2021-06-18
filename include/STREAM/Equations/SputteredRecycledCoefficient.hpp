#ifndef _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP
#define _STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP

namespace STREAM{
    class SputteredRecycledCoefficient {
        private: 
            std::unordered_map< len_t /*iIon1*/, std::unordered_map<len_t /*iIon2*/, real_t /*Y*/> > 
                coefficientTable = { //Är detta korekt syntax för unordered_map?
                    {/* D */ 1, { //Är detta ett rimligt sätt att göra unordered_map på? Dvs kalla D, C, O för 1, 6, 8? 
                                {/* C */ 6, 0},
                                {/* O */ 8, 0}
                               }
                    },
                    {/* C */ 6, {
                                {/* D */ 1, 0.03},
                                {/* C */ 6, 0},
                                {/* O */ 8, 1}
                               }
                    },
                    {/* O */ 8, {
                                {/* D */ 1, 0},
                                {/* C */ 6, 0},
                                {/* O */ 8, 1}
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
            SputteredRecycledCoefficient(); // Denna ska vara tom va? Det enda vi vill ha är tabellen
            void AddSRCoefficient(len_t upper, len_t lower, real_t coefficient);
            real_t GetSRCoefficient(len_t upper, len_t lower);
    };
}
#endif/*_STREAM_EQUATIONS_SPUTTERED_RECYCLED_COEFFICIENT_HPP*/
