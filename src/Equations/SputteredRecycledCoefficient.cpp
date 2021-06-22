#include "STREAM/Equations/SputteredRecycledCoefficient.hpp"

using namespace STREAM;
//using namespace DREAM;
using namespace std;

/**
 * Constructor
 */
SputteredRecycledCoefficient::SputteredRecycledCoefficient(
    real_t **coefficientTable
) : coefficientTable(coefficientTable) {}

/**
 * Destructor.
 */
SputteredRecycledCoefficient::~SputteredRecycledCoefficient() {
    delete [] coefficientTable[0];
    delete [] coefficientTable;
}

/**
 * Access the sputtering-recycling coefficient
 *
 *   Y^upper_lower
 *
 * for incident ion species 'lower' and sputtered ion species 'upper'.
 */
real_t SputteredRecycledCoefficient::GetSRCoefficient(len_t upper, len_t lower){
    /*if (coefficientTable->find(upper) != coefficientTable->end() && coefficientTable->find(upper)->second.find(lower) != coefficientTable->find(upper)->second.end()){ 
        return coefficientTable->find(upper)->second.find(lower)->second;
    } else {
        return 0; 
    }*/
    return coefficientTable[lower][upper];
}
