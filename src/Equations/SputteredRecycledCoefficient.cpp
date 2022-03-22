#include "STREAM/Equations/SputteredRecycledCoefficient.hpp"

using namespace STREAM;
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
    return coefficientTable[lower][upper];
}

/*

/ **
 * Constructor
 * /
SputteredRecycledCoefficient::SputteredRecycledCoefficient(
    DREAM::FVM::Interpolator1D **coefficientTable, DREAM::IonHandler *ihdl
) : coefficientTable(coefficientTable), ions(ihdl) {
    
}

/ **
 * Destructor.
 * /
SputteredRecycledCoefficient::~SputteredRecycledCoefficient() {
    delete [] coefficientTable[0];
    delete [] coefficientTable;
    
    delete [] currCoefficientTable[0];
    delete [] currCoefficientTable;
}

SputteredRecycledCoefficient::AllocateForCurrTable() {
    len_t nZ = ions->GetNz();
    currCoefficientTable = new real_t*[nZ];
    for (len_t z = 0; z < nZ; z++) {
    	currCoefficientTable[z] = new real_t*[nZ];
    }
}

/ **
 * Access the sputtering-recycling coefficient
 *
 *   Y^upper_lower
 *
 * for incident ion species 'lower' and sputtered ion species 'upper'.
 * /
real_t SputteredRecycledCoefficient::GetSRCoefficient(len_t upper, len_t lower){
    return currCoefficientTable[lower][upper];
}

/ **
 * This method indicates whether or not the coefficient table needs to
 * be rebuilt.
 *
 * t: Time at which to check if a coefficient table rebuild is needed.
 * /
bool SputteredRecycledCoefficient::NeedsRebuild(const real_t t) const {
    len_t nZ = ions->GetNz();
    real_t Y;
    bool needsRebuild = false;
    len_t lIon=0;
    while (lIon < nZ && !needsRebuild) {
    	len_t uIon=0;
    	while (uIon < nZ && !needsRebuild) {
    	    Y = *this->coefficientTable[lIon][uIon]->Eval(t);
    	    if (Y != currCoefficientTable[lIon][uIon]) {
    	    	needsRebuild = true;
    	    }
    	}
    }
    return (needsRebuild);
}

/ **
 * Rebuild this grid.
 * /
bool SputteredRecycledCoefficient::Rebuild(const real_t t) {
    len_t nZ = ions->GetNz();
    real_t Y;
    for (len_t lIon = 0; lIon < nZ; lIon++) {
    	for (len_t uIon = 0; uIon < nZ; uIon++) {
    	    currCoefficientTable[lIon][uIon] = *this->coefficientTable[lIon][uIon]->Eval(t);
    	}
    }
    return true;
}


*/
