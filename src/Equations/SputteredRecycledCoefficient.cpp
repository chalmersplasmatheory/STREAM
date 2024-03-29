#include "STREAM/Equations/SputteredRecycledCoefficient.hpp"

using namespace STREAM;
using namespace std;

/**
 * Constructor
 */
SputteredRecycledCoefficient::SputteredRecycledCoefficient(
    DREAM::FVM::Interpolator1D ***coefficientTable, DREAM::IonHandler *ihdl
) : coefficientTable(coefficientTable), ions(ihdl) {
    AllocateForCurrTable();
}

/**
 * Destructor.
 */
SputteredRecycledCoefficient::~SputteredRecycledCoefficient() {
    delete [] coefficientTable[0];
    delete [] coefficientTable;
    
    delete [] currCoefficientTable[0];
    delete [] currCoefficientTable;
}

void SputteredRecycledCoefficient::AllocateForCurrTable() {
    len_t nZ = ions->GetNZ();
    currCoefficientTable = new real_t*[nZ];
    currCoefficientTable[0] = new real_t[nZ*nZ];
    for (len_t z = 1; z < nZ; z++) 
    	currCoefficientTable[z] = currCoefficientTable[z-1] + nZ;
}

/**
 * Access the sputtering-recycling coefficient
 *
 *   Y^upper_lower
 *
 * for incident ion species 'lower' and sputtered ion species 'upper'.
 */
real_t SputteredRecycledCoefficient::GetSRCoefficient(len_t upper, len_t lower){
    return currCoefficientTable[lower][upper];
}

/**
 * This method indicates whether or not the coefficient table needs to
 * be rebuilt.
 *
 * t: Time at which to check if a coefficient table rebuild is needed.
 */
bool SputteredRecycledCoefficient::NeedsRebuild(const real_t t) {
    len_t nZ = ions->GetNZ();
    real_t Y;
    bool needsRebuild = false;
    len_t lIon=0;
    while (lIon < nZ && !needsRebuild) {
    	len_t uIon=0;
    	while (uIon < nZ && !needsRebuild) {
    	    Y = *this->coefficientTable[lIon][uIon]->Eval(t);
    	    if (Y != currCoefficientTable[lIon][uIon]) {
    	    	return true;
    	    }
    	}
    }
    return false;
}

/**
 * Rebuild this grid.
 */
bool SputteredRecycledCoefficient::Rebuild(const real_t t) {
    len_t nZ = ions->GetNZ();
    for (len_t lIon = 0; lIon < nZ; lIon++) {
    	for (len_t uIon = 0; uIon < nZ; uIon++) {
    	    if (this->coefficientTable[lIon][uIon] == nullptr) 
    	        currCoefficientTable[lIon][uIon] = 0;
            else
        	currCoefficientTable[lIon][uIon] = *this->coefficientTable[lIon][uIon]->Eval(t);
    	}
    }
    return true;
}
