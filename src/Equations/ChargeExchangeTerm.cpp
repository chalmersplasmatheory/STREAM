#include "STREAM/Equations/ChargeExchangeTerm.hpp"

using namespace STREAM;
using namespace DREAM;

/**
 * Constructor
 */
ChargeExchangeTerm::ChargeExchangeTerm(FVM::Grid *g, IonHandler *ions, const len_t iIon, FVM::UnknownQuantityHandler *u, ADAS *adas, PlasmaVolume *pv) 
    : DiagonalLinearTerm(g), ions(ions), iIon(iIon), unknowns(u), adas(adas), pv(pv) {
    
    this->DiagonalTerm::SetName("ChargeExchangeTerm");
    
    this->id_Wi    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void ChargeExchangeTerm::SetMatrixElements(FVM::Matrix *mat, real_t*) { 
    len_t nMultiples = ions->GetNzs();
    len_t n = ions->GetIndex(iIon,0);
    mat->SetElement(iIon, n, weights[n]); 
}

/**
 * Set function vector for this term.
 */
void ChargeExchangeTerm::SetVectorElements(real_t *vec, const real_t *x) { 
    len_t nMultiples = ions->GetNzs();
    len_t n = ions->GetIndex(iIon,0);
    vec[iIon] += weights[n] * x[n]; 
}

/**
 * Implementation of weights for this diagonal term
 */
void ChargeExchangeTerm::SetWeights(){ 
    real_t V_p  = pv->GetPlasmaVolume();
    real_t V_ni = pv->GetNeutralVolume(iIon);
    
    real_t W_i = unknowns->GetUnknownData(id_Wi)[0];
    real_t N_i = unknowns->GetUnknownData(id_Ni)[0];
    real_t T_i = 2/3 * W_i / N_i;
        
    len_t NZ = ions->GetNZ();
    len_t n = ions->GetIndex(iIon,0); 
    
    real_t R_icx = 0;
    
    for(len_t iz=0; iz<NZ; iz++) {
        real_t n_i = ions->GetIonDensity(0, iz, 1);
        if(ions->IsTritium(iIon)){
            R_icx = adas->GetCCD(1,3);
        } else {
            R_icx = adas->GetCCD(1);
        }
        weights[n] += V_ni/V_p * 3/2 * (T_i-T_0) * R_icx * n_i; 
    }
}


