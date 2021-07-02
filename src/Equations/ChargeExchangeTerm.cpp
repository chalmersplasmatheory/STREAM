#include "STREAM/Equations/ChargeExchangeTerm.hpp"

using namespace STREAM;
using namespace DREAM;

/**
 * Constructor
 */
ChargeExchangeTerm::ChargeExchangeTerm(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ions, const len_t iIon, ADAS *adas, PlasmaVolume *pv, FVM::Grid *operandGrid) 
    : DiagonalComplexTerm(g, u, operandGrid), ions(ions), iIon(iIon), adas(adas), pv(pv) {
    
    this->DiagonalTerm::SetName("ChargeExchangeTerm");
    
    this->id_Tcold = u->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD); 
    this->id_ncold = u->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    this->id_Wi    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    this->id_ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES);
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void ChargeExchangeTerm::SetMatrixElements(FVM::Matrix *mat, real_t*) { 
    len_t n = ions->GetIndex(iIon,0);
    mat->SetElement(iIon, n, weights[n]); 
}

/**
 * Set function vector for this term.
 */
void ChargeExchangeTerm::SetVectorElements(real_t *vec, const real_t *x) { 
    len_t n = ions->GetIndex(iIon,0);
    vec[iIon] += weights[n] * x[n]; 
}

/**
 * Set of weights for this diagonal term
 */
void ChargeExchangeTerm::SetWeights(){ 
    real_t V_p  = pv->GetPlasmaVolume();
    real_t V_ni = pv->GetNeutralVolume(iIon);
    
    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[0];
    real_t n_cold = unknowns->GetUnknownData(id_ncold)[0];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[0];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[0];    
        
    len_t nZ = ions->GetNZ();
    len_t n = ions->GetIndex(iIon,0); 
    
    real_t R_icx;
    
    for(len_t iz=0; iz<nZ; iz++) {
        real_t n_i = ions->GetIonDensity(0, iz, 1);
        if(ions->IsTritium(iz)){
            R_icx = adas->GetCCD(1,3)->Eval(1, n_cold, T_cold);
        } else { 
            len_t Z  = ions->GetZ(iz);
            R_icx = adas->GetCCD(Z)->Eval(1, n_cold, T_cold);
        }
        weights[n] += V_ni/V_p * 3/2 * (2/3 * W_i / N_i - T_0) * R_icx * n_i; 
    }
}

/**
 * Set of derivatives of weights for this diagonal term
 */
void ChargeExchangeTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
    real_t V_p  = pv->GetPlasmaVolume();
    real_t V_ni = pv->GetNeutralVolume(iIon);
    real_t dV_nidT = pv->GetNeutralVolume_dT(iIon);
    real_t dV_nidn = pv->GetNeutralVolume_dn(iIon);
    
    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[0];
    real_t n_cold = unknowns->GetUnknownData(id_ncold)[0];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[0];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[0];
    
        
    len_t nZ = ions->GetNZ();
    len_t n = ions->GetIndex(iIon,0); 
    
    real_t R_icx;
    
    if(derivId == id_Tcold) {
        for(len_t iz=0; iz<nZ; iz++) {
            real_t dR_icxdT;
            real_t n_i = ions->GetIonDensity(0, iz, 1);
            if(ions->IsTritium(iz)){
                R_icx = adas->GetCCD(1,3)->Eval(1, n_cold, T_cold);
                dR_icxdT = adas->GetCCD(1,3)->Eval_deriv_T(1, n_cold, T_cold);
            } else {
                len_t Z = ions->GetZ(iz);
                R_icx = adas->GetCCD(Z)->Eval(1, n_cold, T_cold);
                dR_icxdT = adas->GetCCD(Z)->Eval_deriv_T(1, n_cold, T_cold);
            }
            diffWeights[iIon] += dV_nidT/V_p * 3/2 * (2/3 * W_i / N_i - T_0) * R_icx * n_i + V_ni/V_p * 3/2 * (2/3 * W_i / N_i-T_0) * dR_icxdT * n_i; 
        }
    } else if(derivId == id_ncold) {
        for(len_t iz=0; iz<nZ; iz++) {
            real_t dR_icxdn;
            real_t n_i = ions->GetIonDensity(0, iz, 1);
            if(ions->IsTritium(iz)){
                R_icx = adas->GetCCD(1,3)->Eval(1, n_cold, T_cold);
                dR_icxdn = adas->GetCCD(1,3)->Eval_deriv_n(1, n_cold, T_cold);
            } else {
                len_t Z = ions->GetZ(iz);
                R_icx = adas->GetCCD(Z)->Eval(1, n_cold, T_cold);
                dR_icxdn = adas->GetCCD(Z)->Eval_deriv_n(1, n_cold, T_cold);
            }
            diffWeights[iIon] += dV_nidn/V_p * 3/2 * (2/3 * W_i / N_i - T_0) * R_icx * n_i + V_ni/V_p * 3/2 * (2/3 * W_i / N_i-T_0) * dR_icxdn * n_i; 
        }
    } else if(derivId == id_Wi) {
        for(len_t iz=0; iz<nZ; iz++) {
            real_t n_i = ions->GetIonDensity(0, iz, 1);
            if(ions->IsTritium(iz)){ 
                R_icx = adas->GetCCD(1,3)->Eval(1, n_cold, T_cold);
            } else {
                len_t Z = ions->GetZ(iz);
                R_icx = adas->GetCCD(Z)->Eval(1, n_cold, T_cold);
            }
            diffWeights[iIon*nZ+n] += V_ni/V_p * 3/2 * (2/3 * 1 / N_i) * R_icx * n_i; 
        }
    } else if(derivId == id_Ni) {
        for(len_t iz=0; iz<nZ; iz++) {
            real_t n_i = ions->GetIonDensity(0, iz, 1);
            if(ions->IsTritium(iz)){
                R_icx = adas->GetCCD(1,3)->Eval(1, n_cold, T_cold);
            } else {
                len_t Z = ions->GetZ(iz);
                R_icx = adas->GetCCD(Z)->Eval(1, n_cold, T_cold);
            }
            diffWeights[iIon*nZ+n] += V_ni/V_p * 3/2 * (-2/3 * W_i / (N_i*N_i)) * R_icx * n_i; 
        }
    } else if(derivId == id_ni) {
        for(len_t iz=0; iz<nZ; iz++) {
            len_t n_iz = ions->GetIndex(iz,1);
            if(ions->IsTritium(iz)){
                R_icx = adas->GetCCD(1,3)->Eval(1, n_cold, T_cold);
            } else {
                len_t Z = ions->GetZ(iz);
                R_icx = adas->GetCCD(Z)->Eval(1, n_cold, T_cold);
            }
            diffWeights[iIon*nMultiples+n_iz] = V_ni/V_p * 3/2 * (2/3 * W_i / N_i - T_0) * R_icx; 
        }
    }
}


/**
 * Set all diffweights to 0.
 */
void ChargeExchangeTerm::ResetDiffWeights(){
    len_t nMultiples = ions->GetNzs();
    len_t nZ = ions->GetNZ();

    for(len_t i = 0; i<nMultiples*nZ; i++){
        diffWeights[i] = 0;
    }
}

/**
 * Allocate differentiation coefficients.
 */
void ChargeExchangeTerm::AllocateDiffWeights() {
    DeallocateDiffWeights();
    len_t nMultiples = ions->GetNzs();
    len_t nZ = ions->GetNZ();

    diffWeights = new real_t[nMultiples*nZ];
    ResetDiffWeights();
}

/**
 * Deallocate differentiation coefficients.
 */
void ChargeExchangeTerm::DeallocateDiffWeights() {
    if(diffWeights != nullptr)
        delete [] diffWeights;
}


