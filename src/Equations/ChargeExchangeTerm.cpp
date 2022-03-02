
#include "DREAM/Constants.hpp"
#include "STREAM/Equations/ChargeExchangeTerm.hpp"

using namespace STREAM;
using namespace DREAM;

/**
 * Constructor
 */
ChargeExchangeTerm::ChargeExchangeTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ions, const len_t iIon, 
    ADAS *adas, PlasmaVolume *pv, EllipticalRadialGridGenerator *r, FVM::Grid *operandGrid,  
    len_t D_index
    ) : DiagonalComplexTerm(g, u, operandGrid), ions(ions), iIon(iIon), adas(adas), pv(pv), radials(r), D_index(D_index) {
    
    this->DiagonalTerm::SetName("ChargeExchangeTerm");
    
    this->id_Wi      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    this->id_ni      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES);
    this->id_lambdai = unknowns->GetUnknownID(STREAM::OptionConstants::UQTY_LAMBDA_I);
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
    
    len_t nr = radials->GetNr();
    len_t n = ions->GetIndex(iIon,0); 

    // Reset weights
    weights[n] = 0;
    
    real_t W_i, N_i, T_i, n_i, R_icx;
    
    W_i = unknowns->GetUnknownData(id_Wi)[iIon*nr]; 
    N_i = unknowns->GetUnknownData(id_Ni)[iIon*nr];
    if(N_i == 0) {
        T_i=0;
    } else {
        T_i=2.0/3.0 * W_i / N_i;
    }
    n_i = ions->GetIonDensity(0, iIon, 1);
    if(ions->IsTritium(iIon)){
        R_icx = adas->GetCCD(1,3)->Eval(0, n_i, T_i /DREAM::Constants::ec);
    } else if (ions->IsHydrogen(iIon)) {
        R_icx = adas->GetCCD(1,1)->Eval(0, n_i, T_i /DREAM::Constants::ec);
    } else { 
        len_t Z  = ions->GetZ(iIon);
        R_icx = adas->GetCCD(Z)->Eval(0, n_i, T_i /DREAM::Constants::ec);
    }
    weights[n] -= V_ni/V_p * 3.0/2.0 * (T_i - DREAM::Constants::ec*T_0) * R_icx * n_i;
}

/**
 * Set of derivatives of weights for this diagonal term
 */
void ChargeExchangeTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
    real_t V_p  = pv->GetPlasmaVolume();
    real_t V_ni = pv->GetNeutralVolume(iIon);
    real_t dV_nidlambdai = pv->GetNeutralVolume_dLambdai(iIon);
    
    len_t nr = radials->GetNr();
    
    ResetDiffWeights();
        
    len_t nZ = ions->GetNZ();
    len_t n = ions->GetIndex(iIon,0); 
    
    real_t R_icx, dR_icxdT, dR_icxdn;
    real_t W_i, N_i, T_i, n_i;

    W_i = unknowns->GetUnknownData(id_Wi)[iIon*nr]; 
    N_i = unknowns->GetUnknownData(id_Ni)[iIon*nr];

    if (N_i == 0)
        return;
    
    if(derivId == id_Wi) {
        T_i=2.0/3.0 * W_i / N_i;
	n_i = ions->GetIonDensity(0, iIon, 1);
	if(ions->IsTritium(iIon)){ 
	    R_icx = adas->GetCCD(1,3)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdT = adas->GetCCD(1,3)->Eval_deriv_T(0, n_i, T_i /DREAM::Constants::ec);
	} else if(ions->IsHydrogen(iIon)){ 
	    R_icx = adas->GetCCD(1,1)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdT = adas->GetCCD(1,1)->Eval_deriv_T(0, n_i, T_i /DREAM::Constants::ec);
	} else {
	    len_t Z = ions->GetZ(iIon);
	    R_icx = adas->GetCCD(Z)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdT = adas->GetCCD(Z)->Eval_deriv_T(0, n_i, T_i /DREAM::Constants::ec);
	}
	diffWeights[iIon*nZ+n] -= V_ni/V_p * 3.0/2.0 * (2.0/3.0 * 1 / N_i) * R_icx * n_i 
                                       + V_ni/V_p * 3.0/2.0 * (T_i - DREAM::Constants::ec*T_0) * dR_icxdT * 2.0/3.0 * 1 / N_i * n_i;
    } else if(derivId == id_Ni) {
	T_i=2.0/3.0 * W_i / N_i;
	n_i = ions->GetIonDensity(0, iIon, 1);
	if(ions->IsTritium(iIon)){
	    R_icx = adas->GetCCD(1,3)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdT = adas->GetCCD(1,3)->Eval_deriv_T(0, n_i, T_i /DREAM::Constants::ec);
	} else if(ions->IsHydrogen(iIon)){
	    R_icx = adas->GetCCD(1,1)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdT = adas->GetCCD(1,1)->Eval_deriv_T(0, n_i, T_i /DREAM::Constants::ec);
	} else {
	    len_t Z = ions->GetZ(iIon);
	    R_icx = adas->GetCCD(Z)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdT = adas->GetCCD(Z)->Eval_deriv_T(0, n_i, T_i /DREAM::Constants::ec);
	}
	diffWeights[iIon*nZ+n] -= V_ni/V_p * 3.0/2.0 * (-2.0/3.0 * W_i / (N_i*N_i)) * R_icx * n_i
							+V_ni/V_p * 3.0/2.0 * (T_i - DREAM::Constants::ec*T_0) * (-2.0/3.0 * W_i / (N_i*N_i) ) * dR_icxdT * n_i; 
    } else if(derivId == id_lambdai) {
	T_i=2.0/3.0 * W_i / N_i;
	n_i = ions->GetIonDensity(0, iIon, 1);
	if(ions->IsTritium(iIon)){
	    R_icx = adas->GetCCD(1,3)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	} else if(ions->IsHydrogen(iIon)){
	    R_icx = adas->GetCCD(1,1)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	} else {
	    len_t Z = ions->GetZ(iIon);
	    R_icx = adas->GetCCD(Z)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	}
	diffWeights[iIon*nZ+n] -= dV_nidlambdai/V_p * 3.0/2.0 * (T_i - DREAM::Constants::ec*T_0) * R_icx * n_i; 
    } else if(derivId == id_ni) {
	T_i=2.0/3.0 * W_i / N_i;
	len_t n_iz = ions->GetIndex(iIon,1);
	n_i = ions->GetIonDensity(0, iIon, 1);
	if(ions->IsTritium(iIon)){
	    R_icx = adas->GetCCD(1,3)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdn = adas->GetCCD(1,3)->Eval_deriv_n(0, n_i, T_i /DREAM::Constants::ec);
	} else if(ions->IsHydrogen(iIon)){
	    R_icx = adas->GetCCD(1,1)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdn = adas->GetCCD(1,1)->Eval_deriv_n(0, n_i, T_i /DREAM::Constants::ec);
	} else {
	    len_t Z = ions->GetZ(iIon);
	    R_icx = adas->GetCCD(Z)->Eval(0, n_i, T_i /DREAM::Constants::ec);
	    dR_icxdn = adas->GetCCD(Z)->Eval_deriv_n(0, n_i, T_i /DREAM::Constants::ec);
	}
	diffWeights[iIon*nMultiples+n_iz] -= V_ni/V_p * 3.0/2.0 * (T_i - DREAM::Constants::ec*T_0) * R_icx
											+V_ni/V_p * 3.0/2.0 * (T_i - DREAM::Constants::ec*T_0) * dR_icxdn * n_i;  
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


