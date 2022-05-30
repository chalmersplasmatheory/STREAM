#include "STREAM/Equations/NeutralTransport.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor. 
 */
NeutralTransport::NeutralTransport(FVM::Grid *g, IonHandler *ihdl,
	const len_t iIon, FVM::UnknownQuantityHandler *u, NeutralInflux *NI, PlasmaVolume *PV) : 
	IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon), NI(NI), PV(PV) {
    SetName("NeutralTransport");
    
    this->unknowns   = u;
    this->id_Ip      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    this->id_Tcold   = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_Wi      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    this->id_ncold   = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    this->id_lambdai = unknowns->GetUnknownID(STREAM::OptionConstants::UQTY_LAMBDA_I);
    
    len_t nZ = ions->GetNZ();
    for (len_t k = 0; k < nZ; k++) {
        len_t Z = ions->GetZ(k);
        for (len_t l = 1; l <= Z; l++) {
            sum_derivs++;
        }
    }
    this->dn_kj = new real_t[nZ]; //[nZ+1];
    
}

void NeutralTransport::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*){
    this->id_Iwall = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_WALL);
    real_t Gamma0        = this->NI->EvaluateNeutralInflux(iIon, t); // (t, iIon);
    real_t dGamma0dIp    = this->NI->EvaluateNeutralInflux_dIp(iIon); // (t, iIon); 
    real_t dGamma0dIwall = this->NI->EvaluateNeutralInflux_dIwall(iIon); // (t, iIon); 
    real_t dGamma0dTcold = this->NI->EvaluateNeutralInflux_dTcold(iIon); // (t, iIon); 
    real_t dGamma0dWi    = this->NI->EvaluateNeutralInflux_dWi(iIon); // (t, iIon); 
    real_t dGamma0dNi    = this->NI->EvaluateNeutralInflux_dNi(iIon); // (t, iIon);
            
    real_t V_tot         = PV->GetTotalNeutralVolume(iIon);
    real_t dVtotdlambdai = PV->GetTotalNeutralVolume_dLambdai(iIon);
        
    this->wall_term = Gamma0/V_tot;
    this->dI_p      = dGamma0dIp/V_tot;
    this->dI_wall   = dGamma0dIwall/V_tot;
    this->dT_cold   = dGamma0dTcold/V_tot;
    this->dW_i      = dGamma0dWi/V_tot;
    this->dN_i      = dGamma0dNi/V_tot;
    this->dlambda_i = - Gamma0/(V_tot*V_tot)*dVtotdlambdai;
    
    len_t nZ = ions->GetNZ();
    real_t dGamma0dnkj=0;
    for (len_t kIon = 0; kIon < nZ; kIon++) {
        dGamma0dnkj = this->NI->EvaluateNeutralInflux_dnkj(iIon, kIon); // (t, iIon, kIon);
        this->dn_kj[kIon] = dGamma0dnkj/V_tot; // this->dn_kj[kIon+1]
    }
}

bool NeutralTransport::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if (Z0 != 0)
        return false;

    if(derivId==uqtyId){
        len_t nZ = ions->GetNZ();
        for (len_t kIon = 0, idx = 0; kIon < nZ; kIon++) {
            len_t Z_k = ions->GetZ(kIon);
            idx++;
            for (len_t j = 1; j <= Z_k; j++, idx++) {
                jac->SetElement(rOffset, idx, this->dn_kj[kIon]);  // this->dn_kj[kIon+1]
            }
        }
        
		return true;
    } else if(derivId==id_Ip){
		jac->SetElement(rOffset, 0, this->dI_p);
		return true;
	} else if(derivId==id_Iwall){
		jac->SetElement(rOffset, 0, this->dI_wall);
		return true;
	} else if(derivId==id_Tcold){
		jac->SetElement(rOffset, 0, this->dT_cold);
		return true;
	} else if(derivId==id_Wi){
		jac->SetElement(rOffset, iIon, this->dW_i);
		return true;
	} else if(derivId==id_Ni){
		jac->SetElement(rOffset, iIon, this->dN_i);
		return true;
	} else if(derivId==id_lambdai){
		jac->SetElement(rOffset, iIon, this->dlambda_i);
		return true;
	}
	else {
	    return false;
    }
}
void NeutralTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t, const len_t Z0, const len_t rOffset
) {
    if (Z0 != 0)
        return;

    mat->SetElement(rOffset, rOffset, wall_term); 
} 


void NeutralTransport::SetCSVectorElements(
    real_t* vec, const real_t*, const len_t , const len_t Z0, const len_t rOffset
) {
    if (Z0 != 0)
        return;

    vec[rOffset]+=wall_term; 
}
