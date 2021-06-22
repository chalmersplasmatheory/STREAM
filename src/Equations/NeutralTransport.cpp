#include "STREAM/Equations/NetralTransport.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
NeutralTransport::NeutralTransport(FVM::Grid *g, IonHandler *ihdl,
	const len_t iIon, NeutralInflux *NI, PlasmaVolume *PV, 
	real_t vessel_vol) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon), 
	NI(NI), PV(PV), vessel_vol(vessel_vol) {
    SetName("NeutralTransport");
}

void NeutralTransport::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*){
    real_t Gamma0 = this->NI->EvaluateNeutralInflux(t, iIon);
    real_t V_p = PV->GetPlasmaVolume();
    real_t V_ni= PV->GetNeutralVolume(iIon);
    real_t dGamma0dnij   = this->NI->EvaluateNeutralInflux_dnij(t, iIon);
    real_t dGamma0dIp    = this->NI->EvaluateNeutralInflux_dIp(t, iIon); 
    real_t dGamma0dIwall = this->NI->EvaluateNeutralInflux_dIwall(t, iIon); 
    real_t dGamma0dTcold = this->NI->EvaluateNeutralInflux_dTcold(t, iIon); 
    real_t dGamma0dWi    = this->NI->EvaluateNeutralInflux_dWi(t, iIon); 
    real_t dGamma0dNi    = this->NI->EvaluateNeutralInflux_dNi(t, iIon);
        
    this->wall_term = Gamma0/(vessel_vol-V_p+V_ni);
    this->dn_ij     = dGamma0dni/(vessel_vol-V_p+V_nij);
    this->dI_p      = dGamma0dIp/(vessel_vol-V_p+V_ni);
    this->dI_wall   = dGamma0dIwall/(vessel_vol-V_p+V_ni);
    this->dT_cold   = dGamma0dTcold/(vessel_vol-V_p+V_ni);
    this->dW_i      = dGamma0dWi/(vessel_vol-V_p+V_ni);
    this->dN_i      = dGamma0dNi/(vessel_vol-V_p+V_ni);
}

bool IonTransport::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t, const len_t rOffset
) {    
    if(derivId==uqtyId){
        nZ = ions->GetNZ();
        for (len_t k = 0, idx = 0; k < nZ; k++) {
            len_t Z0 = ions->GetZ0(k);
            for (len_t l = 0; l <= Z0; l++, idx++) {
                jac->SetElement(rOffset+idx, rOffset+idx,this->dn_ij);
            }
        }  
		return true;
    } else 
    if(derivId==id_Ip){
		jac->SetElement(rOffset, 0,this->dI_p);
		return true;
	} else if(derivId==id_Iwall){
		jac->SetElement(rOffset, 0,this->dI_wall);
		return true;
	} else if(derivId==id_Tcold){
		jac->SetElement(rOffset, 0,this->dT_cold);
		return true;
	} else if(derivId==id_Wi){
		jac->SetElement(rOffset, iIon,this->dW_i);
		return true;
	} else if(derivId==id_Ni){
		jac->SetElement(rOffset, iIon,this->dN_i);
		return true;
	}
	else {
	    return false;
    }
}
void IonTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t, const len_t, const len_t rOffset
) {
    mat->SetElement(rOffset, rOffset, wall_term); // Är det en operator som ska verka på en parameter? Vilken parameter?
} 


void IonTransport::SetCSVectorElements(
    real_t* vec, const real_t*, const len_t iIon, const len_t, const len_t rOffset
) {
    vec[rOffset]=wall_term; 
}
