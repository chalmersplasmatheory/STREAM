#include "STREAM/Equations/NetralTransport.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
NeutralTransport::NeutralTransport(FVM::Grid *g, IonHandler *ihdl,
	const len_t iIon, NeutralInflux *NI, PlasmaVolume *PV, ConfinementTime *CT, 
	real_t vessel_vol) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon), 
	NI(NI), PV(PV), CT(CT), vessel_vol(vessel_vol) {
    SetName("NeutralTransport");
}

void NeutralTransport::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*){
    real_t Gamma0 = this->NI->EvaluateNeutralInflux(t, iIon);
    real_t V_p = PV->GetPlasmaVolume();
    real_t V_ni= PV->GetNeutralVolume(iIon);
    real_t tauinv = coefftauinv->CT(0);
    real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(0); 
    real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    real_t dtauinvdNi    = this->coefftauinv->EvaluateConfinementTime_dNi(0);
        
    this->wall_term = Gamma0/(vessel_vol-V_p+V_ni);
    this->dI_p      = - wall_term/tauinv * dtauinvdIp;
    this->dI_wall   = - wall_term/tauinv * dtauinvdIwall;
    this->dT_cold   = - wall_term/tauinv * dtauinvdTcold;
    this->dW_i      = - wall_term/tauinv * dtauinvdWi;
    this->dN_i      = - wall_term/tauinv * dtauinvdNi;
}

bool IonTransport::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t, const len_t rOffset
) {    
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
