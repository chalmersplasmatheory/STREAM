#include "DREAM/Constants.hpp"
#include "STREAM/Equations/IonTransport.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"
 
using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
IonTransport::IonTransport(FVM::Grid *g, IonHandler *ihdl,
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u, PlasmaVolume *PV
	) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon), coefftauinv(tauinv), ions(ihdl), PV(PV) {
	
    SetName("IonTransport");

    this->unknowns = u;
    this->id_Ip    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    
}

/**
 * Destructor.

IonTransport::~IonTransport() {
    delete this->tauinv;
    delete this->dI_p;
    delete this->dI_wall;
    delete this->dT_cold;
    delete this->dW_i;
    delete this->dN_i;
    delete this->dn_i;
}
 */

void IonTransport::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler* 
) {
    this->id_Iwall = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_WALL);
    real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(0); 
    real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    real_t dtauinvdNi    = this->coefftauinv->EvaluateConfinementTime_dNi(0);
    
    this->tauinv  = coefftauinv->EvaluateConfinementTime(0);
    this->dn_i    = - tauinv; 
    this->dI_p    = - dtauinvdIp;
    this->dI_wall = - dtauinvdIwall;
    this->dT_cold = - dtauinvdTcold;
    this->dW_i    = - dtauinvdWi;
    this->dN_i    = - dtauinvdNi;
}

bool IonTransport::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
	if(Z0==0)
		return false;
    real_t n_i = ions->GetIonDensity(0, iIon, Z0);
    
    if(derivId==uqtyId){
		jac->SetElement(rOffset, rOffset,this->dn_i); 
		return true;
    } else if(derivId==id_Ip){
		jac->SetElement(rOffset, 0,this->dI_p*n_i);
		return true;
	} else if(derivId==id_Iwall){
		jac->SetElement(rOffset, 0,this->dI_wall*n_i);
		return true;
	} else if(derivId==id_Tcold){
		jac->SetElement(rOffset, 0,this->dT_cold*n_i);
		return true;
	} else if(derivId==id_Wi){
		jac->SetElement(rOffset, iIon,this->dW_i*n_i);
		return true;
	} else if(derivId==id_Ni){
		jac->SetElement(rOffset, iIon,this->dN_i*n_i);
		return true;
	}
	else {
	    return false;
    }
}
void IonTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t, const len_t Z0, const len_t rOffset
) {
	if(Z0==0)
		return;
    mat->SetElement(rOffset, rOffset, -tauinv);
} 


void IonTransport::SetCSVectorElements(
    real_t* vec, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
) {
	if(Z0==0)
		return;
    real_t n_i      = ions->GetIonDensity(0, iIon, Z0);
    vec[rOffset] += -n_i*tauinv;
}

