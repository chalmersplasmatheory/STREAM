#include "DREAM/Constants.hpp"
#include "STREAM/Equations/IonHeatTransport.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
IonHeatTransport::IonHeatTransport(FVM::Grid *g, IonHandler *ihdl,
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u, 
	EllipticalRadialGridGenerator *r, len_t D_index
	) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon), coefftauinv(tauinv), ions(ihdl), radials(r), D_index(D_index) {
	
    SetName("IonHeatTransport");

    this->unknowns = u;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
}

/**
 * Destructor.

IonHeatTransport::~IonHeatTransport() {
    delete this->tauinv;
    delete this->dI_p;
    delete this->dI_wall;
    delete this->dT_cold;
    delete this->dW_i;
    delete this->dN_i;
    delete this->dn_i;
}
 */

void IonHeatTransport::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler* 
) {    
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(0); 
    real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    real_t dtauinvdNi    = this->coefftauinv->EvaluateConfinementTime_dNi(0);
    len_t nr = radials->GetNr();
    this->W_i     = unknowns->GetUnknownData(id_Wi)[D_index*nr];
    this->tauinv  = coefftauinv->EvaluateConfinementTime(0);
    this->dI_p    = - 2.0/3.0 * dtauinvdIp * W_i;
    this->dI_wall = - 2.0/3.0 * dtauinvdIwall * W_i;
    this->dT_cold = - 2.0/3.0 * dtauinvdTcold * W_i;
    this->dW_i    = - 2.0/3.0 * ( dtauinvdWi * W_i + tauinv);
    this->dN_i    = - 2.0/3.0 * dtauinvdNi * W_i;
}

bool IonHeatTransport::SetCSJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t Z0, const len_t rOffset
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
void IonHeatTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t, const len_t Z0, const len_t rOffset
) {
    mat->SetElement(rOffset, rOffset, -2.0/3.0 * tauinv);
} 


void IonHeatTransport::SetCSVectorElements(
    real_t* vec, const real_t*, const len_t, const len_t Z0, const len_t rOffset
) {
    vec[rOffset]-=2.0/3.0 * W_i * tauinv; 
}

