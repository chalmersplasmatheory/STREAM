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
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u
	) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon, allocCoefficients), coefftauinv(tauinv), ions(ihdl) {
	
    SetName("IonHeatTransport");

    this->unknowns = u;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
	Allocate();
}

/**
 * Destructor.
 */
IonHeatTransport::~IonHeatTransport() {
    delete this->tauinv;
    delete this->dI_p;
    delete this->dI_wall;
    delete this->dT_cold;
    delete this->dW_i;
    delete this->dN_i;
    delete this->dn_i;
}

void IonHeatTransport::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler* 
) {
    
    real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTe(0); 
    real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    real_t dtauinvdNi    = this->coefftauinv->EvaluateConfinementTime_dni(0);
    
    this->tauinv  = coefftauinv->EvaluateConfinementTime(0);
    this->dn_i    = - 3/2 * Constants::ec * tauinv; 
    this->dI_p    = - 3/2 * Constants::ec * dtauinvdIp;
    this->dI_wall = - 3/2 * Constants::ec * dtauinvdIwall;
    this->dT_cold = - 3/2 * Constants::ec * dtauinvdTcold;
    this->dW_i    = - 3/2 * Constants::ec * dtauinvdWi;
    this->dN_i    = - 3/2 * Constants::ec * dtauinvdNi;
    
    }
}

bool IonHeatTransport::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    real_t n_i = ions->GetIonDensity(0, iIon, Z0);
    
    if(derivId==uqtyId){
		jac->SetElement(rOffset+Z0, rOffset+Z0,this->dn_i); 
		return true;
    } else if(derivId==id_Ip){
		jac->SetElement(rOffset+Z0, 0,this->dI_p*n_i);
		return true;
	} else if(derivId==id_Iwall){
		jac->SetElement(rOffset+Z0, 0,this->dI_wall*n_i);
		return true;
	} else if(derivId==id_Tcold){
		jac->SetElement(rOffset+Z0, 0,this->dT_cold*n_i);
		return true;
	} else if(derivId==id_Wi){
		jac->SetElement(rOffset+Z0, iIon,this->dW_i*n_i);
		return true;
	} else if(derivId==id_Ni){
		jac->SetElement(rOffset+Z0, iIon,this->dN_i*n_i);
		return true;
	}
	else {
	    return false;
    }
}
void IonHeatTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    mat->SetElement(rOffset+Z0, rOffset+Z0, -3/2 * Constants::ec*tauinv);
} 


void IonHeatTransport::SetCSVectorElements(
    real_t* vec, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    real_t n_i     = ions->GetIonDensity(0, iIon, Z0);
    vec[rOffset+Z0]=-3/2 * Constants::ec*n_i*tauinv; 
}

