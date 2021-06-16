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
IonTransport::IonTransport(FVM::Grid *g, IonHandler *ihdl/*, bool allocCoefficients*/, len_t iz,
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u
	) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon, allocCoefficients), coefftauinv(tauinv), ions(ihdl), iz(iz) {
	
    SetName("IonTransport");

    this->unknowns = u;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
	Allocate();
}

/**
 * Destructor.
 */
IonTransport::~IonTransport() {
    delete this->tauinv;
    delete this->dI_p;
    delete this->dI_wall;
    delete this->dT_cold;
    delete this->dW_i;
    delete this->dn_i;
    delete this->dIS;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void IonTransport::Allocate() {
    this->tauinv = new real_t*;
    this->dI_p = new real_t*;
    this->dI_wall = new real_t*;
    this->dT_cold = new real_t*;
    this->dW_i = new real_t*;
    this->dn_i = new real_t*;
}

void IonTransport::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler* 
) {
    
    real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(0, t); 
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0, t); 
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTe(0, t); 
    real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(0, t); 
    real_t dtauinvdni    = this->coefftauinv->EvaluateConfinementTime_dni(0, t);
    
    this->tauinv  = coefftauinv->EvaluateConfinementTime(0, t);
    this->dIS    = - tauinv; 
    this->dI_p    = - dtauinvdIp * n_i;
    this->dI_wall = - dtauinvdIwall * n_i;
    this->dT_cold = - dtauinvdTcold * n_i;
    this->dW_i    = - dtauinvdWi * n_i;
    this->dn_i    = - dtauinvdni * n_i;
    
    }
}

bool IonTransport::SetCSJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if (derivId != this->id_IS && derivId != this->id_Ip && derivId != this->id_Iwall && derivId != this->id_Tcold && derivId != this->id_Wi && derivId != this->id_ni)
        return false;
    
    real_t n_i = ions->GetIonDensity(0, iIon, Z0);
    
    if(derivId==id_IS){ // Detta känns fel, ska det vara id_IS för n_i? // ska index vara Z0,Z0 och Z0,0?
		jac->SetElement(rOffset+Z0, rOffset+Z0,this->dIS); 
    } else if(derivId==id_Ip){
		jac->SetElement(rOffset+Z0, 0,this->dI_p*n_i)
	} else if(derivId==id_Iwall){
		jac->SetElement(rOffset+Z0, 0,this->dI_wall*n_i)
	} else if(derivId==id_Tcold){
		jac->SetElement(rOffset+Z0, 0,this->dT_cold*n_i)
	} else if(derivId==id_Wi){
		jac->SetElement(rOffset+Z0, iz,this->dW_i*n_i)
	} else if(derivId==id_ni){
		jac->SetElement(rOffset+Z0, iz,this->dn_i*n_i)
	}
	
	return true;
}
void IonTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t* /* rhs */, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    mat[rOffset+Z0, rOffset+Z0]=-tauinv; // Är detta rätt? 
} 


void IonTransport::SetCSVectorElements(
    real_t* vec, const real_t* /* x */, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    real_t n_i     = ions->GetIonDensity(0, iIon, Z0);
    vec[rOffset+Z0]=-n_i*tauinv; // Är detta rätt?
}

