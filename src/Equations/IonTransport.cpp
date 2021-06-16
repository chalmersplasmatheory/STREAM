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
IonTransport::IonTransport(FVM::Grid *g, IonHandler *ihdl, bool allocCoefficients,
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u
	) : IonEquationTerm<DREAM::FVM::EquationTerm>(g, ihdl, iIon, allocCoefficients), coefftauinv(tauinv), ions(ihdl) {
	
    SetName("IonTransport");

    this->unknowns = u;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
    tauinv = coefftauinv->EvaluateConfinementTime(0, t); 
    n_i    = ions->GetIonDensity(ir, iIon, Z0);

    
	//Allocate();
}

/**
 * Destructor.
 */
IonTransport::~IonTransport(){}

bool IonTransport::SetCSJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t *t/* ska det vara t här?*/,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if (derivId != this->id_IS && derivId != this->id_Ip && derivId != this->id_Iwall && derivId != this->id_Tcold && derivId != this->id_Wi && derivId != this->id_ni)
        return false;
    
    real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(0, t); 
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0, t); 
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTe(0, t); 
    real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(0, t); 
    real_t dtauinvdni    = this->coefftauinv->EvaluateConfinementTime_dni(0, t);
        
    real_t dn_i    = - tauinv; 
    real_t dI_p    = - dtauinvdIp * n_i;
    real_t dI_wall = - dtauinvdIwall * n_i;
    real_t dT_cold = - dtauinvdTcold * n_i;
    real_t dW_i    = - dtauinvdWi * n_i;
    real_t dn_i    = - dtauinvdni * n_i;
    
    if(derivId==id_IS){ // Detta känns fel, ska det vara id_IS för n_i? // ska index vara Z0,Z0 och Z0,0?
		jac->SetElement(Z0,Z0,this->dn_i); 
    } else if(derivId==id_Ip){
		jac->SetElement(Z0,0,this->dI_p)
	} else if(derivId==id_Iwall){
		jac->SetElement(Z0,0,this->dI_wall)
	} else if(derivId==id_Tcold){
		jac->SetElement(Z0,0,this->dT_cold)
	} else if(derivId==id_Wi){
		jac->SetElement(Z0,Z0,this->dW_i)
	} else if(derivId==id_ni){
		jac->SetElement(Z0,Z0,this->dn_i)
	}
	
	return true;
}
void IonTransport::SetCSMatrixElements(
    FVM::Matrix *mat, real_t* /* rhs */, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    mat[Z0]=-n_i*tauinv; // Är detta rätt? Eller ska man använda rhs istället för mat? Ska de andra variablerna användas?
} 


void IonTransport::SetCSVectorElements(
    real_t* vec, const real_t* /* x */, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    vec[Z0]=-tauinv; // Är detta rätt? Eller ska man använda x istället för mat? Ska de andra variablerna användas?
}

