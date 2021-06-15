#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonTransportDiffusion.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
IonTransportDiffusion::IonTransportDiffusion(FVM::Grid *g, IonHandler *ihdl, bool allocCoefficients,
	const len_t *iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u
	) : IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm>(g, ihdl, iIon, allocCoefficients), coefftauinv(tauinv), ions(ihdl) {
	
    SetName("IonTransportDiffusion");

    this->unknowns = unknowns;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
    this->radials = radials;
    
    AddUnknownForJacobian(unknowns, this->id_Ip);
    AddUnknownForJacobian(unknowns, this->id_Iwall);
    AddUnknownForJacobian(unknowns, this->id_Tcold);
    AddUnknownForJacobian(unknowns, this->id_Wi); 
    AddUnknownForJacobian(unknowns, this->id_ni);
	
	Allocate();
}

/**
 * Destructor.
 */
IonTransportDiffusion::~IonTransportDiffusion(){
	Deallocate();
}

/**
 * Allocate memory for the differentiation coefficient.
 */
void IonTransportDiffusion::Allocate(){
    Deallocate();
    
    len_t nzs=ions->GetNzs();

    const len_t nr = this->grid->GetNr();

    this->dI_p = new real_t[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dI_p[Z0-1] = new real_t[(nr+1)];
    
    this->dI_wall = new real_t[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dI_wall[Z0-1] = new real_t[(nr+1)];
    
    this->dT_cold = new real_t[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dT_cold[Z0-1] = new real_t[(nr+1)];
    
    this->dW_i = new real_t[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dW_i[Z0-1] = new real_t[(nr+1)];
    
    this->dn_i = new real_t[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dn_i[Z0-1] = new real_t[(nr+1)];
}

void IonTransportDiffusion::Deallocate(){
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dI_p[Z0-1];
    delete [] this->dI_p;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dI_wall;
    delete [] this->dI_wall;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dT_cold;
    delete [] this->dT_cold;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dW_i;
    delete [] this->dW_i;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dn_i;
    delete [] this->dn_i;

}

void IonTransportDiffusion::SetCoeffs(const len_t Z0){
	if(Z0<1)
		return;
	
	real_t a = radials->GetMinorRadius();
	
	const len_t nr = this->grid->GetNr();
	
    for(ir=0; ir<nr+1; ir++){
	    real_t tauinv        = this->coefftauinv->EvaluateConfinementTime(ir, t); 
        real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(ir, t); 
        real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(ir, t); 
        real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTe(ir, t); 
        real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(ir, t); 
        real_t dtauinvdni    = this->coefftauinv->EvaluateConfinementTime_dni(ir, t);
        
        // Ska det vara d...[Z0-1][ir] eller bara d...[ir]?
        this->dI_p[Z0-1][ir]    = a * a * dtauinvdIp;
        this->dI_wall[Z0-1][ir] = a * a * dtauinvdIwall;
        this->dT_cold[Z0-1][ir] = a * a * dtauinvdTcold;
        this->dW_i[Z0-1][ir]    = a * a * dtauinvdWi;
        this->dn_i[Z0-1][ir]    = a * a * dtauinvdni;
		
		Drr(ir,0,0) += a * a * tauinv * n_i; 
	}

}


void IonTransportDiffusion::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples){
	if (derivId != this->id_Ip && derivId != this->id_Iwall && derivId != this->id_Tcold && derivId != this->id_Wi && derivId != this->id_ni)
        return;
    
    // ResetDifferentiationCoefficients(); //Ska vara med?
    
    const len_t nr = this->grid->GetNr();
	
	if(derivId==id_Ip)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dI_p[Z0ForPartials-1][ir]	
					
	else if(derivId==id_Iwall)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dI_wall[Z0ForPartials-1][ir]	
	
	else if(derivId==id_Tcold)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dTcold[Z0ForPartials-1][ir]	
				
	else if(derivId==id_Wi)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dW_i[Z0ForPartials-1][ir]	
					
	else if(derivId==id_ni)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dn_i[Z0ForPartials-1][ir]		
}
