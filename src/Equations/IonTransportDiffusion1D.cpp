#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "STREAM/Equations/IonTransportDiffusion3D.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
IonTransportDiffusion::IonTransportDiffusion(FVM::Grid *g, IonHandler *ihdl, bool allocCoefficients,
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u
	) : IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm>(g, ihdl, iIon, allocCoefficients), coefftauinv(tauinv), ions(ihdl) {
	
    SetName("IonTransportDiffusion");

    this->unknowns = u;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
    this->radials = radials;
    
    AddUnknownForJacobian(unknowns, this->id_Ip);
    AddUnknownForJacobian(unknowns, this->id_Iwall);
    AddUnknownForJacobian(unknowns, this->id_Tcold);
    AddUnknownForJacobian(unknowns, this->id_Wi); 
    AddUnknownForJacobian(unknowns, this->id_Ni);
	
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

    const len_t nr = this->grid->GetNr();

    this->dI_p = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dI_p[Z0-1] = new real_t[(nr+1)];
    
    this->dI_wall = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dI_wall[Z0-1] = new real_t[(nr+1)];
    
    this->dT_cold = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dT_cold[Z0-1] = new real_t[(nr+1)];
    
    this->dW_i = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dW_i[Z0-1] = new real_t[(nr+1)];
    
    this->dN_i = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dN_i[Z0-1] = new real_t[(nr+1)];
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
        delete [] this->dN_i;
    delete [] this->dN_i;

}

void IonTransportDiffusion::SetDiffusionTerm(const len_t Z0, real_t){ 
	if(Z0<1)
		return;
	
	real_t a = radials->GetMinorRadius();
	
	const len_t nr = this->grid->GetNr();
	
    for(len_t ir=0; ir<nr+1; ir++){
	    real_t tauinv        = this->coefftauinv->EvaluateConfinementTime(ir); 
        real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(ir); 
        real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(ir); 
        real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(ir); 
        real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(ir); 
        real_t dtauinvdNi    = this->coefftauinv->EvaluateConfinementTime_dNi(ir);
        
        this->dI_p[Z0-1][ir]    = a * a * dtauinvdIp;
        this->dI_wall[Z0-1][ir] = a * a * dtauinvdIwall;
        this->dT_cold[Z0-1][ir] = a * a * dtauinvdTcold;
        this->dW_i[Z0-1][ir]    = a * a * dtauinvdWi;
        this->dN_i[Z0-1][ir]    = a * a * dtauinvdNi;
		
		Drr(ir,0,0) += a * a * tauinv; 
	}

}


void IonTransportDiffusion::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples){
	if (derivId != this->id_Ip && derivId != this->id_Iwall && derivId != this->id_Tcold && derivId != this->id_Wi && derivId != this->id_Ni)
        return;
    
    const len_t nr = this->grid->GetNr();
	
	if(derivId==id_Ip){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dI_p[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}	
					
	else if(derivId==id_Iwall){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dI_wall[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}
	
	else if(derivId==id_Tcold){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dT_cold[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}
				
	else if(derivId==id_Wi){
		for(len_t n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dW_i[Z0ForPartials-1][ir];	
	            }
	        }
	    }
	}
					
	else if(derivId==id_Ni){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dN_i[Z0ForPartials-1][ir];		
	            }
	        }
	    }
	}
}
