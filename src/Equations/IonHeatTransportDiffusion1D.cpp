#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "STREAM/Equations/IonHeatTransportDiffusion3D.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
IonHeatTransportDiffusion::IonHeatTransportDiffusion(FVM::Grid *g, IonHandler *ihdl, bool allocCoefficients,
	const len_t iIon, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *u
	) : IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm>(g, ihdl, iIon, allocCoefficients), coefftauinv(tauinv), ions(ihdl) {
	
    SetName("IonHeatTransportDiffusion");
    
    this->unknowns = u;
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
    this->radials = radials;
    
    AddUnknownForJacobian(unknowns, this->id_ni);
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
IonHeatTransportDiffusion::~IonHeatTransportDiffusion(){
	Deallocate();
}

/**
 * Allocate memory for the differentiation coefficient.
 */
void IonHeatTransportDiffusion::Allocate(){
    Deallocate();
    
    len_t nzs=ions->GetNzs();

    const len_t nr = this->grid->GetNr();

    this->dn_i = new real_t*[Zion];
    this->dn_i[0]= new real_t[Zion*(nr+1)*nzs];
    for (len_t Z0=1; Z0<Zion; Z0++)
        this->dn_i[Z0]=this->dn_i[Z0-1]+(nr+1)*nzs;
    
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

/**
 * Deallocate memory for the differentiation coefficient.
 */
void IonHeatTransportDiffusion::Deallocate(){
    delete [] this->dn_i[0];
    delete [] this->dn_i;
    
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

void IonHeatTransportDiffusion::SetDiffusionTerm(const len_t Z0, real_t){ 
	if(Z0<1)
		return;
	
	real_t a = radials->GetMinorRadius();
	
	const len_t nr = this->grid->GetNr();
	
    for(len_t  ir=0; ir<nr+1; ir++){
	    real_t n_i = ions->GetIonDensity(ir, iIon, Z0);
        real_t tauinv        = this->coefftauinv->EvaluateConfinementTime(ir); 
        real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(ir); 
        real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(ir); 
        real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(ir); 
        real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(ir); 
        real_t dtauinvdNi    = this->coefftauinv->EvaluateConfinementTime_dNi(ir);
        
        this->dn_i[Z0-1][ir]    = 3/2 * Constants::ec * a * a * tauinv; 
        this->dI_p[Z0-1][ir]    = 3/2 * Constants::ec * a * a * dtauinvdIp * n_i;
        this->dI_wall[Z0-1][ir] = 3/2 * Constants::ec * a * a * dtauinvdIwall * n_i;
        this->dT_cold[Z0-1][ir] = 3/2 * Constants::ec * a * a * dtauinvdTcold * n_i;
        this->dW_i[Z0-1][ir]    = 3/2 * Constants::ec * a * a * dtauinvdWi * n_i;
        this->dN_i[Z0-1][ir]    = 3/2 * Constants::ec * a * a * dtauinvdNi * n_i;
		
		Drr(ir,0,0) += 3/2 * Constants::ec * a * a * tauinv * n_i; 
	}
}


void IonHeatTransportDiffusion::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples){
	if (derivId != this->id_ni && derivId != this->id_Ip && derivId != this->id_Iwall && derivId != this->id_Tcold && derivId != this->id_Wi && derivId != this->id_Ni)
        return;
    
    const len_t nr = this->grid->GetNr();
	
	if(derivId==id_ni){ 
		for(len_t  n=0; n<nMultiples; n++){
			for(len_t  ir=0; ir<nr+1; ir++){
				dDrr(ir,0,0,n)=dn_i[Z0ForPartials-1][ir+nr*n];
            }
        }
    }
	
	else if(derivId==id_Ip){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t  ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dI_p[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}
					
	else if(derivId==id_Iwall){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t  ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dI_wall[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}
	
	else if(derivId==id_Tcold){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t  ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dT_cold[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}
				
	else if(derivId==id_Wi){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t  ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dW_i[Z0ForPartials-1][ir];	
	            }
	        }
	    }
	}
					
	else if(derivId==id_Ni){
		for(len_t  n=0; n<nMultiples; n++){
			if(n==iIon){
				for(len_t  ir=0; ir<nr+1; ir++){
					dDrr(ir,0,0,n)=dN_i[Z0ForPartials-1][ir];
	            }
	        }
	    }
	}
}
