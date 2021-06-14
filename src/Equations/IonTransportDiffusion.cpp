#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;
using namespace STREAM;

IonTransportDiffusion::IonTransportDiffusion(FVM::Grid *g, IonHandler *ihdl,
	const len_t iIon, FVM::Interpolator1D* coefftauinv, FVM::MultiInterpolator1D* DrrHat, FVM::UnknownQuantityHandler *u
	) : IonEquationTerm<FVM::DiffusionTerm>(FVM::Grid *g, 
	IonHandler *ihdl, const len_t iIon), coefftauinv(coefftauinv), DrrHat(DrrHat) {
	
    SetName("IonTransportsDiffusion");

    this->id_ni = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    AddUnknownForJacobian(u, id_ni);
    
    this->id_Wi = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    AddUnknownForJacobian(u, id_Wi);
    
    this->id_Ni = u->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    AddUnknownForJacobian(u, id_Ni);
	
	Allocate();
}

IonTransportDiffusion::~IonTransportDiffusion(){
	Deallocate();
}

void IonTransportDiffusion::Allocate(){
    Deallocate();
    
    len_t nzs=ihdl->GetNzs(); // Vad?

    const len_t nr = this->g->GetNr();

	// Derivatives wrt ion (charge state) densities
    this->dDrrdni = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dDrrdni[Z0-1] = new real_t[(nr+1)*nzs];
        
    // Derivatives wrt ion thermal energy (of this species only)
    this->dDrrdWi = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dDrrdWi[Z0-1] = new real_t[(nr+1)];
        
    // Derivatives wrt the total ion density (of this species only)
    this->dDrrdNi = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dDrrdNi[Z0-1] = new real_t[(nr+1)];
}

void IonTransportDiffusion::Deallocate(){
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dDrrdni[Z0-1];
    delete [] this->dDrrdni;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dDrrdWi[Z0-1];
    delete [] this->dDrrdWi;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dDrrdNi[Z0-1];
    delete [] this->dDrrdNi;
}

/*
void IonTransportDiffusion::SetCoeffsAllCS(const real_t t){
}

void IonTransportDiffusion::SetDiffCoeffsAllCS(const real_t t){
}
*/

void IonTransportDiffusion::SetCoeffs(const len_t Z0){
	if(Z0<1)
		return;
	
	real_t a = radials->GetMinorRadius()->Eval(t);
    const real_t *tauinv = this->coefftauinv->Eval(t);
	const len_t nr = this->g->GetNr();
	const real_t *W_i = unknowns->GetUnknownData(this->id_Wi); // Hitta W_i, eller T_i och n_i för givet Z0??
    const real_t *N
    _i = unknowns->GetUnknownData(this->id_Ni);
	const real_t *T_i = 2/3*W_i/N_i;
	for(ir=0; ir<nr+1; ir++)
        // Fel interpolation va?
    	real_t T=0; 
    	real_t n=0; 
        if(ir<nr)
            T += deltaRadialFlux[ir] * T_i[ir];
            n += deltaRadialFlux[ir] * n_i[ir];
            // W += deltaRadialFlux[ir] * W_i[ir];
        if(ir>0)
            T += (1-deltaRadialFlux[ir]) * T_i[ir-1];
            n += (1-deltaRadialFlux[ir]) * n_i[ir-1];
            W += (1-deltaRadialFlux[ir]) * W_i[ir-1];
		Drr(ir,0,0)=3/2 * n * T /* * W */ * a * a * tauinv;
}


void IonTransportDiffusion::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples){
	if(derivId==id_ni)
		for(n=0; n<nMultiples; n++)
			for(ir=0; ir<nr+1; ir++)
				dDrr(ir,0,0,n)=dDrrdni[Z0ForPartials-1][ir+nr*n]
				
	if(derivId==id_Wi)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dDrrdWi[Z0ForPartials-1][ir]	
					
	if(derivId==id_Ni)
		for(n=0; n<nMultiples; n++)
			if(n==iIon)
				for(ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dDrrdNi[Z0ForPartials-1][ir]			
}
