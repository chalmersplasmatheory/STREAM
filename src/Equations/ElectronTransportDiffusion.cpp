/**
 * Implementation of a heat diffusion operator taking the form
 *
 *   d/dr ( V'*D*n_cold * d/dr (T_cold) )
 *
 * This operator should be applied to 'T_cold'.
 */

#include "DREAM/Constants.hpp"
#include "STREAM/Equations/ElectronTransportDiffusion.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;
using namespace STREAM;


/**
 * Constructor.
 */

ElectronTransportDiffusion::ElectronTransportDiffusion(
    FVM::Grid *grid, enum OptionConstants::momentumgrid_type mgtype,
    EllipticalRadialGridGenerator *radials, ConfinementTime *invtau, FVM::UnknownQuantityHandler *unknowns
) : FVM::DiffusionTerm(grid), mgtype(mgtype) { //Rätt med confinementtime?

    SetName("ElectronTransportDiffusion");

    this->unknowns = unknowns;
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    
    this->radials = radials;
    this->a       = radials->GetMinorRadius()->Eval(t);
    this->coeffD  = a*a*invtau;
    
    AddUnknownForJacobian(unknowns, this->id_n_cold);
    AddUnknownForJacobian(unknowns, this->id_T_cold);

    AllocateDiffCoeff();
}

/**
 * Destructor.
 */
ElectronTransportDiffusion::~ElectronTransportDiffusion() {
    delete this->coeffD;
    delete [] this->dD;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void ElectronTransportDiffusion::AllocateDiffCoeff() {
    const len_t nr = this->grid->GetNr();
    this->dD = new real_t[nr+1];
}

/**
 * Called whenever the grid is rebuilt.
 */
bool ElectronTransportDiffusion::GridRebuilt() {
    this->FVM::DiffusionTerm::GridRebuilt();

    delete [] this->dD;
    AllocateDiffCoeff();
    
    return true;
}

/**
 * Rebuild the coefficients for this equation term.
 */
void ElectronTransportDiffusion::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const real_t *D = this->coeffD->Eval(t);
    const len_t nr = this->grid->GetNr();

    const real_t *ncold = unknowns->GetUnknownData(this->id_n_cold);

    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t n=0;
        if(ir<nr)
            n += deltaRadialFlux[ir] * ncold[ir];
        if(ir>0)
            n += (1-deltaRadialFlux[ir]) * ncold[ir-1];

        // Factor ec (=elementary charge) to convert from
        // eV to joule
        this->dD[ir] = 1.5 * Constants::ec * D[ir];
        
        Drr(ir, 0, 0) += 1.5 * Constants::ec * D[ir] * n;
    }
}

/**
 * Set jacobian of diffusion coefficients for this diffusion term.
 *
 * derivId:    ID of the quantity with respect to which the coefficient should
 *             be differentiated.
 * nMultiples: (not used).
 */
void ElectronTransportDiffusion::SetPartialDiffusionTerm(
    len_t derivId, len_t
) {
    if (derivId != this->id_n_cold)
        return;

    ResetDifferentiationCoefficients();

    const len_t nr = this->grid->GetNr();
    for (len_t ir = 0; ir < nr+1; ir++)
        dDrr(ir, 0, 0) = this->dD[ir];
}

