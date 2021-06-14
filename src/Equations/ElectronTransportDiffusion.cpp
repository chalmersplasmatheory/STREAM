/**
 * Implementation of a transport operator taking the form
 *
 *   
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
    EllipticalRadialGridGenerator *radials, FVM::Interpolator1D *tauinv, FVM::UnknownQuantityHandler *unknowns
) : FVM::DiffusionTerm(grid), mgtype(mgtype), coefftauinv(tauinv) { //Rätt datatyper?

    SetName("ElectronTransportDiffusion"); // Behövs?

    this->unknowns  = unknowns;
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    
    this->radials = radials;
    
    AddUnknownForJacobian(unknowns, this->id_n_cold); // Behövs?

    AllocateDiffCoeff(); // Behövs?
}

/**
 * Destructor.
 */
ElectronTransportDiffusion::~ElectronTransportDiffusion() {
    delete this->coefftauinv;
    delete [] this->dtauinv;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void ElectronTransportDiffusion::AllocateDiffCoeff() {
    const len_t nr = this->grid->GetNr();
    this->dtauinv = new real_t[nr+1];
}

/**
 * Called whenever the grid is rebuilt.
 */
bool ElectronTransportDiffusion::GridRebuilt() {
    this->FVM::DiffusionTerm::GridRebuilt();

    delete [] this->dtauinv;
    AllocateDiffCoeff();
    
    return true;
}

/**
 * Rebuild the coefficients for this equation term.
 */
  // Korrekt?
void ElectronTransportDiffusion::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns //, EllipticalRadialGridGenerator *radials Behövs denna
) {
    real_t a = radials->GetMinorRadius()->Eval(t);
    
    const real_t *tauinv = this->coefftauinv->Eval(t);
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
        this->dtauinv[ir] = 1.5 * Constants::ec * a * a * tauinv[ir];
        
        M(ir, 0, 0) += 1.5 * Constants::ec * a * a * tauinv[ir] * n; // Rätt?
    }
}

/**
 * Set jacobian of diffusion coefficients for this diffusion term.
 *
 * derivId:    ID of the quantity with respect to which the coefficient should
 *             be differentiated.
 * nMultiples: (not used).
 */
 // Korrekt?
void ElectronTransportDiffusion::SetPartialDiffusionTerm(
    len_t derivId, len_t
) {
    if (derivId != this->id_n_cold)
        return;

    ResetDifferentiationCoefficients();

    const len_t nr = this->grid->GetNr();
    for (len_t ir = 0; ir < nr+1; ir++)
        dM(ir, 0, 0) = this->dtauinv[ir];
}

