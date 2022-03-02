/**
 * Electron heat transport term implemented as
 *
 *   dW/dt = P_inj(1 - f_o * exp(-eta_{pol, o}) - f_x * exp(-eta_{pol, x}))
 *
 * where eta is poloidal optical thickness 
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "STREAM/Equations/ElectronCyclotronHeating.hpp"
#include <cmath>

// TODO : Add derivatives to rebuild, fix SetJacobianBlock, fix SetMatrixElements. fix SetVectorElements

using namespace STREAM;
using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
ElectronCyclotronHeating::ElectronCyclotronHeating(
    DREAM::FVM::Grid *g, EllipticalRadialGridGenerator *rGrid, 
    DREAM::FVM::UnknownQuantityHandler *unknowns
) : EquationTerm(g), radials(rGrid) {

    SetName("ElectronCyclotronHeating"); 

    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    
    this->radials = radials;
}

/**
 * Destructor.
 */
ElectronCyclotronHeating::~ElectronCyclotronHeating() {
}


/**
 * Rebuild coefficients for this term.
 */
void ElectronCyclotronHeating::Rebuild(
    const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler *unknowns
) {
    
    this->n_polo = this->confinementTime->EvaluateOpticalThickness_o(0);
    this->n_polx = this->confinementTime->EvaluateOpticalThickness_x(0);
    this->parentheses_ECH = P_inj * (1 - f_o*exp(-eta_polo) - f_x*exp(-eta_polx));
}

/**
 * Set jacobian blocks.
 */
bool ElectronCyclotronHeating::SetJacobianBlock(
    const len_t, const len_t derivId,
    DREAM::FVM::Matrix *jac, const real_t *Wcold
) {
    if (derivId == id_Ip) {
        jac->SetElement(0, 0, this->dI_p);
    } else if (derivId == id_Iwall) {
        jac->SetElement(0, 0, this->dI_wall * Wcold[0]);
    } else if (derivId == id_Ni) {
        jac->SetElement(0, 0, this->dN_i * Wcold[0]);
    } else if (derivId == id_Tcold) {
        jac->SetElement(0, 0, this->dT_cold * Wcold[0]);
    } else if (derivId == id_Wcold) {
        this->SetMatrixElements(jac, nullptr);
    } else if (derivId == id_Wi) {
        jac->SetElement(0, 0, this->dW_i * Wcold[0]);
    } else
        return false;

    return true;
}

/**
 * Set elements of the linear operator matrix.
 */
void ElectronCyclotronHeating::SetMatrixElements(
    DREAM::FVM::Matrix *mat, real_t*
) {
    mat->SetElement(0, 0, this->invtau);
}

/**
 * Set elements of the residual vector.
 */
void ElectronCyclotronHeating::SetVectorElements(
    real_t *vec, const real_t *Wcold
) {
    vec[0] += this->invtau * Wcold[0];
}

