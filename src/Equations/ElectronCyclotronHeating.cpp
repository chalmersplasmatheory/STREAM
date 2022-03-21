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

using namespace STREAM;
using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
ElectronCyclotronHeating::ElectronCyclotronHeating(
    DREAM::FVM::Grid *g, EllipticalRadialGridGenerator *rGrid, 
    DREAM::FVM::UnknownQuantityHandler *unknowns, OpticalThickness *OT, 
    real_t P_inj, real_t f_o, real_t f_x, real_t theta
) : EquationTerm(g), radials(rGrid), OT(OT), P_inj(P_inj), f_o(f_o), f_x(f_x), theta(theta) {

    SetName("ElectronCyclotronHeating"); 
	
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    
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
    const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*
) {
    
    real_t eta_polo = this->OT->EvaluateOpticalThickness_o(0) / cos(theta);
    real_t eta_polx = this->OT->EvaluateOpticalThickness_x(0) / cos(theta);
    this->parentheses_ECH = P_inj * (1 - f_o*exp(-eta_polo) - f_x*exp(-eta_polx));
    
    real_t dnpolo_dTe = this->OT->EvaluateOpticalThickness_o_dTe(0) / cos(theta);
    real_t dnpolx_dTe = this->OT->EvaluateOpticalThickness_x_dTe(0) / cos(theta);
    this->dpECH_dTe = P_inj * (f_o*exp(-eta_polo) * dnpolo_dTe + f_x*exp(-eta_polx) * dnpolx_dTe);
    
    real_t dnpolo_dne = this->OT->EvaluateOpticalThickness_o_dne(0) / cos(theta);
    real_t dnpolx_dne = this->OT->EvaluateOpticalThickness_x_dne(0) / cos(theta);
    this->dpECH_dne = P_inj * (f_o*exp(-eta_polo) * dnpolo_dne + f_x*exp(-eta_polx) * dnpolx_dne);
}

/**
 * Set jacobian blocks.
 */
bool ElectronCyclotronHeating::SetJacobianBlock(
    const len_t, const len_t derivId,
    DREAM::FVM::Matrix *jac, const real_t*
) {
    if (derivId == id_Tcold) {
        jac->SetElement(0, 0, this->dpECH_dTe);
    } else if (derivId == id_ncold) {
        jac->SetElement(0, 0, this->dpECH_dne);
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
    mat->SetElement(0, 0, this->parentheses_ECH);
}

/**
 * Set elements of the residual vector.
 */
void ElectronCyclotronHeating::SetVectorElements(
    real_t *vec, const real_t*
) {
    vec[0] += this->parentheses_ECH;
}

