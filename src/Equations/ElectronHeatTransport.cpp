/**
 * Electron heat transport term implemented as
 *
 *   dW/dt = 3/2 * n_e*T_e/tau
 *
 * where tau is the particle confinement time.
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "STREAM/Equations/ElectronHeatTransport.hpp"


using namespace STREAM;


/**
 * Constructor.
 */
ElectronHeatTransport::ElectronHeatTransport(
    DREAM::FVM::Grid *g, ConfinementTime *tau,
    EllipticalRadialGridGenerator *rGrid, DREAM::FVM::UnknownQuantityHandler *unknowns
) : EquationTerm(g), confinementTime(tau), radials(rGrid) {

    SetName("ElectronHeatTransportDiffusion"); 

    this->id_Ip    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_Wcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_W_COLD);
    this->id_Wi    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    
    this->radials = radials;
}

/**
 * Destructor.
 */
ElectronHeatTransport::~ElectronHeatTransport() {
}


/**
 * Rebuild coefficients for this term.
 */
void ElectronHeatTransport::Rebuild(
    const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler *unknowns
) {
    if (this->id_Iwall == 0)
        this->id_Iwall = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_WALL);

    this->invtau = this->confinementTime->EvaluateConfinementTime(0);
    this->dI_p = this->confinementTime->EvaluateConfinementTime_dIp(0);
    this->dI_wall = this->confinementTime->EvaluateConfinementTime_dIwall(0);
    this->dT_cold = this->confinementTime->EvaluateConfinementTime_dTcold(0);
    this->dW_i = this->confinementTime->EvaluateConfinementTime_dWi(0);
    this->dN_i = this->confinementTime->EvaluateConfinementTime_dNi(0);
}

/**
 * Set jacobian blocks.
 */
bool ElectronHeatTransport::SetJacobianBlock(
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
void ElectronHeatTransport::SetMatrixElements(
    DREAM::FVM::Matrix *mat, real_t*
) {
    mat->SetElement(0, 0, this->invtau);
}

/**
 * Set elements of the residual vector.
 */
void ElectronHeatTransport::SetVectorElements(
    real_t *vec, const real_t *Wcold
) {
    vec[0] += this->invtau * Wcold[0];
}

