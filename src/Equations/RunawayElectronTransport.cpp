/**
 * Implementation of the runaway electron transport term.
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "STREAM/Equations/RunawayElectronTransport.hpp"


using namespace STREAM;

/**
 * Constructor.
 */
RunawayElectronTransport::RunawayElectronTransport(
    DREAM::FVM::Grid *g, EllipticalRadialGridGenerator *r,
    RunawayElectronConfinementTime *retau, DREAM::FVM::UnknownQuantityHandler *uqh
) : EquationTerm(g), radials(r), reConfinementTime(retau) {
    
    this->id_Ip     = uqh->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    this->id_Iwall  = uqh->GetUnknownID(DREAM::OptionConstants::UQTY_I_WALL);
    this->id_Efield = uqh->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
}

/**
 * Destructor.
 */
RunawayElectronTransport::~RunawayElectronTransport() { }


/**
 * Rebuild this term.
 */
void RunawayElectronTransport::Rebuild(
    const real_t, const real_t,
    DREAM::FVM::UnknownQuantityHandler*
) {
    this->invtau = this->reConfinementTime->EvaluateInverse(0);

    // Derivatives
    this->dI_p = this->reConfinementTime->Evaluate_dIp(0);
    this->dI_wall = this->reConfinementTime->Evaluate_dIp(0);
    this->dE_field = this->reConfinementTime->Evaluate_dIp(0);
}

/**
 * Set jacobian block for this term.
 */
bool RunawayElectronTransport::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    DREAM::FVM::Matrix *jac, const real_t *nre
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(jac, nullptr);
    else if (derivId == id_Ip)
        jac->SetElement(0, 0, nre[0]*this->dI_p);
    else if (derivId == id_Iwall)
        jac->SetElement(0, 0, nre[0]*this->dI_wall);
    else if (derivId == id_Efield)
        jac->SetElement(0, 0, nre[0]*this->dE_field);
    else
        return false;

    return true;
}

/**
 * Set the linear operator matrix elements corresponding
 * to this term.
 */
void RunawayElectronTransport::SetMatrixElements(
    DREAM::FVM::Matrix *mat, real_t*
) {
    mat->SetElement(0, 0, this->invtau);
}

/**
 * Set the residual vector corresponding to this term.
 */
void RunawayElectronTransport::SetVectorElements(
    real_t *vec, const real_t *nre
) {
    vec[0] += nre[0] * this->invtau;
}

