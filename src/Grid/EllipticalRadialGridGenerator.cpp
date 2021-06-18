/**
 * Implementation of the elliptical radial grid generator.
 */

#include "DREAM/NotImplementedException.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"


using namespace STREAM;


/**
 * Constructor.
 */
EllipticalRadialGridGenerator::EllipticalRadialGridGenerator(
    DREAM::FVM::Interpolator1D *a, DREAM::FVM::Interpolator1D *B0,
    DREAM::FVM::Interpolator1D *kappa, DREAM::FVM::Interpolator1D *delta
) : RadialGridGenerator(1), a(a), B0(B0), kappa(kappa), delta(delta) {
}

/**
 * Destructor.
 */
EllipticalRadialGridGenerator::~EllipticalRadialGridGenerator() {
    delete this->delta;
    delete this->kappa;
    delete this->B0;
    delete this->a;
}


/**
 * This method indicates whether or not the grid needs to
 * be rebuilt.
 *
 * t: Time at which to check if a grid rebuild is needed.
 */
bool EllipticalRadialGridGenerator::NeedsRebuild(const real_t t) const {
    real_t
        na     = *this->a->Eval(t),
        nB0    = *this->B0->Eval(t),
        nkappa = *this->kappa->Eval(t),
        ndelta = *this->delta->Eval(t);

    // Only rebuild if values have changed...
    return (
        na != currA || nB0 != currB0 ||
        nkappa != currKappa || ndelta != currTriang
    );
}

/**
 * Rebuild this grid.
 */
bool EllipticalRadialGridGenerator::Rebuild(
    const real_t t, DREAM::FVM::RadialGrid *rg
) {
    this->currA      = *this->a->Eval(t);
    this->currB0     = *this->B0->Eval(t);
    this->currKappa  = *this->kappa->Eval(t);
    this->currTriang = *this->delta->Eval(t);

    // We always consider the point in between r=0 and r=a
    // in STREAM simulations.
    this->r = this->currA / 2;
    this->r_f[0] = 0;
    this->r_f[1] = this->currA;
    this->dr = this->currA;
    
    rg->Initialize(&this->r, this->r_f, &this->dr, &this->dr_f);

    // Reference magnetic field data
    // (TODO make new allocation here unnecessary)
    real_t R0 = std::numeric_limits<real_t>::infinity();
    real_t *BtorGOverR0 = new real_t[1];
    real_t *psiPrimeRef = new real_t[1];
    real_t *BtorGOverR0_f = new real_t[2];
    real_t *psiPrimeRef_f = new real_t[2];

    BtorGOverR0[0] = this->currB0;
    psiPrimeRef[0] = 0;     // no poloidal magnetic field

    BtorGOverR0_f[0] = BtorGOverR0_f[1] = this->currB0;
    psiPrimeRef_f[0] = psiPrimeRef_f[1] = 0;

    // Initialize radial grid with geometric data
    rg->SetReferenceMagneticFieldData(
        BtorGOverR0, BtorGOverR0_f, psiPrimeRef, psiPrimeRef_f, R0
    );

    return true;
}

/**
 * Evaluate Jacobian J/R in the given (r,\theta) point.
 */
real_t EllipticalRadialGridGenerator::JacobianAtTheta(
    const len_t, const real_t
) {
    return this->currKappa * this->r;
}
real_t EllipticalRadialGridGenerator::JacobianAtTheta_f(
    const len_t ir, const real_t
) {
    if (ir == 0) return 0;
    else return this->currKappa*this->currA;
}

/**
 * Evaluate R/R0.
 */
real_t EllipticalRadialGridGenerator::ROverR0AtTheta(
    const len_t, const real_t
) {
    return 1.0;
}
real_t EllipticalRadialGridGenerator::ROverR0AtTheta_f(
    const len_t ir, const real_t theta
) {
    return ROverR0AtTheta(ir, theta);
}

/**
 * Evaluate |grad r|^2 at the given (r,\theta).
 */
real_t EllipticalRadialGridGenerator::NablaR2AtTheta(
    const len_t, const real_t theta
) {
    real_t s=sin(theta), c=cos(theta), k=this->currKappa;
    return c*c + s*s/(k*k);
}
real_t EllipticalRadialGridGenerator::NablaR2AtTheta_f(
    const len_t ir, const real_t theta
) {
    return NablaR2AtTheta(ir, theta);
}

/**
 * Evaluate several important geometrical quantities at once.
 */
void EllipticalRadialGridGenerator::EvaluateGeometricQuantities(
    const len_t ir, const real_t theta,
    real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2
) {
    B = this->currB0;
    Jacobian = JacobianAtTheta(ir, theta);
    ROverR0 = 1.0;
    NablaR2 = NablaR2AtTheta(ir, theta);
}

/**
 * Evaluate several important geometrical quantities at once,
 * on the flux grid.
 */
void EllipticalRadialGridGenerator::EvaluateGeometricQuantities_fr(
    const len_t ir, const real_t theta,
    real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2
) {
    B = this->currB0;
    Jacobian = JacobianAtTheta_f(ir, theta);
    ROverR0 = 1.0;
    NablaR2 = NablaR2AtTheta_f(ir, theta);
}

/**
 * Calculate the radial flux label 'r' corresponding to the given
 * point in the cartesian SPI coordinate system (centred on the
 * magnetic axis).
 */
void EllipticalRadialGridGenerator::GetRThetaFromCartesian(
    real_t *r, real_t *theta, real_t x, real_t y, real_t, real_t, real_t
) {
    *r = hypot(x*this->currKappa, y) / this->currKappa;
    *theta = std::atan2(y, this->currKappa*x);
}

/**
 * Calculate the gradient of the radial flux label 'r'
 * in cartesian SPI coordinates.
 */
void EllipticalRadialGridGenerator::GetGradRCartesian(
    real_t *gradRCartesian, real_t, real_t theta
) {
    gradRCartesian[0] = cos(theta);
    gradRCartesian[1] = this->currKappa * sin(theta);
    gradRCartesian[2] = 0;
}

/**
 * Calculate the radial flux label 'r' at the point of closest
 * approach to the magnetic axis along the straight line between
 * the points given by (x1, y1, z1) and (x2, y2, z2) in the SPI
 * coordinate system.
 */
real_t EllipticalRadialGridGenerator::FindClosestApproach(
    /*real_t x1, real_t y1, real_t,
    real_t x2, real_t y2, real_t*/
    real_t, real_t, real_t, real_t, real_t, real_t
) {
    // TODO
    throw DREAM::NotImplementedException(
        "Method 'FindClosestApproach()' for SPI is not implemented in "
        "'EllipticalRadialGridGenerator' yet."
    );
}

