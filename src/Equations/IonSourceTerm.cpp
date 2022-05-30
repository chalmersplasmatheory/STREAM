/**
 * Implementation of an equation term allowing for prescribed ion densities.
 */

#include "STREAM/Equations/IonSourceTerm.hpp"
#include "DREAM/MultiInterpolator1D.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace STREAM;
using namespace std;


/**
 * Constructor.
 */
IonSourceTerm::IonSourceTerm(
    DREAM::FVM::Grid *grid, DREAM::IonHandler *ihdl, const len_t nIons,
    const len_t *ionIndices, DREAM::MultiInterpolator1D *data,
	PlasmaVolume *pv, DREAM::FVM::UnknownQuantityHandler *unknowns
) : DREAM::IonSourceTerm(grid, ihdl, nIons, ionIndices, data), plasmaVolume(pv) {

    SetName("IonSourceTerm");

	this->id_lambda_i = unknowns->GetUnknownID(OptionConstants::UQTY_LAMBDA_I);
}


/**
 * Destructor.
 */
IonSourceTerm::~IonSourceTerm() {}


/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 * x:       Current value of the unknown quantity.
 *
 * (This term represents a constant, and since the derivative
 * with respect to anything of a constant is zero, we don't need
 * to do anything).
 */
bool IonSourceTerm::SetJacobianBlock(
    const len_t, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t*
) {
	// Derivative of constant = 0
	if (derivId == this->id_lambda_i) {
		for (len_t i = 0; i < nIons; i++) {
			real_t V_n_tot = this->plasmaVolume->GetTotalNeutralVolume(i);
			real_t V_n = this->plasmaVolume->GetNeutralVolume(i);
			real_t dV_n_dl = this->plasmaVolume->GetNeutralVolume_dLambdai(i);

			// Only non-trivial for Z0=0
			const len_t idx = this->ions->GetIndex(ionIndices[i], 0);
			real_t *s = currentData[i];

			jac->SetElement(idx, i, s[0]*dV_n_dl/V_n_tot * (1 - V_n/V_n_tot));
		}

		return true;
	}

	return false;
}

/**
 * Set the elements in the matrix and on the RHS corresponding
 * to this quantity.
 *
 * mat: Matrix to set elements in (1 is added to the diagonal)
 * rhs: Right-hand-side. Values will be set to the current value of
 *      this parameter.
 */
void IonSourceTerm::SetMatrixElements(DREAM::FVM::Matrix*, real_t *rhs) {
    for (len_t i = 0; i < nIons; i++) {
        for (len_t Z0 = 0; Z0 <= Z[i]; Z0++) {
			real_t Vfac;
			if (Z0 == 0)
				Vfac = this->plasmaVolume->GetNeutralVolume(i) / this->plasmaVolume->GetTotalNeutralVolume(i);
			else
				Vfac = 1;

            const len_t idx = this->ions->GetIndex(ionIndices[i], Z0);
            real_t *s = currentData[i];

			rhs[idx] += Vfac * s[Z0];
        }
    }
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * ni:  Ion densities in previous iteration.
 */
void IonSourceTerm::SetVectorElements(real_t *vec, const real_t*) {
    for (len_t i = 0; i < nIons; i++) {
		for (len_t Z0 = 0; Z0 <= Z[i]; Z0++) {
			real_t Vfac;
			if (Z0 == 0)
				Vfac = this->plasmaVolume->GetNeutralVolume(i) / this->plasmaVolume->GetTotalNeutralVolume(i);
			else
				Vfac = 1;

            const len_t idx = this->ions->GetIndex(ionIndices[i], Z0);
            real_t *s = currentData[i];

			vec[idx] += Vfac*s[Z0];
        }
    }
}

