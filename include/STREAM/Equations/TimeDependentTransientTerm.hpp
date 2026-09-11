/**
 * Implementation of an Euler backward transient term.
 *
 * df   f_{n+1} - f_n
 * -- ~ -------------
 * dt        dt
 *
 */


#ifndef _STREAM_TIME_DEPENDENT_TRANSIENT_TERM_HPP
#define _STREAM_TIME_DEPENDENT_TRANSIENT_TERM_HPP

#include "FVM/Equation/LinearTransientTerm.hpp"

namespace STREAM {
    class TimeDependentTransientTerm : public DREAM::FVM::LinearTransientTerm {
    private:
        real_t constantScaleFactor;
		DREAM::FVM::Interpolator1D *factor;
		real_t value;

	protected:
		virtual void SetWeights() override {
			for (len_t i = 0; i < this->grid->GetNCells(); i++)
				this->weights[i] = this->constantScaleFactor * value;
		}

    public:
        TimeDependentTransientTerm(
			DREAM::FVM::Grid* g, const len_t unknownId,
			DREAM::FVM::Interpolator1D *f, real_t scaleFactor = 1.0
		) : DREAM::FVM::LinearTransientTerm(g, unknownId),
			constantScaleFactor(scaleFactor), factor(f)
		{SetName("TransientTerm");}

		virtual void Rebuild(const real_t t, const real_t dt, DREAM::FVM::UnknownQuantityHandler *uqty) override {
			this->value = factor->Eval(t)[0];
			this->LinearTransientTerm::Rebuild(t, dt, uqty);
		}
    };
}

#endif/*_STREAM_TIME_DEPENDENT_TRANSIENT_TERM_HPP*/
