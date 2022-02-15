#ifndef _STREAM_ION_SOURCE_TERM_HPP
#define _STREAM_ION_SOURCE_TERM_HPP

#include <DREAM/Equations/Fluid/IonSourceTerm.hpp>
#include "STREAM/Equations/PlasmaVolume.hpp"

namespace STREAM {
	class IonSourceTerm : public DREAM::IonSourceTerm {
	protected:
		PlasmaVolume *plasmaVolume;
		len_t id_lambda_i;

	public:
		IonSourceTerm(
			DREAM::FVM::Grid*, DREAM::IonHandler*, const len_t, const len_t*,
			DREAM::MultiInterpolator1D*, PlasmaVolume*, DREAM::FVM::UnknownQuantityHandler*
		);
		virtual ~IonSourceTerm();

		virtual bool SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
		virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
		virtual void SetVectorElements(real_t*, const real_t*) override;
	};
}

#endif/*_STREAM_ION_SOURCE_TERM_HPP*/
