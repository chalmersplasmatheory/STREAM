#ifndef _STREAM_EQUATIONs_ION_TRANSPORT_DIFFUSION_HPP
#define _STREAM_EQUATIONs_ION_TRANSPORT_DIFFUSION_HPP

#include "DREAM/Equations/Fluid/IonEquation.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"

namespace STREAM {
	class IonTransportDiffusion : public DREAM::IonEquationTerm<FVM::DiffusionTerm> {
	private:
		DREAM::FVM::Interpolator1D *coefftauinv; // Rätt implementation av 1/tau?
		
		real_t **dDrrdni;
		real_t **dDrrdWi;
		real_t **dDrrdNi;
		//real_t **dDrrdTcold; // Behövs?
		
		len_t id_ni, id_Wi, id_Ni;
		
		void Allocate();
		void Deallocate();

	protected:
		virtual void SetCoeffs(const len_t Z0) override; // Behövs?
		//virtual void SetCoeffsAllCS(const real_t dt) override; // Behövs?
		//virtual void SetDiffCoeffsAllCS(const real_t t) override; // Behövs?
		virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/) override;
		
	public:
		IonTransportDiffusion(DREAM::FVM::Grid *g, DREAM::IonHandler *ihdl, const len_t iIon, 
			DREAM::FVM::Interpolator1D*,DREAM::FVM::MultiInterpolator1D*);
		~IonTransportDiffusion();
	}
}
#endif/*_STREAM_EQUATIONs_ION_TRANSPORT_DIFFUSION_HPP*/
