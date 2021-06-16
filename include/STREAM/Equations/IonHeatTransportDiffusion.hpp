#ifndef _STREAM_EQUATIONS_ION_HEAT_TRANSPORT_DIFFUSION_HPP
#define _STREAM_EQUATIONS_ION_HEAT_TRANSPORT_DIFFUSION_HPP

#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"

namespace STREAM {
	class IonHeatTransportDiffusion : public DREAM::IonChargedAdvectionDiffusionTerm<DREAM::FVM::DiffusionTerm> { 
	private:
	    ConfinementTime *coefftauinv; 
		DREAM::IonHandler *ions;
		DREAM::FVM::UnknownQuantityHandler *unknowns;
        EllipticalRadialGridGenerator *radials;
		
		real_t **dIS;
        real_t **dI_p;
        real_t **dI_wall;
        real_t **dT_cold;
        real_t **dW_i;
        real_t **dn_i;
		
		len_t id_IS, id_Ip, id_Iwall, id_Tcold, id_Wi, id_ni;
		
		void Allocate();
		void Deallocate();

	protected:
		virtual void SetDiffusionTerm(const len_t Z0, real_t t);
		virtual void SetPartialDiffusionTerm(len_t /*derivId*/, len_t /*nMultiples*/) override;
		
	public:
		IonHeatTransportDiffusion(DREAM::FVM::Grid *g, DREAM::IonHandler *ihdl, bool allocCoefficients, const len_t iIon, 
			ConfinementTime *tauinv, DREAM::FVM::UnknownQuantityHandler *u);
		~IonHeatTransportDiffusion();
	};
}
#endif/*_STREAM_EQUATIONS_ION_HEAT_TRANSPORT_DIFFUSION_HPP*/
