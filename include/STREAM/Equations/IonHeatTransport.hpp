#ifndef _STREAM_EQUATIONS_ION_HEAT_TRANSPORT_HPP
#define _STREAM_EQUATIONS_ION_HEAT_TRANSPORT_HPP

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"

namespace STREAM {
	class IonHeatTransport : public DREAM::IonEquationTerm<DREAM::FVM::EquationTerm> { 
	private:
	    ConfinementTime *coefftauinv; 
		DREAM::IonHandler *ions;
		DREAM::FVM::UnknownQuantityHandler *unknowns;
		EllipticalRadialGridGenerator *radials;
		
        real_t dI_p;
        real_t dI_wall;
        real_t dT_cold;
        real_t dW_i;
        real_t dN_i;
		
		len_t id_Ip, id_Iwall, id_Tcold, id_Wi, id_Ni;
		
		real_t tauinv, W_i;
		
		len_t D_index;
		
	public:
		IonHeatTransport(DREAM::FVM::Grid *g, DREAM::IonHandler *ihdl, const len_t iIon,
			ConfinementTime *tauinv, DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, len_t D_index=0);
		/*~IonHeatTransport();*/
        
        void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*);

		virtual bool SetCSJacobianBlock(
            const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        virtual void SetCSMatrixElements(
            DREAM::FVM::Matrix*, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; } 

        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 5; }
	};
}
#endif/*_STREAM_EQUATIONS_ION_HEAT_TRANSPORT_HPP*/
