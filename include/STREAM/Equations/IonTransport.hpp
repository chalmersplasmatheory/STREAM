#ifndef _STREAM_EQUATIONS_ION_TRANSPORT_HPP
#define _STREAM_EQUATIONS_ION_TRANSPORT_HPP

#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"

namespace STREAM {
	class IonTransport : public DREAM::IonEquationTerm<DREAM::FVM::EquationTerm> { 
	private:
	    len_t iz;
	    ConfinementTime *coefftauinv; 
		DREAM::IonHandler *ions;
		DREAM::FVM::UnknownQuantityHandler *unknowns;
		
		real_t *dIS;//=nullptr; // Ion species
        real_t *dI_p;//=nullptr;
        real_t *dI_wall;//=nullptr;
        real_t *dT_cold;//=nullptr;
        real_t *dW_i;//=nullptr;
        real_t *dn_i;//=nullptr; // All ions
		
		
		len_t id_IS, id_Ip, id_Iwall, id_Tcold, id_Wi, id_ni;
		
		real_t tauinv;//, n_i; // Bra?
		
	public:
		IonTransport(DREAM::FVM::Grid *g, DREAM::IonHandler *ihdl,/* bool allocCoefficients,*/ const len_t iIon, len_t iz,
			ConfinementTime *tauinv, DREAM::FVM::UnknownQuantityHandler *u);
		~IonTransport();
		
        void IonTransport::Allocate();
        
        void IonTransport::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*)

		virtual bool SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        virtual void SetCSMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
	};
}
#endif/*_STREAM_EQUATIONS_ION_TRANSPORT_HPP*/
