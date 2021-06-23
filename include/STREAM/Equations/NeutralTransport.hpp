#ifndef _STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP
#define _STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP

#include "FVM/Grid/Grid.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"


namespace STREAM {
	class NeutralTransport : public DREAM::IonEquationTerm<DREAM::FVM::EquationTerm> { 
	private:
		NeutralInflux *NI;
		PlasmaVolume *PV;
		DREAM::FVM::UnknownQuantityHandler *unknowns;
		
		real_t wall_term;
		
        real_t dn_ij;
        real_t dI_p;
        real_t dI_wall;
        real_t dT_cold;
        real_t dW_i;
        real_t dN_i;
        real_t dn_cold;
        
        len_t sum_derivs=0;
        
        len_t id_Ip, id_Iwall, id_Tcold, id_Wi, id_Ni, id_ncold;
		
	public:
		NeutralTransport(DREAM::FVM::Grid *g, DREAM::IonHandler *ihdl, const len_t iIon, DREAM::FVM::UnknownQuantityHandler*, NeutralInflux*, PlasmaVolume*);
		//~NeutralTransport();
        
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
        
        virtual len_t GetNumberOfNonZerosPerRow() const override{ return 1; }
        
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return sum_derivs+6; }
	};
}
#endif/*_STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP*/
