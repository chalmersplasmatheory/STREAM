#ifndef _STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP
#define _STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"

//TODO Använda Esmées redan existerande kod för gamma_{n,i}

namespace STREAM {
	class NeutralTransport : public DREAM::IonEquationTerm<DREAM::FVM::EquationTerm> { 
	private:
		NeutralInflux *NI;
		PlasmaVolume *PV;
	    ConfinementTime *CT; 
		
		real_t vessel_vol;
		
		real_t wall_term;
		real_t tauinv;
		
        real_t dI_p;
        real_t dI_wall;
        real_t dT_cold;
        real_t dW_i;
        real_t dN_i;
        
		
	public:
		NeutralTransport(DREAM::FVM::Grid *g, DREAM::IonHandler *ihdl, const len_t iIon, NeutralInflux*, PlasmaVolume*, ConfinementTime*, real_t);
		~NeutralTransport();
        
        void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*);
        
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
        
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 6; }
	};
}
#endif/*_STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP*/
