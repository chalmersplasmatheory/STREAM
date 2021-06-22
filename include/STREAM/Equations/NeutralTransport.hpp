#ifndef _STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP
#define _STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"

namespace STREAM {
	class NeutralTransport { 
	private:
		NeutralInflux *NI;
		PlasmaVolume *PV;
		
		real_t vessel_vol;
		
	public:
		NeutralTransport(NeutralInflux*, PlasmaVolume*, real_t);
		~NeutralTransport();
        
        real_t GetNeutralTransport();
        
        /*
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
        */
	};
}
#endif/*_STREAM_EQUATIONS_NEUTRAL_TRANSPORT_HPP*/
