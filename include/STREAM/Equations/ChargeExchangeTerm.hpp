#ifndef _STREAM_EQUATIONS_CHARGE_EXCHANGE_TERM_HPP
#define _STREAM_EQUATIONS_CHARGE_EXCHANGE_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/ADAS.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"

namespace STREAM {
    class ChargeExchangeTerm : public DREAM::FVM::DiagonalComplexTerm {
        private:
            DREAM::IonHandler *ions;
            const len_t iIon;
            DREAM::ADAS *adas;
            PlasmaVolume *pv;
            EllipticalRadialGridGenerator *radials;
            
            const real_t T_0 = 0.026; // Room temperature in eV
            
            real_t id_Tcold, id_ncold, id_Wi, id_Ni, id_ni, id_lambdai;
            
            len_t D_index;
            
            virtual void AllocateDiffWeights() override;
            virtual void DeallocateDiffWeights() override;
            virtual void ResetDiffWeights() override;
            
        protected:
            virtual len_t GetNumberOfWeightsElements() override
                {return ions->GetNzs()*grid->GetNCells();} 
            
            virtual void SetWeights() override;
            virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;
        public:
            ChargeExchangeTerm(
                DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*, DREAM::IonHandler*, 
                const len_t, DREAM::ADAS*, PlasmaVolume*, EllipticalRadialGridGenerator *r, 
                DREAM::FVM::Grid *operandGrid=nullptr, len_t D_index=0
            );

            virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; } 
            
            virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
            virtual void SetVectorElements(real_t*, const real_t*) override;
            
    };
}
#endif/* _STREAM_EQUATIONS_CHARGE_EXCHANGE_TERM_HPP*/
