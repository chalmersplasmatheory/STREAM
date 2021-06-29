#ifndef _STREAM_EQUATIONS_CHARGE_EXCHANGE_TERM_HPP
#define _STREAM_EQUATIONS_CHARGE_EXCHANGE_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "DREAM/ADAS.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"

namespace STREAM {
    class ChargeExchangeTerm : public DREAM::FVM::DiagonalLinearTerm {
        private:
            DREAM::IonHandler *ions;
            const len_t iIon;
            DREAM::FVM::UnknownQuantityHandler *unknowns; 
            DREAM::ADAS *adas;
            PlasmaVolume *pv;
            
            real_t T_0 = 0.026;
            
            real_t id_Wi, id_Ni;
            
        protected:
            virtual len_t GetNumberOfWeightsElements() override
                {return ions->GetNzs()*grid->GetNCells();} 
            
            virtual void SetWeights() override;
        public:
            ChargeExchangeTerm(DREAM::FVM::Grid*, DREAM::IonHandler*, const len_t, DREAM::FVM::UnknownQuantityHandler*, DREAM::ADAS*, PlasmaVolume*);

            virtual len_t GetNumberOfNonZerosPerRow() const override { return this->ions->1; } 
            
            virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
            virtual void SetVectorElements(real_t*, const real_t*) override;
            
    };
}
#endif/* _STREAM_EQUATIONS_CHARGE_EXCHANGE_TERM_HPP*/
