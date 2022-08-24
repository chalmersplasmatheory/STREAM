#ifndef _STREAM_EQUATIONS_DISTRIBUTION_PARALLEL_TRANSPORT_HPP
#define _STREAM_EQUATIONS_DISTRIBUTION_PARALLEL_TRANSPORT_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"
#include "DREAM/Constants.hpp"
#include "STREAM/Equations/ConnectionLength.hpp"
#include "STREAM/Equations/RunawayElectronConfinementTime.hpp"

namespace STREAM {
    class DistributionParallelTransport : public DREAM::FVM::DiagonalComplexTerm {
        private:
            ConnectionLength *CL;
            RunawayElectronConfinementTime *REC;
            EllipticalRadialGridGenerator *radials;
            
            len_t id_Ip=0, id_Iwall=0, id_Efield=0; 
            
            real_t pcutoff;
                        
            virtual void AllocateDiffWeights() override;
            virtual void DeallocateDiffWeights() override;
            virtual void ResetDiffWeights() override;
            
        protected:
            virtual len_t GetNumberOfWeightsElements() override
                {return grid->GetMomentumGrid(0)->GetNp2()*grid->GetMomentumGrid(0)->GetNp1();} 
            
            virtual void SetWeights() override;
            virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;
        public:
            DistributionParallelTransport(
                DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*,
                ConnectionLength *CL, RunawayElectronConfinementTime *REC, DREAM::FVM::Grid *operandGrid, real_t
            );

            virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; } 
            
            virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
            virtual void SetVectorElements(real_t*, const real_t*) override;
            
            void Initialize();
    };
}
#endif/* _STREAM_EQUATIONS_DISTRIBUTION_PARALLEL_TRANSPORT_HPP*/
