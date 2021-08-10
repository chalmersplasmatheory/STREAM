#ifndef _STREAM_EQUATION_ION_RATE_EQUATION_HPP
#define _STREAM_EQUATION_ION_RATE_EQUATION_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"

namespace STREAM {
    class IonRateEquation : public DREAM::IonEquationTerm<DREAM::FVM::EquationTerm> {
    protected:
        enum SetMode {MATRIX, JACOBIAN};
        DREAM::ADAS *adas;
        DREAM::FVM::UnknownQuantityHandler *unknowns;
        PlasmaVolume *volumes;
        len_t id_ions, id_n_cold, id_n_tot, id_T_cold, id_lambda_i;
        bool addFluidIonization; // the full ADAS ionization rate is added in this equation term
        bool addFluidJacobian;   // only the jacobian of the ionization is set with this term
        real_t
            **Rec,         // Radiative recombination rates (nZs x nr)
            **PartialNRec, // d/dn_cold of radiative recombination rates  (nZs x nr)
            **PartialTRec, // d/dT_cold of radiative recombination rates  (nZs x nr)
            **Ion,         // Ionization rate coefficients (nZs x nr)
            **PartialNIon, // d/dn_cold of ionization rate coefficients (nZs x nr)
            **PartialTIon; // d/dT_cold of ionization rate coefficients (nZs x nr)

        // Diagnostic utilities
        real_t
            **posIonizTerm,
            **negIonizTerm,
            **posRecTerm,
            **negRecTerm,
            **posCXTerm,
            **negCXTerm;

        const bool includeChargeExchange=true;
    public:
        IonRateEquation(
            DREAM::FVM::Grid*, DREAM::IonHandler*, const len_t, DREAM::ADAS*, 
            DREAM::FVM::UnknownQuantityHandler*, PlasmaVolume*, bool,bool,bool
        );
        virtual ~IonRateEquation();

        void AllocateRateCoefficients();
        void DeallocateRateCoefficients();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 3; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            {
                len_t nnz = this->GetNumberOfNonZerosPerRow();
                nnz += 2; // 1 for ncold partial derivative and 1 for Tcold 
                return nnz; 
            }

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*) override;

        virtual bool SetCSJacobianBlock(
            const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;

        virtual void SetCSMatrixElements(
            DREAM::FVM::Matrix *mat, real_t *rhs, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override {
            SetCSMatrixElements(mat, rhs, iIon, Z0, rOffset, MATRIX);
        }

        virtual void SetCSMatrixElements(
            DREAM::FVM::Matrix*, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset, SetMode
        );

        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;

        real_t **GetPositiveIonizationTerm() { return this->posIonizTerm; }
        real_t **GetNegativeIonizationTerm() { return this->negIonizTerm; }
        real_t **GetPositiveRecombinationTerm() { return this->posRecTerm; }
        real_t **GetNegativeRecombinationTerm() { return this->negRecTerm; }
        real_t **GetPositiveChargeExchangeTerm() { return this->posCXTerm; }
        real_t **GetNegativeChargeExchangeTerm() { return this->negCXTerm; }

        len_t GetZ() { return this->Zion; }
        len_t GetIon() { return this->iIon; }
    };
}

#endif/*_STREAM_EQUATION_ION_RATE_EQUATION_HPP*/

