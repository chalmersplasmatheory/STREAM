#include "STREAM/Equations/RadiatedPowerTerm.hpp"
#include "STREAM/Settings/OptionConstants.hpp"

/**
* Implementation of class representing the radiated power. 
*/

using namespace STREAM;

    RadiatedPowerTerm::RadiatedPowerTerm(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *u, DREAM::IonHandler *ionHandler, DREAM::ADAS *adas, DREAM::NIST *nist, DREAM::AMJUEL *amjuel, DREAM::OptionConstants::ion_opacity_mode *opacity_modes, bool includePRB, PlasmaVolume *volumes) : DREAM::RadiatedPowerTerm(g, u, ionHandler, adas, nist, amjuel, opacity_modes, includePRB), unknowns(u), ions(ionHandler), volumes(volumes){ 
        id_T_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD); 
        id_n_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
        id_lambda_i = u->GetUnknownID(OptionConstants::UQTY_LAMBDA_I);

        this->Volfac = new real_t[ionHandler->GetNzs()];
        this->jacWeights = new real_t[GetNumberOfWeightsElements()];
    }
    RadiatedPowerTerm::~RadiatedPowerTerm() {
        delete [] this->jacWeights;
        delete [] this->Volfac;
    }
    void RadiatedPowerTerm::SetWeights() {
        this->SetWeights(true);
    }
    void RadiatedPowerTerm::SetWeights(bool includeIonized, real_t *w){
        const len_t NZ = ions->GetNZ();
        const real_t V_p = volumes->GetPlasmaVolume();
        for (len_t iz = 0, idx = 0; iz<NZ; iz++){
            const len_t Z = ions->GetZ(iz);
            for (len_t Z0 = 0; Z0 <= Z; Z0++, idx++) {
                if (Z0 > 0) {
                    if (includeIonized)
                        this->Volfac[idx] = 1;
                    else
                        this->Volfac[idx] = 0;
                } else
                    this->Volfac[idx] = volumes->GetNeutralVolume(iz)/V_p; 
            }
        } 
        this->DREAM::RadiatedPowerTerm::SetWeights(this->Volfac, w);
    }
    void RadiatedPowerTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
        const len_t NZ = ions->GetNZ();
        const real_t V_p = volumes->GetPlasmaVolume();
        for (len_t iz = 0, idx = 0; iz<NZ; iz++){
            const len_t Z = ions->GetZ(iz);
            for (len_t Z0 = 0; Z0 <= Z; Z0++, idx++) {
                if (Z0 > 0)
                    this->Volfac[idx] = 1;
                else
                    this->Volfac[idx] = volumes->GetNeutralVolume(iz)/V_p; 
            }
        } 
        this->DREAM::RadiatedPowerTerm::SetDiffWeights(derivId, nMultiples, this->Volfac);   
    }
    bool RadiatedPowerTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t *x){
        bool contrib = false;
        if (derivId == id_lambda_i) {
            // Rebuild weights 
            this->SetWeights(false, jacWeights);
            
            for (len_t iz = 0; iz<ions->GetNZ(); iz++){
                real_t dVdlambda = volumes->GetNeutralVolume_dLambdai(iz);

                if (jacWeights[iz] != 0)
                    jac->SetElement(iz, 0, jacWeights[iz] * dVdlambda / this->Volfac[iz]);
            }
        }
        contrib |= this->DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, jac, x);
        return contrib;
    }
