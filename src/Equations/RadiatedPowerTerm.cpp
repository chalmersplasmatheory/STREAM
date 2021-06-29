#include "STREAM/Equations/RadiatedPowerTerm.hpp"

/**
* Implementation of class representing the radiated power.
*/

using namespace STREAM;

    RadiatedPowerTerm::RadiatedPowerTerm(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *u, DREAM::IonHandler *ionHandler, DREAM::ADAS *adas, DREAM::NIST *nist, DREAM::AMJUEL *amjuel, DREAM::OptionConstants::ion_opacity_mode *opacity_modes, bool includePRB, PlasmaVolume *volumes) : DREAM::RadiatedPowerTerm(g, u, ionHandler, adas, nist, amjuel, opacity_modes, includePRB), unknowns(u), ions(ionHandler), volumes(volumes){ 
        id_T_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD); 
        id_n_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    }
    void RadiatedPowerTerm::SetWeights(){
        const real_t V_p = volumes->GetPlasmaVolume();
        real_t *Vol = new real_t[ions->GetNZ()];
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            Vol[iz] = volumes->GetNeutralVolume(iz)/V_p; 
        } 
        this->DREAM::RadiatedPowerTerm::SetWeights(Vol);   
    }
    void RadiatedPowerTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
        const real_t V_p = volumes->GetPlasmaVolume();
        real_t *Vol = new real_t[ions->GetNZ()]; 
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            Vol[iz] = volumes->GetNeutralVolume(iz)/V_p; 
        } 
        this->DREAM::RadiatedPowerTerm::SetDiffWeights(derivId, nMultiples, Vol);   
    }
    bool RadiatedPowerTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t *x){
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            if (derivId == id_n_cold){
                real_t dVdn = volumes->GetNeutralVolume_dn(iz);
                jac->SetElement(iz,0,weights[iz] * dVdn); 
            }else if (derivId == id_T_cold){
                real_t dVdT = volumes->GetNeutralVolume_dT(iz);
                jac->SetElement(iz,0,weights[iz] * dVdT);
            }
        }     
        return this->DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, jac, x);
    }
