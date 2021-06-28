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
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            const real_t V_n = volumes->GetNeutralVolume(iz); 
            DREAM::RadiatedPowerTerm::SetWeights(V_n/V_p);
        }    
    }
    void RadiatedPowerTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
        const real_t V_p = volumes->GetPlasmaVolume();
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            const real_t V_n = volumes->GetNeutralVolume(iz); 
            DREAM::RadiatedPowerTerm::SetDiffWeights(derivId, nMultiples, V_n/V_p);
        }    
    }
    bool RadiatedPowerTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t *x){
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            if (derivId == id_n_cold){
                real_t dVdn = volumes->GetNeutralVolume_dn(iz);
                jac->SetElement(iz,0,weights[iz] * dVdn); //We want to add this to the derivative of the term * V_n. Right now it looks like we overwrite it, how add instead of overwrite?
                return DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, jac, x);
            }else if (derivId == id_T_cold){
                real_t dVdT = volumes->GetNeutralVolume_dT(iz);
                jac->SetElement(iz,0,weights[iz] * dVdT);
                return DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, jac, x);
            }else{
                return DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, jac, x);
            }
        }     
    }
