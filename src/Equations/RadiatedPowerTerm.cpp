#include "STREAM/Equations/RadiatedPowerTerm.hpp"

/**
* Implementation of class representing the radiated power.
*/

using namespace STREAM;

    RadiatedPowerTerm::RadiatedPowerTerm(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *u, DREAM::IonHandler *ionHandler, DREAM::ADAS *adas, DREAM::NIST *nist, DREAM::AMJUEL *amjuel, bool includePRB, PlasmaVolume *volumes) : DREAM::RadiatedPowerTerm(g, u, ionHandler, adas, nist, amjuel, nullptr, includePRB){ 
        //Right now we have nullptr instead of opacity, change to something like DREAM::OptionConstants::OPACITY_MODE_TRANSPARENT (but not exactly this because then it doesn't compile). We need to loop over ion type.
        this->unknowns = u;
        this->volumes = volumes;
        id_T_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD); 
        id_n_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    }
    void RadiatedPowerTerm::SetWeights(){
        const real_t V_p = this->volumes->GetPlasmaVolume();
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            const real_t V_n = this->volumes->GetNeutralVolume(iz); 
            DREAM::RadiatedPowerTerm::SetWeights(V_n/V_p);
        }    
    }
    void RadiatedPowerTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
        const real_t V_p = this->volumes->GetPlasmaVolume();
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            const real_t V_n = this->volumes->GetNeutralVolume(iz); 
            DREAM::RadiatedPowerTerm::SetDiffWeights(derivId, nMultiples, V_n/V_p);
        }    
    }
    bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x){
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            if (derivId == id_n_cold){
                real_t dV_ndN = this->volumes->GetNeutralVolume_dn(iz);
                return DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, *jac, x);
            }else if (derivId == id_T_cold){
                real_t dV_ndT = this->volumes->GetNeutralVolume_dT(iz);
                return DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, *jac, x);
            }else{
                return DREAM::RadiatedPowerTerm::SetJacobianBlock(uqtyId, derivId, *jac, x);
            }
        }     
    }
