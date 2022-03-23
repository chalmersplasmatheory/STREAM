#include "STREAM/Equations/NeutralInflux.hpp"

using namespace STREAM;
using namespace DREAM;

/**
 * Constructor
 */
NeutralInflux::NeutralInflux(DREAM::IonHandler *ihdl, SputteredRecycledCoefficient *SRC, ConfinementTime *coefftauinv, PlasmaVolume *PV) : ions(ihdl), SRC(SRC), coefftauinv(coefftauinv), PV(PV) {
} 

/**
 * Evaluates the neutral influx
 */
real_t NeutralInflux::EvaluateNeutralInflux(const len_t iIon){
    real_t tauinv = coefftauinv->EvaluateConfinementTime(0);
    real_t V_p    = PV->GetPlasmaVolume(); 
    len_t nZ = ions->GetNZ();
    
    real_t Gamma0=0; 
    real_t n_kj = 0;
    real_t Y = 0;
    for (len_t kIon=0; kIon<nZ; kIon++) { 
        len_t Z_k   = ions->GetZ(kIon); 
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
        if (Y > 0) {
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                Gamma0 += V_p * Y * n_kj * tauinv;
            }
        }
    }
     
    return Gamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the ion species density
 */
real_t NeutralInflux::EvaluateNeutralInflux_dnkj(const len_t iIon, const len_t kIon){
    real_t tauinv = coefftauinv->EvaluateConfinementTime(0);
    real_t V_p    = PV->GetPlasmaVolume(); 
    
    real_t Y = 0;
    Y=this->SRC->GetSRCoefficient(iIon,kIon);
    
    real_t dGamma0 = V_p * Y * tauinv;
    return dGamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the plasma current
 */
real_t NeutralInflux::EvaluateNeutralInflux_dIp(const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 
    real_t dtauinvdIp = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    for (len_t kIon=0; kIon<nZ; kIon++) { 
        len_t Z_k   = ions->GetZ(kIon); 
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
        if (Y > 0) {
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                dGamma0 += V_p * Y * n_kj * dtauinvdIp;
            }
        }
    }
    
    return dGamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the wall current
 */
real_t NeutralInflux::EvaluateNeutralInflux_dIwall(const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    for (len_t kIon=0; kIon<nZ; kIon++) { 
        len_t Z_k   = ions->GetZ(kIon); 
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
        if (Y > 0) {
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                dGamma0 += V_p * Y * n_kj * dtauinvdIwall;
            }
        }
    }
    return dGamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the electron temperature
 */
real_t NeutralInflux::EvaluateNeutralInflux_dTcold(const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(0); 
    
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    for (len_t kIon=0; kIon<nZ; kIon++) { 
        len_t Z_k   = ions->GetZ(kIon); 
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
        if (Y > 0) {
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                dGamma0 += V_p * Y * n_kj * dtauinvdTcold;
            }
        }
    }
    
    return dGamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the ion energy
 */
real_t NeutralInflux::EvaluateNeutralInflux_dWi(const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdWi = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    for (len_t kIon=0; kIon<nZ; kIon++) { 
        len_t Z_k   = ions->GetZ(kIon); 
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
        if (Y > 0) {
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                dGamma0 += V_p * Y * n_kj * dtauinvdWi;
            }
        }
    }
    
    return dGamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the total ion density
 */
real_t NeutralInflux::EvaluateNeutralInflux_dNi(const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdNi = this->coefftauinv->EvaluateConfinementTime_dNi(0); 
    
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    for (len_t kIon=0; kIon<nZ; kIon++) { 
        len_t Z_k   = ions->GetZ(kIon); 
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
        if (Y > 0) {
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                dGamma0 += V_p * Y * n_kj * dtauinvdNi;
            }
        }
    }
    
    return dGamma0;
}
