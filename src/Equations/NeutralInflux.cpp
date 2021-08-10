#include "STREAM/Equations/NeutralInflux.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

/**
 * Constructor
 */
NeutralInflux::NeutralInflux(DREAM::IonHandler *ihdl, SputteredRecycledCoefficient *SRC, ConfinementTime *coefftauinv, PlasmaVolume *PV, real_t c1, real_t c2, real_t c3) : ions(ihdl), SRC(SRC), coefftauinv(coefftauinv), PV(PV), c1(c1), c2(c2), c3(c3) {
} 

real_t NeutralInflux::DeuteriumRecyclingCoefficient(real_t t){
    return c1-c2*(1-exp(-t/c3));
}

/**
 * Evaluates the neutral influx
 */
real_t NeutralInflux::EvaluateNeutralInflux(real_t t, const len_t iIon){
    real_t tauinv = coefftauinv->EvaluateConfinementTime(0);
    real_t V_p    = PV->GetPlasmaVolume(); 
    len_t Z_i   = ions->GetZ(iIon);  
    len_t nZ = ions->GetNZ();
    
    real_t Gamma0=0; 
    real_t n_kj = 0;
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
        n_kj = ions->GetIonDensity(0, iIon, 1);
        Gamma0 = V_p * Y * n_kj * tauinv;
    } else {
        for (len_t kIon=0; kIon<nZ; kIon++) { 
            len_t Z_k   = ions->GetZ(kIon); 
            Y=this->SRC->GetSRCoefficient(iIon,kIon);
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
real_t NeutralInflux::EvaluateNeutralInflux_dnkj(real_t t, const len_t iIon, const len_t kIon){
    real_t tauinv = coefftauinv->EvaluateConfinementTime(0);
    real_t V_p    = PV->GetPlasmaVolume(); 
    len_t Z_i   = ions->GetZ(iIon); 
    
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
    } else {
        Y=this->SRC->GetSRCoefficient(iIon,kIon);
    }
    real_t dGamma0 = V_p * Y * tauinv;
    return dGamma0;
}

/**
 * Evaluates the derivative of the neutral influx with respect to the plasma current
 */
real_t NeutralInflux::EvaluateNeutralInflux_dIp(real_t t, const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 
    real_t dtauinvdIp = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    
    len_t Z_i   = ions->GetZ(iIon); 
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
        n_kj = ions->GetIonDensity(0, iIon, 1);
        dGamma0 = V_p * Y * n_kj * dtauinvdIp;
    } else {
        for (len_t kIon=0; kIon<nZ; kIon++) { 
            len_t Z_k   = ions->GetZ(kIon); 
            Y=this->SRC->GetSRCoefficient(iIon,kIon);
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
real_t NeutralInflux::EvaluateNeutralInflux_dIwall(real_t t, const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    
    len_t Z_i   = ions->GetZ(iIon);
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
        n_kj = ions->GetIonDensity(0, iIon, 1);
        dGamma0 = V_p * Y * n_kj * dtauinvdIwall;
    } else {
        for (len_t kIon=0; kIon<nZ; kIon++) { 
            len_t Z_k   = ions->GetZ(kIon); 
            Y=this->SRC->GetSRCoefficient(iIon,kIon);
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
real_t NeutralInflux::EvaluateNeutralInflux_dTcold(real_t t, const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(0); 
    
    len_t Z_i   = ions->GetZ(iIon);
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
        n_kj = ions->GetIonDensity(0, iIon, 1);
        dGamma0 = V_p * Y * n_kj * dtauinvdTcold;
    } else {
        for (len_t kIon=0; kIon<nZ; kIon++) { 
            len_t Z_k   = ions->GetZ(kIon); 
            Y=this->SRC->GetSRCoefficient(iIon,kIon);
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
real_t NeutralInflux::EvaluateNeutralInflux_dWi(real_t t, const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdWi = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    
    len_t Z_i   = ions->GetZ(iIon);
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
        n_kj = ions->GetIonDensity(0, iIon, 1);
        dGamma0 = V_p * Y * n_kj * dtauinvdWi;
    } else {
        for (len_t kIon=0; kIon<nZ; kIon++) { 
            len_t Z_k   = ions->GetZ(kIon); 
            Y=this->SRC->GetSRCoefficient(iIon,kIon);
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
real_t NeutralInflux::EvaluateNeutralInflux_dNi(real_t t, const len_t iIon){
    real_t V_p = PV->GetPlasmaVolume(); 

    real_t dtauinvdNi = this->coefftauinv->EvaluateConfinementTime_dNi(0); 
    
    len_t Z_i   = ions->GetZ(iIon);
    len_t nZ = ions->GetNZ();
    
    real_t dGamma0=0;
    real_t n_kj = 0;
    real_t Y = 0;
    if (Z_i==1 && !ions->IsTritium(iIon)) {
        Y=DeuteriumRecyclingCoefficient(t);
        n_kj = ions->GetIonDensity(0, iIon, 1);
        dGamma0 = V_p * Y * n_kj * dtauinvdNi;
    } else {
        for (len_t kIon=0; kIon<nZ; kIon++) { 
            len_t Z_k   = ions->GetZ(kIon); 
            Y=this->SRC->GetSRCoefficient(iIon,kIon);
            for (len_t j=1; j<=Z_k; j++) {
                n_kj = ions->GetIonDensity(0, kIon, j);
                dGamma0 += V_p * Y * n_kj * dtauinvdNi;
            }
        }
    }
    return dGamma0;
}
