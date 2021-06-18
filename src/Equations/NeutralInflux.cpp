#include "STREAM/Equations/NeutralInflux.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

/**
 * Constructor
 */
NeutralInflux::NeutralInflux(DREAM::IonHandler *ihdl, SputteredRecycledCoefficient *SRC, ConfinementTime *coefftauinv, /*PlasmaVolume *PV, ***Add when have PlasmaVolume class */ real_t c1, real_t c2, real_t c3) : ions(ihdl), SRC(SRC), coefftauinv(coefftauinv), /*PV(PV), ***Add when have PlasmaVolume class */ c1(c1), c2(c2), c3(c3) {
    this->tauinv = coefftauinv->EvaluateConfinementTime(0);
    /*this->V_p    = PV->GetPlasmaVolume(); ***Add when have PlasmaVolume class */
} // Korrekt sätt att hantera constructor?

real_t DeuteriumRecyclingCoefficient(real_t t){
    return c1-c2*(1-exp(-t/c3));
}

real_t EvaluateNeutralInflux(real_t t, const len_t iIon){
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * n_ij * tauinv;
        }
    }
    return Gamma0;
}

real_t EvaluateNeutralInflux_dni(real_t t, const len_t iIon){
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * tauinv;
        }
    }
    return Gamma0;
}

real_t EvaluateNeutralInflux_dIp(real_t t, const len_t iIon){
    real_t dtauinvdIp = this->coefftauinv->EvaluateConfinementTime_dIp(0); 
    
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * n_ij * dtauinvdIp;
        }
    }
    return Gamma0;
}

real_t EvaluateNeutralInflux_dIwall(real_t t, const len_t iIon){
    real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(0); 
    
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * n_ij * dtauinvdIwall;
        }
    }
    return Gamma0;
}

real_t EvaluateNeutralInflux_dTcold(real_t t, const len_t iIon){
    real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTcold(0); 
    
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * n_ij * dtauinvdTcold;
        }
    }
    return Gamma0;
}

real_t EvaluateNeutralInflux_dWi(real_t t, const len_t iIon){
    real_t dtauinvdWi = this->coefftauinv->EvaluateConfinementTime_dWi(0); 
    
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * n_ij * dtauinvdWi;
        }
    }
    return Gamma0;
}

real_t EvaluateNeutralInflux_dNi(real_t t, const len_t iIon){
    real_t dtauinvdNi = this->coefftauinv->EvaluateConfinementTime_dNi(0); 
    
    len_t Z   = ions->GetZ(const len_t iIon); // Bör denna vara en variabel i funktionen eller klassen?
    len_t *Zs = ions->GetZs(); // Är detta en array med alla aktiva ämnens Z? 
    real_t Gamma0=0; // Är detta rätt sätt att börja en summa?
    real_t n_ij = 0;
    real_t Y 0;
    for (const len_t &k : Zs) { // Är det såhär man loopar genom array?
        if (Z==1 && k==1) {
            if (ions->IsTritium(iIon)) {
                Y=0; // Ska Y=0 vid Tritium?
            } else {
                Y=DeuteriumRecyclingCoefficient(t);
            }
        } else {
            Y=this->SRC->GetSRCoefficient(Z,k);
        }
        for (Z0=1; Z0<=Z; Z0++) {
            n_ij = ions->GetIonDensity(0, iIon, Z0);
            Gamma0 += /* V_p *  ***Add when have PlasmaVolume class */Y * n_ij * dtauinvdNi;
        }
    }
    return Gamma0;
}
