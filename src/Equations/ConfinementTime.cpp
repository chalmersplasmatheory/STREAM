#include "STREAM/Equations/ConfinementTime.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

 
/**
 * Constructor
 */
ConfinementTime::ConfinementTime(
	FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
	IonHandler *ions, real_t l_MK2, real_t B_v, real_t I_ref, real_t connectionLengthFactor, len_t D_index
) {
    unknowns = u;
    radials  = r;
    this->ions = ions;
    this->D_index = D_index;
    this->l_MK2=l_MK2;
    this->B_v = B_v;
    this->I_ref = I_ref;
    this->connectionLengthFactor = connectionLengthFactor;
}

/**
 * Evaluates the inverted confinement time
 */
real_t ConfinementTime::EvaluateConfinementTime(len_t ir){ 
    return EvaluatePerpendicularConfinementTime(ir)
        + EvaluateParallelConfinementTime(ir);
}

/**
 * Evaluates the parallel confinement time.
 */
real_t ConfinementTime::EvaluateParallelConfinementTime(len_t ir) {
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    len_t nr = radials->GetNr();
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
	real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t ec = DREAM::Constants::ec;

    real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);
    
    return 4/(connectionLengthFactor*a*B) * exp(-I_p/I_ref) *
		sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*(B_v*B_v + Beddy*Beddy)/mi);
}

/**
 * Evaluates the perpendicular confinement time.
 */
real_t ConfinementTime::EvaluatePerpendicularConfinementTime(len_t ir) {
    real_t a      = radials->GetMinorRadius();
    real_t B      = radials->GetMagneticField();
    real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
    return T_cold/(8*a*a*B);
}

/**
 * Evaluates the effective connection length.
 */
real_t ConfinementTime::EvaluateConnectionLength(len_t ir) {
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();

    real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return 0.25 * connectionLengthFactor * a*B * exp(I_p/I_ref) / (sqrt(B_v*B_v + Beddy*Beddy));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    len_t nr = radials->GetNr();
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
	real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t ec = DREAM::Constants::ec;
    
	real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return -4/(connectionLengthFactor*a*B*I_ref) * exp(-I_p/I_ref) * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*(B_v*B_v+Beddy*Beddy)/mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIwall(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    len_t nr = radials->GetNr();
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
	real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t ec = DREAM::Constants::ec;
    
	real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return 4/(connectionLengthFactor*a*B) *Constants::mu0*Constants::mu0*I_wall/ (2*2*M_PI*M_PI*l_MK2*l_MK2) * exp(-I_p/I_ref) * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/((B_v*B_v+Beddy*Beddy)*mi));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the electron temperature
 */
real_t ConfinementTime::EvaluateConfinementTime_dTcold(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    len_t nr = radials->GetNr();
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
	real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); 
    real_t ec = DREAM::Constants::ec;
    
	real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return 1.0/(8*a*a*connectionLengthFactor*B) + 2*ec/(connectionLengthFactor*a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Beddy*Beddy)/((ec*T_cold+2.0/3.0*W_i/N_i)*mi));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the ion energy
 */
real_t ConfinementTime::EvaluateConfinementTime_dWi(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    len_t nr = radials->GetNr();
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
	real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); 
    real_t ec = DREAM::Constants::ec;
    
	real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return 4/3.0*1/(connectionLengthFactor*a*B)*1/N_i * exp(-I_p/I_ref) * sqrt((B_v*B_v+Beddy*Beddy)/((ec*T_cold+2.0/3.0*W_i/N_i)*mi));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the total ion density
 */
real_t ConfinementTime::EvaluateConfinementTime_dNi(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    len_t nr = radials->GetNr();
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
	real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); 
    real_t ec = DREAM::Constants::ec;
    
	real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return -4/3.0*1/(connectionLengthFactor*a*B)*W_i/(N_i*N_i) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Beddy*Beddy)/((ec*T_cold+2.0/3.0*W_i/N_i)*mi));
}

/**
 * Get IDs for unknowns (since I_wall is defined late during the construction of
 * the equation system this must be done seperately from other initialization)
 */
void ConfinementTime::Initialize() {
    id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
}

