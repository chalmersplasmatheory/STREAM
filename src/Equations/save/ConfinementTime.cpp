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
	IonHandler *ions, ConnectionLength* CL, len_t D_index
) {
    unknowns = u;
    radials  = r;
    this->ions = ions;
    this->CL = CL;
    this->D_index = D_index;
    
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
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
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);
    
    return Lfinv * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/mi);
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
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t ec = DREAM::Constants::ec;

    real_t dLfinv_dIp = this->CL->EvaluateInverseConnectionLength_dIp(ir);
    
    return dLfinv_dIp * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIwall(len_t ir){
    len_t nr = radials->GetNr();
    
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t ec = DREAM::Constants::ec;

    real_t dLfinv_dIwall = this->CL->EvaluateInverseConnectionLength_dIwall(ir);
    
    return dLfinv_dIwall * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the electron temperature
 */
real_t ConfinementTime::EvaluateConfinementTime_dTcold(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); 
    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);

    return 1.0/(8*a*a*B) + 1.0/2.0*ec * Lfinv * 1.0/sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the ion energy
 */
real_t ConfinementTime::EvaluateConfinementTime_dWi(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);
    
    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);

    return Lfinv * 1.0/3.0*1.0/N_i * 1.0 / sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the total ion density
 */
real_t ConfinementTime::EvaluateConfinementTime_dNi(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);
    
    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);

    return -Lfinv * 1.0/3.0*W_i/(N_i*N_i) * 1.0 / sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*mi);
}


