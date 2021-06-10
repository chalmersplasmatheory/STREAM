#include "STREAM/Equations/ConfinementTime.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
/* Behövs alla dessa? */

using namespace STREAM
/* Rätt namespace? */

/**
 * Constructor
 */
 
ConfinementTime::ConfinementTime(FVM::UnknownQuantityHandler *u, real_t a, real_t B, real_t l_MK2) {
    unknowns = u;
    this->a = a;
    this->B = B;
    this->l_MK2=l_MK2;

    id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_Imk2    = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Ti = unknowns->GetUnknownID(OptionConstants::UQTY_T_I);
    /* Är dessa rätt? */
}

/**
 * Evaluates the inverted confinement time
 */

real_t ConfinementTime::EvaluateConfinementTime(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */

real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return -4/(a*B*I_ref) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */

real_t ConfinementTime::EvaluateConfinementTime_dIMK2(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    
    return 4/(a*B) *Constants::mu0*Constants::mu0*I_MK2/ (M_PI*M_PI*l_MK2*l_MK2) * exp(-I_p/I_ref) * sqrt((T_e+T_i)/((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)*(k_B*Constants::mD)));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the electron temperature
 */

real_t ConfinementTime::EvaluateConfinementTime_dTe(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return 1/(8*a*a*B) + 2/(a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/((T_e+T_i)*(k_B*Constants::mD)));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the ion temperature
 */

real_t ConfinementTime::EvaluateConfinementTime_dTi(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return 2/(a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/((T_e+T_i)*(k_B*Constants::mD)));
}

/*

/**
 * Evaluates the confinement time
 * /

real_t ConfinementTime::EvaluateConfinementTime(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD))),-1);
}

/**
 * Evaluates the derivative of the confinement time with respect to the plasma current
 * /

real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD))),-2)*4/(a*B*I_ref) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD));
}

/**
 * Evaluates the derivative of the confinement time with respect to the wall current
 * /

real_t ConfinementTime::EvaluateConfinementTime_dIMK2(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    
    return -pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD))),-2)*4/(a*B) *Constants::mu0*Constants::mu0*I_MK2/ (M_PI*M_PI*l_MK2*l_MK2) * exp(-I_p/I_ref) * sqrt((T_e+T_i)/((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)*(k_B*Constants::mD)));
}

/**
 * Evaluates the derivative of the confinement time with respect to the electron temperature
 * /

real_t ConfinementTime::EvaluateConfinementTime_dTe(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return -pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD))),-2)*(1/(8*a*a*B) + 2/(a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/((T_e+T_i)*(k_B*Constants::mD))));
}

/**
 * Evaluates the derivative of the confinement time with respect to the ion temperature
 * /

real_t ConfinementTime::EvaluateConfinementTime_dTi(len_t ir){
    real_t I_p = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_MK2 = unknowns->GetUnknownData(id_Imk2)[ir];
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t T_i = unknowns->GetUnknownData(id_Tcold)[ir];
    
    return -pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+T_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/(k_B*Constants::mD))),-2)*2/(a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_MK2*I_MK2)/((T_e+T_i)*(k_B*Constants::mD)));
}

*/
