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
ConfinementTime::ConfinementTime(FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, real_t l_MK2) {
    unknowns = u;
    radials  = r;
    this->l_MK2=l_MK2;

    id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    id_Ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
}

/**
 * Evaluates the inverted confinement time
 */
real_t ConfinementTime::EvaluateConfinementTime(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return -4/(a*B*I_ref) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIwall(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    
    return 4/(a*B) *Constants::mu0*Constants::mu0*I_wall/ (M_PI*M_PI*l_MK2*l_MK2) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)/((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)*(Constants::mD)));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the electron temperature
 */
real_t ConfinementTime::EvaluateConfinementTime_dTe(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return 1/(8*a*a*B) + 2/(a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/((T_e+2/3*W_i/N_i)*(Constants::mD)));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the ion energy
 */
real_t ConfinementTime::EvaluateConfinementTime_dWi(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return 4/3*1/(a*B)*1/N_i * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/((T_e+2/3*W_i/N_i)*(Constants::mD)));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the ion density
 */
real_t ConfinementTime::EvaluateConfinementTime_dNi(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return -4/3*1/(a*B)*W_i/(N_i*N_i) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/((T_e+2/3*W_i/N_i)*(Constants::mD)));
}

/*

/**
 * Evaluates the confinement time
 * /
real_t ConfinementTime::EvaluateConfinementTime(len_t ir, real_t t){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig* /
    
    return pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD))),-1);
}

/**
 * Evaluates the derivative of the confinement time with respect to the plasma current
 * /
real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir, real_t t){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig* /
    
    return pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD))),-2)*4/(a*B*I_ref) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD));
}

/**
 * Evaluates the derivative of the confinement time with respect to the wall current
 * /
real_t ConfinementTime::EvaluateConfinementTime_dIwall(len_t ir, real_t t){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig* /
    
    
    return -pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_ii)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD))),-2)*4/(a*B) *Constants::mu0*Constants::mu0*I_wall/ (M_PI*M_PI*l_MK2*l_MK2) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)/((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)*(Constants::mD)));
}

/**
 * Evaluates the derivative of the confinement time with respect to the electron temperature
 * /
real_t ConfinementTime::EvaluateConfinementTime_dTe(len_t ir, real_t t){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig* /
    
    return -pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD))),-2)*(1/(8*a*a*B) + 2/(a*B) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/((T_e+2/3*W_i/N_i)*(Constants::mD))));
}

/**
 * Evaluates the derivative of the confinement time with respect to the ion energy
 * /
real_t ConfinementTime::EvaluateConfinementTime_dWi(len_t ir, real_t t){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig* /
    
    return -pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD))),-2)*4/3*1/(a*B)*1/N_i * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/((T_e+2/3*W_i/N_i)*(Constants::mD)));
}

/**
 * Evaluates the derivative of the confinement time with respect to the ion density
 * /
real_t ConfinementTime::EvaluateConfinementTime_dNi(len_t ir, real_t t){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t T_e    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig* /
    
    return pow((T_e/(8*a*a*B) + 4/(a*B) * exp(-I_p/I_ref) * sqrt((T_e+2/3*W_i/N_i)*(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/(Constants::mD))),-2)*(-4/3*1/(a*B)*W_i/(N_i*N_i) * exp(-I_p/I_ref) * sqrt((B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall)/((T_e+2/3*W_i/N_i)*(Constants::mD))));
}

*/
