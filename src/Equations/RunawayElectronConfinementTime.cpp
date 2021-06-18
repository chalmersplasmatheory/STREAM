#include "STREAM/Equations/RunawayElectronConfinementTime.hpp"
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
}

/**
 * Evaluates the inverted confinement time
 */
real_t ConfinementTime::EvaluateConfinementTime(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return Constants::c/a + 4 * Constants::c/(a*B) * exp(-I_p/I_ref) * sqrt(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    return -4 * Constants::c/(a*B*I_ref) * exp(-I_p/I_ref) * sqrt(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIwall(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); /*Kanske ändrar sig*/
    
    
    return 4 * Constants::c/(a*B) *Constants::mu0*Constants::mu0*I_wall/ (M_PI*M_PI*l_MK2*l_MK2) * exp(-I_p/I_ref) * sqrt(1/(B_v*B_v+Constants::mu0*Constants::mu0/ (M_PI*M_PI*l_MK2*l_MK2)*I_wall*I_wall));
}


