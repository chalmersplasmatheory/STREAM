#include "STREAM/Equations/RunawayElectronConfinementTime.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
RunawayElectronConfinementTime::RunawayElectronConfinementTime(
	FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
	ConnectionLength *CL, real_t I_ref
) {
    unknowns = u;
    radials  = r;
    this->I_ref = I_ref;
    this->CL = CL;
    
    id_Ip     = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_Efield = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
}

/**
 * Evaluates the inverted confinement time.
 */
real_t RunawayElectronConfinementTime::EvaluateInverse(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t tau1 = EvaluateRunawayElectronConfinementTime1(ir);
    real_t tau2 = EvaluateRunawayElectronConfinementTime2(ir);

    return exp(-I_p/I_ref)/tau1 + (1-exp(-I_p/I_ref))/tau2;
}

/**
 * Evaluates the runaway confinement time that is dominant early
 * during plasma startup.
 */
real_t RunawayElectronConfinementTime::EvaluateRunawayElectronConfinementTime1(len_t ir) {
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);

    real_t Lfinv  = this->CL->EvaluateInverseConnectionLength(ir);
    real_t m   = Constants::me;
    real_t mc  = Constants::me * Constants::c;
    real_t mc2 = mc * Constants::c;
    real_t e   = Constants::ec;

    return sqrt((2*m/(Lfinv*e*E)) * (1 + e*E/(Lfinv*2*mc2)));
}

/**
 * Evaluates the runaway confinement time that is dominant later
 * during the plasma startup. This confinement time is estimated
 * based on the drift orbit shift of the electron.
 */
real_t RunawayElectronConfinementTime::EvaluateRunawayElectronConfinementTime2(len_t ir) {
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);
    real_t a      = radials->GetMinorRadius();
    real_t R0     = radials->GetMajorRadius();

    return R0/(10*a) * (I_p/1e6 / E);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */
real_t RunawayElectronConfinementTime::Evaluate_dIp(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);
    
    real_t Lfinv  = this->CL->EvaluateInverseConnectionLength(ir);
    real_t dLfinv_dIp  = this->CL->EvaluateInverseConnectionLength_dIp(ir);
    real_t m   = Constants::me;
    real_t mc  = Constants::me * Constants::c;
    real_t mc2 = mc * Constants::c;
    real_t e   = Constants::ec;
    
    real_t tau1   = EvaluateRunawayElectronConfinementTime1(ir);
    real_t tau2   = EvaluateRunawayElectronConfinementTime2(ir);
    
    real_t dTau1_dIp = m*(-1/(Lfinv*Lfinv*e*E) - 1/(Lfinv*Lfinv*Lfinv*mc2))/tau1 * dLfinv_dIp; 
    real_t dTau2_dIp = tau2 / I_p;

    return exp(-I_p/I_ref) * (1/I_ref * (1/tau2 - 1/tau1) - dTau1_dIp/(tau1*tau1) - dTau2_dIp * (exp(I_p/I_ref) - 1) / (tau2*tau2));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t RunawayElectronConfinementTime::Evaluate_dIwall(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);
    
    real_t Lfinv  = this->CL->EvaluateInverseConnectionLength(ir);
    real_t dLfinv_dIwall = this->CL->EvaluateInverseConnectionLength_dIwall(ir);
    real_t m   = Constants::me;
    real_t mc  = Constants::me * Constants::c;
    real_t mc2 = mc * Constants::c;
    real_t e   = Constants::ec;
    
    real_t tau1   = EvaluateRunawayElectronConfinementTime1(ir);
    
    real_t dTau1_dIwall = m*(-1/(Lfinv*Lfinv*e*E) - 1/(Lfinv*Lfinv*Lfinv*mc2))/tau1 * dLfinv_dIwall; 

    return - exp(-I_p/I_ref) / (tau1*tau1) * dTau1_dIwall;
}

real_t RunawayElectronConfinementTime::Evaluate_dE(len_t ir) {
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);

    real_t Lfinv  = this->CL->EvaluateInverseConnectionLength(ir);
    real_t m   = Constants::me;
    real_t e   = Constants::ec;
    
    real_t tau1   = EvaluateRunawayElectronConfinementTime1(ir);
    real_t tau2   = EvaluateRunawayElectronConfinementTime2(ir);
    
    real_t dTau1_dE = - m / (Lfinv*e*E*E*tau1);
    real_t dTau2_dE = - tau2 / E;

    return - exp(-I_p/I_ref) / (tau1*tau1) * dTau1_dE - (1-exp(-I_p/I_ref)) / (tau2*tau2) * dTau2_dE;
}
