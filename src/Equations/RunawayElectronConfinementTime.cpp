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
	real_t l_MK2, real_t B_v, real_t I_ref, real_t connectionLengthFactor
) {
    unknowns = u;
    radials  = r;
    this->l_MK2=l_MK2;
    this->B_v = B_v;
    this->I_ref = I_ref;
    this->connectionLengthFactor = connectionLengthFactor;
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
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);

    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t Beddy = Constants::mu0 / (M_PI*l_MK2) * I_wall;
    real_t Bz = hypot(B_v, Beddy);

    real_t Lf  = connectionLengthFactor*a/4 * B/Bz * exp(I_p / I_ref);
    real_t m   = Constants::me;
    real_t mc  = Constants::me * Constants::c;
    real_t mc2 = mc * Constants::c;
    real_t e   = Constants::ec;

    return sqrt((2*m*Lf/(e*E)) * (1 + e*E*Lf/(2*mc2)));
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
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t Beddy = Constants::mu0 / (M_PI*l_MK2) * I_wall;
    real_t Bz = hypot(B_v, Beddy);

    real_t Lf  = connectionLengthFactor*a/4 * B/Bz * exp(I_p / I_ref);
    real_t m   = Constants::me;
    real_t mc  = Constants::me * Constants::c;
    real_t mc2 = mc * Constants::c;
    real_t e   = Constants::ec;
    
    real_t tau1   = EvaluateRunawayElectronConfinementTime1(ir);
    real_t tau2   = EvaluateRunawayElectronConfinementTime2(ir);
    
    real_t dTau1_dIp = m*(1/(e*E) + Lf/mc2)/tau1 * 1/I_ref * Lf; 
    real_t dTau2_dIp = tau2 / I_p;

    return exp(-I_p/I_ref) * (1/I_ref * (1/tau2 - 1/tau1) - dTau1_dIp/(tau1*tau1) - dTau2_dIp * (exp(I_p/I_ref) - 1) / (tau2*tau2));
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t RunawayElectronConfinementTime::Evaluate_dIwall(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t Beddy = Constants::mu0 / (M_PI*l_MK2) * I_wall;
    real_t Bz = hypot(B_v, Beddy);
    
    real_t Lf  = connectionLengthFactor*a/4 * B/Bz * exp(I_p / I_ref);
    real_t m   = Constants::me;
    real_t mc  = Constants::me * Constants::c;
    real_t mc2 = mc * Constants::c;
    real_t e   = Constants::ec;
    
    real_t tau1   = EvaluateRunawayElectronConfinementTime1(ir);
    
    real_t dTau1_dIwall = m*(1/(e*E) + Lf/mc2)/tau1 * (-Constants::mu0 * Lf / (M_PI*l_MK2) * Beddy/(Bz*Bz));

    return - exp(-I_p/I_ref) / (tau1*tau1) * dTau1_dIwall;
}

real_t RunawayElectronConfinementTime::Evaluate_dE(len_t ir) {
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    real_t E      = std::abs(unknowns->GetUnknownData(id_Efield)[ir]);
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    real_t Beddy = Constants::mu0 / (M_PI*l_MK2) * I_wall;
    real_t Bz = hypot(B_v, Beddy);
    
    real_t Lf  = connectionLengthFactor*a/4 * B/Bz * exp(I_p / I_ref);
    real_t m   = Constants::me;
    real_t e   = Constants::ec;
    
    real_t tau1   = EvaluateRunawayElectronConfinementTime1(ir);
    real_t tau2   = EvaluateRunawayElectronConfinementTime2(ir);
    
    real_t dTau1_dE = - m * Lf / (e*E*E*tau1);
    real_t dTau2_dE = - tau2 / E;

    return - exp(-I_p/I_ref) / (tau1*tau1) * dTau1_dE - (1-exp(-I_p/I_ref)) / (tau2*tau2) * dTau2_dE;
}

void RunawayElectronConfinementTime::Initialize() {
    id_Ip     = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_Iwall  = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    id_Efield = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
}
