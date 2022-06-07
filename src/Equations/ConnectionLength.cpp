#include "STREAM/Equations/ConnectionLength.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

 
/**
 * Constructor
 */
ConnectionLength::ConnectionLength(
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
 * Evaluates the inverted connection length
 */
real_t ConnectionLength::EvaluateInverseConnectionLength(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    
    real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return 4.0 * sqrt(B_v*B_v + Beddy*Beddy) / (connectionLengthFactor*a*B) * exp(-I_p/I_ref);
}

/**
 * Evaluates the effective connection length.
 */
real_t ConnectionLength::EvaluateConnectionLength(len_t ir) {
    return 1 / EvaluateInverseConnectionLength(ir);
}

/**
 * Evaluates the derivative of the inverted connection length with respect to the plasma current
 */
real_t ConnectionLength::EvaluateInverseConnectionLength_dIp(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    
    real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return -1.0/I_ref * 4.0 * sqrt(B_v*B_v + Beddy*Beddy) / (connectionLengthFactor*a*B) * exp(-I_p/I_ref);
}

/**
 * Evaluates the derivative of the inverted connection length with respect to the wall current
 */
real_t ConnectionLength::EvaluateInverseConnectionLength_dIwall(len_t ir){
    real_t I_p    = unknowns->GetUnknownData(id_Ip)[ir];
    real_t I_wall = unknowns->GetUnknownData(id_Iwall)[ir];
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField();
    
    real_t Beddy = Constants::mu0*I_wall / (2*M_PI*l_MK2);

    return Constants::mu0*Constants::mu0*I_wall / (M_PI*M_PI*l_MK2*l_MK2) / (connectionLengthFactor*a*B*sqrt(B_v*B_v + Beddy*Beddy)) * exp(-I_p/I_ref);
}

/**
 * Get IDs for unknowns (since I_wall is defined late during the construction of
 * the equation system this must be done seperately from other initialization)
 */
void ConnectionLength::Initialize() {
    id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
}

