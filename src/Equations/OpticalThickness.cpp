#include "STREAM/Equations/OpticalThickness.hpp"
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
OpticalThickness::OpticalThickness(
	FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, len_t N, real_t theta, real_t phi
) : N(N), theta(theta), phi(phi) {
    unknowns = u;
    radials  = r;
    
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
}

/**
 * Evaluates the optical thickness of the O mode 
 */
real_t OpticalThickness::EvaluateOpticalThickness_o(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t B  = radials->GetMagneticField();
    
    real_t me   = DREAM::Constants::me;
    real_t ec   = DREAM::Constants::ec;
    real_t eps0 = DREAM::Constants::eps0;
    real_t c    = DREAM::Constants::c;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = pow(M_PI,2) * R0 / (lambda * me);
    
    real_t eta_o;
    
    if (N == 1) {
	eta_o = A * T_e * n_e / n_c * sqrt(1.0 - theta) / ( 1.0 + pow(n_parallel, 2) * ( 0.5 - theta ) );
    } else {
    	real_t G = A / pow(c,2) * pow(N,3) / tgamma(N) * pow( pow(N,2) / ( 2.0 * me * pow(c,2) ), N - 2.0) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpo2 = (- E + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpo  = sqrt(n_perpo2);
    	real_t n_o2     = n_perpo2 + pow(n_parallel,2);
    	
    	real_t F_1o = pow( (S - D - n_o2) * (P - n_perpo2) , 2);
    	real_t F_2o = pow(D,2) * pow( P - n_perpo2 , 2) + pow(n_parallel,2) * P * pow( S - n_o2 , 2);
    	
    	eta_o = G * pow(T_e, N - 1) * pow(n_perpo, 2 * N - 3) * F_1o / F_2o;
    }
    
    return eta_o;
}

/**
 * Evaluates the optical thickness of the X mode 
 */
real_t OpticalThickness::EvaluateOpticalThickness_x(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t B  = radials->GetMagneticField();
    
    real_t me   = DREAM::Constants::me;
    real_t ec   = DREAM::Constants::ec;
    real_t eps0 = DREAM::Constants::eps0;
    real_t c    = DREAM::Constants::c;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = pow(M_PI,2) * R0 / (lambda * me);
    
    real_t eta_x;
    
    if (N == 1) {
	eta_x = A * T_e * pow(n_parallel,2) * pow(2.0 - theta, 1.5) * pow(1 + theta, 2) / theta;
    } else {
    	real_t G = A / pow(c,2) * pow(N,3) / tgamma(N) * pow( pow(N,2) / ( 2.0 * me * pow(c,2) ), N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpx2 = (- E - sqrt( pow(E,2) - 4.0 * S * C)) / (2 * S);
    	real_t n_perpx  = sqrt(n_perpx2);
    	real_t n_x2     = n_perpx2 + pow(n_parallel,2);
    	
    	real_t F_1x = pow( (S - D - n_x2) * (P - n_perpx2) , 2);
    	real_t F_2x = pow(D,2) * pow( P - n_perpx2 , 2) + pow(n_parallel,2) * P * pow( S - n_x2 , 2);
    	
    	eta_x = G * pow(T_e, N - 1) * pow(n_perpx, 2 * N - 3) * F_1x / F_2x;
    }
    
    return eta_x;
}

/**
 * Evaluates the derivative of the optical thickness of the O mode with regards to the electron temperature
 */
real_t OpticalThickness::EvaluateOpticalThickness_o_dTe(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t B  = radials->GetMagneticField();
    
    real_t me   = DREAM::Constants::me;
    real_t ec   = DREAM::Constants::ec;
    real_t eps0 = DREAM::Constants::eps0;
    real_t c    = DREAM::Constants::c;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = pow(M_PI,2) * R0 / (lambda * me);
    
    real_t detao_dTe;
    
    if (N == 1) {
	detao_dTe = A * n_e / n_c * sqrt(1.0 - theta) / ( 1.0 + pow(n_parallel, 2) * ( 0.5 - theta ) );
    } else {
    	real_t G = A / pow(c,2) * pow(N,3) / tgamma(N) * pow( pow(N,2) / ( 2.0 * me * pow(c,2) ), N - 2.0) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpo2 = (- E + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpo  = sqrt(n_perpo2);
    	real_t n_o2     = n_perpo2 + pow(n_parallel,2);
    	
    	real_t F_1o = pow( (S - D - n_o2) * (P - n_perpo2) , 2);
    	real_t F_2o = pow(D,2) * pow( P - n_perpo2 , 2) + pow(n_parallel,2) * P * pow( S - n_o2 , 2);
    	
    	detao_dTe = G * (N - 1.0) * pow(T_e, N - 2) * pow(n_perpo, 2 * N - 3) * F_1o / F_2o;
    }
    
    return detao_dTe;
}

/**
 * Evaluates the derivative of the optical thickness of the O mode with regards to the electron temperature
 */
 real_t OpticalThickness::EvaluateOpticalThickness_x_dTe(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t B  = radials->GetMagneticField();
    
    real_t me   = DREAM::Constants::me;
    real_t ec   = DREAM::Constants::ec;
    real_t eps0 = DREAM::Constants::eps0;
    real_t c    = DREAM::Constants::c;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = pow(M_PI,2) * R0 / (lambda * me);
    
    real_t detax_dTe;
    
    if (N == 1) {
	detax_dTe = A * pow(n_parallel,2) * pow(2.0 - theta, 1.5) * pow(1 + theta, 2) / theta;
    } else {
    	real_t G = A / pow(c,2) * pow(N,3) / tgamma(N) * pow( pow(N,2) / ( 2.0 * me * pow(c,2) ), N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpx2 = (- E - sqrt( pow(E,2) - 4.0 * S * C)) / (2 * S);
    	real_t n_perpx  = sqrt(n_perpx2);
    	real_t n_x2     = n_perpx2 + pow(n_parallel,2);
    	
    	real_t F_1x = pow( (S - D - n_x2) * (P - n_perpx2) , 2);
    	real_t F_2x = pow(D,2) * pow( P - n_perpx2 , 2) + pow(n_parallel,2) * P * pow( S - n_x2 , 2);
    	
    	detax_dTe = G * (N - 1.0) * pow(T_e, N - 2) * pow(n_perpx, 2 * N - 3) * F_1x / F_2x;
    }
    
    return detax_dTe;
}

/**
 * Evaluates the derivative of the optical thickness of the O mode with regards to the electron density
 */
 real_t OpticalThickness::EvaluateOpticalThickness_o_dne(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t B  = radials->GetMagneticField();
    
    real_t me   = DREAM::Constants::me;
    real_t ec   = DREAM::Constants::ec;
    real_t eps0 = DREAM::Constants::eps0;
    real_t c    = DREAM::Constants::c;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = pow(M_PI,2) * R0 / (lambda * me);
    
    real_t detao_dne;
    
    if (N == 1) {
	detao_dne = A * T_e * 1.0 / n_c * sqrt(1.0 - theta) / ( 1.0 + pow(n_parallel, 2) * ( 0.5 - theta ) );
    } else {
    	real_t G = A / pow(c,2) * pow(N,3) / tgamma(N) * pow( pow(N,2) / ( 2.0 * me * pow(c,2) ), N - 2.0) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t dS_dne = - 1.0 / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t dP_dne = - 1.0 / n_c;
    	real_t dD_dne = - 1.0 / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t dE_dne = (pow(n_parallel,2) - 2.0 * S - P) * dS_dne + (pow(n_parallel,2) - S) * dP_dne + 2.0 * D * dD_dne;
    	real_t dC_dne = (S - pow(n_parallel,2)) * (2.0 * P * dS_dne + (S - pow(n_parallel,2)) * dP_dne) - D * (D * dP_dne + 2 * P *dD_dne);
    	
    	real_t n_perpo2 = (- E + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpo  = sqrt(n_perpo2);
    	real_t n_o2     = n_perpo2 + pow(n_parallel,2);
    	
    	real_t dnperpo2_dne = 1.0 / sqrt(pow(E,2) - 4.0 * S * C) * ( (2.0 * S * C + E * sqrt(pow(E,2) - 4.0 * S * C) - pow(E,2)) / (2.0 * pow(S,2)) * dS_dne + (E - sqrt(pow(E,2) - 4.0 * S * C)) / (2.0 * S) * dE_dne - dC_dne );
    	real_t dnperpo_dne  = 1.0 / (2.0 * n_perpo) * dnperpo2_dne;
    	
    	real_t F_1o = pow( (S - D - n_o2) * (P - n_perpo2) , 2);
    	real_t F_2o = pow(D,2) * pow( P - n_perpo2 , 2) + pow(n_parallel,2) * P * pow( S - n_o2 , 2);
    	
    	real_t dF1o_dne = 2.0 * (S - D - n_o2) * (P - n_perpo2) * ( (P - n_perpo2) * (dS_dne - dD_dne - dnperpo2_dne) + (S - D - n_o2) * (dP_dne - dnperpo2_dne) );
    	real_t dF2o_dne = pow(n_parallel,2) * (S - n_o2) * ( 2.0 * P * (dS_dne - dnperpo2_dne) + (S - n_o2) * dP_dne ) + 2 * D * (P - n_perpo2) * ( D * (dP_dne - dnperpo2_dne) + (P - n_perpo2) * dD_dne ); 
    	
    	detao_dne = G * pow(T_e, N - 1) * pow(n_perpo, 2 * N - 4) * ( (2.0 * N - 3.0) * F_1o / F_2o * dnperpo_dne + n_perpo * (1.0 / F_2o * dF1o_dne - F_1o / pow(F_2o,2) * dF2o_dne) );
    }
    
    return detao_dne;
}

/**
 * Evaluates the derivative of the optical thickness of the X mode with regards to the electron density
 */
 real_t OpticalThickness::EvaluateOpticalThickness_x_dne(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t B  = radials->GetMagneticField();
    
    real_t me   = DREAM::Constants::me;
    real_t ec   = DREAM::Constants::ec;
    real_t eps0 = DREAM::Constants::eps0;
    real_t c    = DREAM::Constants::c;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = pow(M_PI,2) * R0 / (lambda * me);
    
    real_t detax_dne;
    
    if (N == 1) {
	detax_dne = 0.0;
    } else {
    	real_t G = A / pow(c,2) * pow(N,3) / tgamma(N) * pow( pow(N,2) / ( 2.0 * me * pow(c,2) ), N - 2.0) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t dS_dne = - 1.0 / n_c * ( pow(N,2) / ( pow(N,2) - 1.0 ) );
    	real_t dP_dne = - 1.0 / n_c;
    	real_t dD_dne = - 1.0 / n_c * ( N / ( pow(N,2) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t dE_dne = (pow(n_parallel,2) - 2.0 * S - P) * dS_dne + (pow(n_parallel,2) - S) * dP_dne + 2.0 * D * dD_dne;
    	real_t dC_dne = (S - pow(n_parallel,2)) * (2.0 * P * dS_dne + (S - pow(n_parallel,2)) * dP_dne) - D * (D * dP_dne + 2 * P *dD_dne);
    	
    	real_t n_perpx2 = (- E + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpx  = sqrt(n_perpx2);
    	real_t n_x2     = n_perpx2 + pow(n_parallel,2);
    	
    	real_t dnperpx2_dne = 1.0 / sqrt(pow(E,2) - 4.0 * S * C) * ( (-2.0 * S * C + E * sqrt(pow(E,2) - 4.0 * S * C) + pow(E,2)) / (2.0 * pow(S,2)) * dS_dne - (E + sqrt(pow(E,2) - 4.0 * S * C)) / (2.0 * S) * dE_dne + dC_dne );
    	real_t dnperpx_dne  = 1.0 / (2.0 * n_perpx) * dnperpx2_dne;
    	
    	real_t F_1x = pow( (S - D - n_x2) * (P - n_perpx2) , 2);
    	real_t F_2x = pow(D,2) * pow( P - n_perpx2 , 2) + pow(n_parallel,2) * P * pow( S - n_x2 , 2);
    	
    	real_t dF1x_dne = 2.0 * (S - D - n_x2) * (P - n_perpx2) * ( (P - n_perpx2) * (dS_dne - dD_dne - dnperpx2_dne) + (S - D - n_x2) * (dP_dne - dnperpx2_dne) );
    	real_t dF2x_dne = pow(n_parallel,2) * (S - n_x2) * ( 2.0 * P * (dS_dne - dnperpx2_dne) + (S - n_x2) * dP_dne ) + 2 * D * (P - n_perpx2) * ( D * (dP_dne - dnperpx2_dne) + (P - n_perpx2) * dD_dne ); 
    	
    	detax_dne = G * pow(T_e, N - 1) * pow(n_perpx, 2 * N - 4) * ( (2.0 * N - 3.0) * F_1x / F_2x * dnperpx_dne + n_perpx * (1.0 / F_2x * dF1x_dne - F_1x / pow(F_2x,2) * dF2x_dne) );
    }
    
    return detax_dne;
}
/*
void ConfinementTime::Initialize() {
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
}*/
