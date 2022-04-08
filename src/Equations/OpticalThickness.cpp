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

real_t OpticalThickness::EvaluateTau(len_t ir) {
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R0 = radials->GetMajorRadius();
    real_t a  = radials->GetMinorRadius();
    
    real_t R_a = R0 + a;
    real_t N_par2 = sin(phi) * R_a / 5.4 * sin(phi) * R_a / 5.4;
    real_t tau = 6.6e-22 * n_e * T_e * pow(1 - N_par2, 1.5) / cos(theta);
    
    //printf("\nSCENPLINT:   tau=%.4e", tau*cos(theta));
    return tau;
}

real_t OpticalThickness::EvaluateTau_dTe(len_t ir) {
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t R_a = 1.8;
    real_t N_par2 = sin(phi) * R_a / 5.4 * sin(phi) * R_a / 5.4;
    real_t tau = 6.6e-22 * n_e * pow(1 - N_par2, 1.5) / cos(theta);
    
    return tau;
}

real_t OpticalThickness::EvaluateTau_dne(len_t ir) {
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    
    real_t R_a = 1.8;
    real_t N_par2 = sin(phi) * R_a / 5.4 * sin(phi) * R_a / 5.4;
    real_t tau = 6.6e-22 * T_e * pow(1 - N_par2, 1.5) / cos(theta);
    
    return tau;
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
    real_t mc2  = DREAM::Constants::mc2inEV;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * (omega*omega) / (ec*ec);
    real_t n_par = cos(phi);
    real_t n_par2 = n_par*n_par;
    real_t n_par4 = n_par2*n_par2;
    len_t N2 = N*N;
    len_t N3 = N2*N;
    real_t alpha = n_e / n_c;
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t eta_o;
    
    if (N == 1) {
	eta_o = A * (T_e / mc2) * alpha * sqrt(1.0 - alpha) / ( 1.0 + n_par2 * ( 0.5 - alpha ) );
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2);
    	
    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	real_t P = 1.0 - alpha;
    	real_t D = - alpha * (N / (N2 - 1.0 ));
    	
    	real_t denom = (N2 - 1.0 - alpha * N2);
    	real_t Q_sqrt = sqrt(alpha2 * (1.0 + n_par4 + (4.0 * N2 - 2.0) * n_par2) - 4.0 * alpha3 * N2 * n_par2);
    	real_t Q = - Q_sqrt / (2.0 * denom);
    	
    	real_t R_1 = alpha * (N2 - 1.0) * (2.0 * N - 1.0 + n_par2);
    	real_t R_2 = - 2.0 * alpha2 * N2 * (N - 1.0);
    	real_t R = (R_1 + R_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t T_0 = 2.0 * (N2 - 1.0) * n_par2;
    	real_t T_1 = alpha * (1.0 - (2.0 * N2 - 1.0) * n_par2);
    	real_t T = (T_0 + T_1) / (2.0 * denom);
    	
    	real_t V_1 = - alpha * (N2 - 1.0) * (1.0 - n_par2);
    	real_t V_2 = 2.0 * alpha2 * N2;
    	real_t V = (V_1 + V_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t n_perp_0 = 2.0 * (N2 - 1) * (1.0 - n_par2);
    	real_t n_perp_1 = - alpha * (4.0 * N2 - 1 - (2.0 * N2 - 1.0) * n_par2);
    	real_t n_perp_2 = 2 * alpha2 * N2;
    	real_t n_perp2 = (n_perp_0 + n_perp_1 + n_perp_2) / (2.0 * denom) - Q;
    	real_t n_perp = sqrt(n_perp2);
    	
    	real_t F_1 = (R + Q) * (R + Q) * (T + Q) * (T + Q);
    	real_t F_2 = (D*D) * ((T + Q)*(T + Q)) + n_par2 * P *((V + Q)*(V + Q)); 
    	
    	eta_o = G * pow((T_e / mc2), N - 1.0) * pow(n_perp, 2 * N - 3) * alpha * F_1 / F_2;
    	//printf("\nDYON:      eta_o=%.4e", eta_o);
    	//printf("\n For o:     eta_o=%.3e, n_perp=%.3e, F_1=%.3e, F_2=%.3e, R=%.3e, Q=%.3e", eta_o, n_perp, F_1, F_2, R, Q);
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
    real_t mc2  = DREAM::Constants::mc2inEV;
    
    real_t omega_ce = ec * B / me; 
    real_t omega = N * omega_ce; 
    real_t lambda = 2.0 * M_PI * c / omega; 
    real_t n_c = me * eps0 * (omega*omega) / (ec*ec); 
    real_t n_par = cos(phi);
    real_t n_par2 = n_par*n_par;
    real_t n_par4 = n_par2*n_par2;
    len_t N2 = N*N;
    len_t N3 = N2*N;
    real_t alpha = n_e / n_c;
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t eta_x;
    
    if (N == 1) {
	eta_x = A * (T_e / mc2) * n_par2 * pow(2.0 - alpha, 1.5) * ((1.0 + alpha)*(1.0 + alpha)) / alpha;
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2); 

    	real_t alpha = n_e / n_c;
    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	real_t P = 1.0 - alpha;
    	real_t D = - alpha * (N / (N2 - 1.0 ));
    	
    	real_t denom = (N2 - 1.0 - alpha * N2);
    	real_t Q_sqrt = sqrt(alpha2 * (1.0 + n_par4 + (4.0 * N2 - 2.0) * n_par2) - 4.0 * alpha3 * N2 * n_par2);
    	real_t Q = Q_sqrt / (2.0 * denom);
    	
    	real_t R_1 = alpha * (N2 - 1.0) * (2.0 * N - 1.0 + n_par2);
    	real_t R_2 = - 2.0 * alpha2 * N2 * (N - 1.0);
    	real_t R = (R_1 + R_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t T_0 = 2.0 * (N2 - 1.0) * n_par2;
    	real_t T_1 = alpha * (1.0 - (2.0 * N2 - 1.0) * n_par2);
    	real_t T = (T_0 + T_1) / (2.0 * denom);
    	
    	real_t V_1 = - alpha * (N2 - 1.0) * (1.0 - n_par2);
    	real_t V_2 = 2.0 * alpha2 * N2;
    	real_t V = (V_1 + V_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t n_perp_0 = 2.0 * (N2 - 1) * (1.0 - n_par2);
    	real_t n_perp_1 = - alpha * (4.0 * N2 - 1 - (2.0 * N2 - 1.0) * n_par2);
    	real_t n_perp_2 = 2 * alpha2 * N2;
    	real_t n_perp2 = (n_perp_0 + n_perp_1 + n_perp_2) / (2.0 * denom) - Q;
    	real_t n_perp = sqrt(n_perp2);
    	
    	real_t F_1 = (R + Q) * (R + Q) * (T + Q) * (T + Q);
    	real_t F_2 = (D*D) * ((T + Q)*(T + Q)) + n_par2 * P *((V + Q)*(V + Q)); 
    	
    	eta_x = G * pow((T_e / mc2), N - 1.0) * pow(n_perp, 2 * N - 3) * alpha * F_1 / F_2; 
    	//printf("\nDYON:      eta_x=%.4e", eta_x);
    	//printf("\n For x:     eta_x=%.3e, G=%.3e, n_perp=%.3e, F_1=%.3e, F_2=%.3e, R=%.3e, Q=+%.3e", eta_x,, n_perp, F_1, F_2, R, Q);
    	//printf("\n n_e/n_c=%.4e", n_e/n_c);
	
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
    real_t mc2  = DREAM::Constants::mc2inEV;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * (omega*omega) / (ec*ec);
    real_t n_par = cos(phi);
    real_t n_par2 = n_par*n_par;
    real_t n_par4 = n_par2*n_par2;
    len_t N2 = N*N;
    len_t N3 = N2*N;
    real_t alpha = n_e / n_c;
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detao_dTe;
    
    if (N == 1) {
	detao_dTe = A * (1.0 / mc2) * alpha * sqrt(1.0 - alpha) / ( 1.0 + n_par2 * ( 0.5 - alpha ) );
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2);

    	real_t alpha = n_e / n_c;
    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	real_t P = 1.0 - alpha;
    	real_t D = - alpha * (N / (N2 - 1.0 ));
    	
    	real_t denom = (N2 - 1.0 - alpha * N2);
    	real_t Q_sqrt = sqrt(alpha2 * (1.0 + n_par4 + (4.0 * N2 - 2.0) * n_par2) - 4.0 * alpha3 * N2 * n_par2);
    	real_t Q = - Q_sqrt / (2.0 * denom);
    	
    	real_t R_1 = alpha * (N2 - 1.0) * (2.0 * N - 1.0 + n_par2);
    	real_t R_2 = - 2.0 * alpha2 * N2 * (N - 1.0);
    	real_t R = (R_1 + R_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t T_0 = 2.0 * (N2 - 1.0) * n_par2;
    	real_t T_1 = alpha * (1.0 - (2.0 * N2 - 1.0) * n_par2);
    	real_t T = (T_0 + T_1) / (2.0 * denom);
    	
    	real_t V_1 = - alpha * (N2 - 1.0) * (1.0 - n_par2);
    	real_t V_2 = 2.0 * alpha2 * N2;
    	real_t V = (V_1 + V_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t n_perp_0 = 2.0 * (N2 - 1) * (1.0 - n_par2);
    	real_t n_perp_1 = - alpha * (4.0 * N2 - 1 - (2.0 * N2 - 1.0) * n_par2);
    	real_t n_perp_2 = 2 * alpha2 * N2;
    	real_t n_perp2 = (n_perp_0 + n_perp_1 + n_perp_2) / (2.0 * denom) - Q;
    	real_t n_perp = sqrt(n_perp2);
    	
    	real_t F_1 = (R + Q) * (R + Q) * (T + Q) * (T + Q);
    	real_t F_2 = (D*D) * ((T + Q)*(T + Q)) + n_par2 * P *((V + Q)*(V + Q)); 
    	
    	detao_dTe = G * (N - 1.0) / mc2 * pow((T_e / mc2), N - 2) * pow(n_perp, 2 * N - 3) * alpha * F_1 / F_2;
    }
    
    return detao_dTe;
}

/**
 * Evaluates the derivative of the optical thickness of the X mode with regards to the electron temperature
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
    real_t mc2  = DREAM::Constants::mc2inEV;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * (omega*omega) / (ec*ec);
    real_t n_par = cos(phi);
    real_t n_par2 = n_par*n_par;
    real_t n_par4 = n_par2*n_par2;
    len_t N2 = N*N;
    len_t N3 = N2*N;
    real_t alpha = n_e / n_c;
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detax_dTe;
    
    if (N == 1) {
	detax_dTe = A * 1.0 / mc2 * n_par2 * pow(2.0 - alpha, 1.5) * ((1.0 + alpha)*(1.0 + alpha)) / alpha;
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2);

    	real_t alpha = n_e / n_c;
    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	real_t P = 1.0 - alpha;
    	real_t D = - alpha * (N / (N2 - 1.0 ));
    	
    	real_t denom = (N2 - 1.0 - alpha * N2);
    	real_t Q_sqrt = sqrt(alpha2 * (1.0 + n_par4 + (4.0 * N2 - 2.0) * n_par2) - 4.0 * alpha3 * N2 * n_par2);
    	real_t Q = Q_sqrt / (2.0 * denom);
    	
    	real_t R_1 = alpha * (N2 - 1.0) * (2.0 * N - 1.0 + n_par2);
    	real_t R_2 = - 2.0 * alpha2 * N2 * (N - 1.0);
    	real_t R = (R_1 + R_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t T_0 = 2.0 * (N2 - 1.0) * n_par2;
    	real_t T_1 = alpha * (1.0 - (2.0 * N2 - 1.0) * n_par2);
    	real_t T = (T_0 + T_1) / (2.0 * denom);
    	
    	real_t V_1 = - alpha * (N2 - 1.0) * (1.0 - n_par2);
    	real_t V_2 = 2.0 * alpha2 * N2;
    	real_t V = (V_1 + V_2) / (2.0 * (N2 - 1.0) * denom);
    	
    	real_t n_perp_0 = 2.0 * (N2 - 1) * (1.0 - n_par2);
    	real_t n_perp_1 = - alpha * (4.0 * N2 - 1 - (2.0 * N2 - 1.0) * n_par2);
    	real_t n_perp_2 = 2 * alpha2 * N2;
    	real_t n_perp2 = (n_perp_0 + n_perp_1 + n_perp_2) / (2.0 * denom) - Q;
    	real_t n_perp = sqrt(n_perp2);
    	
    	real_t F_1 = (R + Q) * (R + Q) * (T + Q) * (T + Q);
    	real_t F_2 = (D*D) * ((T + Q)*(T + Q)) + n_par2 * P *((V + Q)*(V + Q)); 
    	
    	detax_dTe = G * (N - 1.0) / mc2 * pow((T_e / mc2), N - 2) * pow(n_perp, 2 * N - 3) * alpha * F_1 / F_2;
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
    real_t mc2  = DREAM::Constants::mc2inEV;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * (omega*omega) / (ec*ec);
    real_t n_par = cos(phi);
    real_t n_par2 = n_par*n_par;
    real_t n_par4 = n_par2*n_par2;
    len_t N2 = N*N;
    len_t N3 = N2*N;
    len_t N4 = N2*N2;
    real_t alpha = n_e / n_c;
    real_t dalpha_dne = 1 / n_c;
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detao_dne;
    
    if (N == 1) {
    	real_t dterm1 = dalpha_dne * sqrt(1.0 - alpha) / (1.0 + n_par2 * ( 0.5 - alpha ));
    	real_t dterm2 = - alpha * 1.0 / (2.0 * sqrt(1.0 - alpha)) / (1.0 + n_par2 * ( 0.5 - alpha ));
    	real_t dterm3 = alpha * sqrt(1.0 - alpha) / ((1.0 + n_par2 * ( 0.5 - alpha ))*(1.0 + n_par2 * ( 0.5 - alpha ))) * n_par2;
	detao_dne = A * (T_e / mc2) * (dterm1 + dterm2 + dterm3);
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2);

    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	
    	real_t P = 1.0 - alpha;
    	real_t dP_dne = - dalpha_dne;
    	
    	real_t D = - alpha * (N / (N2 - 1.0 ));
    	real_t dD_dne = - dalpha_dne * (N / (N2 - 1.0 ));
    	
    	real_t denom = (N2 - 1.0 - alpha * N2);
    	real_t denom2 = denom*denom;
    	real_t Q_sqrt = sqrt(alpha2 * (1.0 + n_par4 + (4.0 * N2 - 2.0) * n_par2) - 4.0 * alpha3 * N2 * n_par2);
    	real_t Q = - Q_sqrt / (2.0 * denom);
    	real_t dQ_dne_1 = alpha * (N2 - 1) * (1 + n_par4 + (4 * N2 - 2.0) * n_par2);
    	real_t dQ_dne_2 = - 6 * alpha2 * N2 * (N2 - 1) * n_par2;
    	real_t dQ_dne_3 = 2 * alpha3 * N4 * n_par2;
    	real_t dQ_dne = - 1.0 / 2.0 * dalpha_dne * (dQ_dne_1 + dQ_dne_2 + dQ_dne_3) / (denom2 * Q_sqrt);
    	
    	real_t R_1 = alpha * (N2 - 1.0) * (2.0 * N - 1.0 + n_par2);
    	real_t R_2 = - 2.0 * alpha2 * N2 * (N - 1.0);
    	real_t R = (R_1 + R_2) / (2.0 * (N2 - 1.0) * denom);
    	real_t dR_dne_0 = (N2 - 1.0)*(N2 - 1.0) * (2 * N - 1 + n_par2);
    	real_t dR_dne_1 = - 4.0 * alpha * N2 * (N2 - 1.0) * (N - 1.0);
    	real_t dR_dne_2 = 2.0 * alpha2 * N4 * (N - 1);
    	real_t dR_dne = 1.0 / (2.0 * (N2 - 1.0)) * dalpha_dne * (dR_dne_0 + dR_dne_1 + dR_dne_2) / (denom2);
    	
    	real_t T_0 = 2.0 * (N2 - 1.0) * n_par2;
    	real_t T_1 = alpha * (1.0 - (2.0 * N2 - 1.0) * n_par2);
    	real_t T = (T_0 + T_1) / (2.0 * denom);
    	real_t dT_dne = 1.0 / 2.0 * dalpha_dne * (N2 - 1.0) * (1.0 + n_par2) / denom2;
    	
    	real_t V_1 = - alpha * (N2 - 1.0) * (1.0 - n_par2);
    	real_t V_2 = 2.0 * alpha2 * N2;
    	real_t V = (V_1 + V_2) / (2.0 * (N2 - 1.0) * denom);
    	real_t dV_dne_0 = - (N2 - 1.0)*(N2 - 1.0) * (1.0 - n_par2);
    	real_t dV_dne_1 = 4.0 * alpha * N2 * (N2 - 1.0);
    	real_t dV_dne_2 = - 2.0 * alpha2 * N4; 
    	real_t dV_dne = 1.0 / (2.0 * (N2 - 1.0)) * dalpha_dne * (dV_dne_0 + dV_dne_1 + dV_dne_2) / denom2;
    	
    	real_t n_perp_0 = 2.0 * (N2 - 1) * (1.0 - n_par2);
    	real_t n_perp_1 = - alpha * (4.0 * N2 - 1 - (2.0 * N2 - 1.0) * n_par2);
    	real_t n_perp_2 = 2 * alpha2 * N2;
    	real_t n_perp2 = (n_perp_0 + n_perp_1 + n_perp_2) / (2.0 * denom) - Q;
    	real_t n_perp = sqrt(n_perp2);
    	real_t dnperp_dne_0 = - (N2 - 1.0) * (2.0 * N2 - 1.0 + n_par2);
    	real_t dnperp_dne_1 = 4.0 * alpha * N2 * (N2 - 1.0);
    	real_t dnperp_dne_2 = - 2.0 * alpha2 * N4;
    	real_t dnperp2_dne = 1.0 / 2.0 * dalpha_dne * (dnperp_dne_0 + dnperp_dne_1 + dnperp_dne_2) / denom2 - dQ_dne;
    	real_t dnperp_dne = 1.0 / (2.0 * n_perp) * dnperp2_dne;
    	
    	real_t F_1 = (R + Q) * (R + Q) * (T + Q) * (T + Q);
    	real_t dF1_dne = 2.0 * (R + Q) * (T + Q) * ((T + Q) * (dR_dne + dQ_dne)  + (R + Q) * (dT_dne + dQ_dne));
    	real_t F_2 = (D*D) * ((T + Q)*(T + Q)) + n_par2 * P *((V + Q)*(V + Q)); 
    	real_t dF2_dne = 2.0 * D * (T + Q) * ((T + Q) * dD_dne + D * (dT_dne + dQ_dne)) + n_par2 * (V + Q) * ((V + Q) * dP_dne + 2 * P * (dV_dne + dQ_dne));
    	
    	detao_dne = G * pow((T_e / mc2), N - 1) * pow(n_perp, 2 * N - 4) * ((2.0 * N - 3.0) * alpha * F_1 / F_2 * dnperp_dne + n_perp * (dalpha_dne * F_1 / F_2 + alpha * 1.0 / F_2 * dF1_dne - alpha * F_1 / (F_2*F_2) * dF2_dne));
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
    real_t mc2  = DREAM::Constants::mc2inEV;
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * (omega*omega) / (ec*ec);
    real_t n_par = cos(phi);
    real_t n_par2 = n_par*n_par;
    real_t n_par4 = n_par2*n_par2;
    len_t N2 = N*N;
    len_t N3 = N2*N;
    len_t N4 = N2*N2;
    real_t alpha = n_e / n_c;
    real_t dalpha_dne = 1 / n_c;
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detax_dne;
    
    if (N == 1) {
        real_t dterm1 = - 1.5 * sqrt(2.0 - alpha) * ((1.0 + alpha)*(1.0 + alpha)) / alpha;
        real_t dterm2 = 2 * pow(2.0 - alpha, 1.5) * (1.0 + alpha) / alpha;
        real_t dterm3 = - pow(2.0 - alpha, 1.5) * ((1.0 + alpha)*(1.0 + alpha)) / (alpha*alpha);
	detax_dne = A * (T_e / mc2) * n_par2 * (dterm1 + dterm2 + dterm3);
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2);

    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	
    	real_t P = 1.0 - alpha;
    	real_t dP_dne = - dalpha_dne;
    	
    	real_t D = - alpha * (N / (N2 - 1.0 ));
    	real_t dD_dne = - dalpha_dne * (N / (N2 - 1.0 ));
    	
    	real_t denom = (N2 - 1.0 - alpha * N2);
    	real_t denom2 = denom*denom;
    	real_t Q_sqrt = sqrt(alpha2 * (1.0 + n_par4 + (4.0 * N2 - 2.0) * n_par2) - 4.0 * alpha3 * N2 * n_par2);
    	real_t Q = Q_sqrt / (2.0 * denom);
    	real_t dQ_dne_1 = alpha * (N2 - 1) * (1 + n_par4 + (4 * N2 - 2.0) * n_par2);
    	real_t dQ_dne_2 = - 6 * alpha2 * N2 * (N2 - 1) * n_par2;
    	real_t dQ_dne_3 = 2 * alpha3 * N4 * n_par2;
    	real_t dQ_dne = 1.0 / 2.0 * dalpha_dne * (dQ_dne_1 + dQ_dne_2 + dQ_dne_3) / (denom2 * Q_sqrt);
    	
    	real_t R_1 = alpha * (N2 - 1.0) * (2.0 * N - 1.0 + n_par2);
    	real_t R_2 = - 2.0 * alpha2 * N2 * (N - 1.0);
    	real_t R = (R_1 + R_2) / (2.0 * (N2 - 1.0) * denom);
    	real_t dR_dne_0 = (N2 - 1.0)*(N2 - 1.0) * (2 * N - 1 + n_par2);
    	real_t dR_dne_1 = - 4.0 * alpha * N2 * (N2 - 1.0) * (N - 1.0);
    	real_t dR_dne_2 = 2.0 * alpha2 * N4 * (N - 1);
    	real_t dR_dne = 1.0 / (2.0 * (N2 - 1.0)) * dalpha_dne * (dR_dne_0 + dR_dne_1 + dR_dne_2) / (denom2);
    	
    	real_t T_0 = 2.0 * (N2 - 1.0) * n_par2;
    	real_t T_1 = alpha * (1.0 - (2.0 * N2 - 1.0) * n_par2);
    	real_t T = (T_0 + T_1) / (2.0 * denom);
    	real_t dT_dne = 1.0 / 2.0 * dalpha_dne * (N2 - 1.0) * (1.0 + n_par2) / denom2;
    	
    	real_t V_1 = - alpha * (N2 - 1.0) * (1.0 - n_par2);
    	real_t V_2 = 2.0 * alpha2 * N2;
    	real_t V = (V_1 + V_2) / (2.0 * (N2 - 1.0) * denom);
    	real_t dV_dne_0 = - (N2 - 1.0)*(N2 - 1.0) * (1.0 - n_par2);
    	real_t dV_dne_1 = 4.0 * alpha * N2 * (N2 - 1.0);
    	real_t dV_dne_2 = - 2.0 * alpha2 * N4; 
    	real_t dV_dne = 1.0 / (2.0 * (N2 - 1.0)) * dalpha_dne * (dV_dne_0 + dV_dne_1 + dV_dne_2) / denom2;
    	
    	real_t n_perp_0 = 2.0 * (N2 - 1) * (1.0 - n_par2);
    	real_t n_perp_1 = - alpha * (4.0 * N2 - 1 - (2.0 * N2 - 1.0) * n_par2);
    	real_t n_perp_2 = 2 * alpha2 * N2;
    	real_t n_perp2 = (n_perp_0 + n_perp_1 + n_perp_2) / (2.0 * denom) - Q;
    	real_t n_perp = sqrt(n_perp2);
    	real_t dnperp_dne_0 = - (N2 - 1.0) * (2.0 * N2 - 1.0 + n_par2);
    	real_t dnperp_dne_1 = 4.0 * alpha * N2 * (N2 - 1.0);
    	real_t dnperp_dne_2 = - 2.0 * alpha2 * N4;
    	real_t dnperp2_dne = 1.0 / 2.0 * dalpha_dne * (dnperp_dne_0 + dnperp_dne_1 + dnperp_dne_2) / denom2 - dQ_dne;
    	real_t dnperp_dne = 1.0 / (2.0 * n_perp) * dnperp2_dne;
    	
    	real_t F_1 = (R + Q) * (R + Q) * (T + Q) * (T + Q);
    	real_t dF1_dne = 2.0 * (R + Q) * (T + Q) * ((T + Q) * (dR_dne + dQ_dne)  + (R + Q) * (dT_dne + dQ_dne));
    	real_t F_2 = (D*D) * ((T + Q)*(T + Q)) + n_par2 * P *((V + Q)*(V + Q)); 
    	real_t dF2_dne = 2.0 * D * (T + Q) * ((T + Q) * dD_dne + D * (dT_dne + dQ_dne)) + n_par2 * (V + Q) * ((V + Q) * dP_dne + 2 * P * (dV_dne + dQ_dne));
    	
    	detax_dne = G * pow((T_e / mc2), N - 1) * pow(n_perp, 2 * N - 4) * ((2.0 * N - 3.0) * alpha * F_1 / F_2 * dnperp_dne + n_perp * (dalpha_dne * F_1 / F_2 + alpha * 1.0 / F_2 * dF1_dne - alpha * F_1 / (F_2*F_2) * dF2_dne));
    }
    
    return detax_dne;
}

