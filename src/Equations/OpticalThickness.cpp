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
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t eta_o;
    
    if (N == 1) {
	eta_o = A * (T_e / mc2) * n_e / n_c * sqrt(1.0 - theta) / ( 1.0 + n_par2 * ( 0.5 - theta ) );
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2) * theta;

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
    	
    	eta_o = G * pow((T_e / mc2), N - 1.0) * pow(n_perp, 2 * N - 3) * F_1 / F_2;
    	//printf("\neta_x=%.4e", eta_o);
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
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t eta_x;
    
    if (N == 1) {
	eta_x = A * (T_e / mc2) * n_par2 * pow(2.0 - theta, 1.5) * ((1.0 + theta)*(1.0 + theta)) / theta;
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2) * theta;

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
    	
    	eta_x = G * pow((T_e / mc2), N - 1.0) * pow(n_perp, 2 * N - 3) * F_1 / F_2;
    	//printf("\neta_x=%.4e", eta_x);
    	//printf("\n For x:     eta_x=%.3e, n_perp=%.3e, F_1=%.3e, F_2=%.3e, R=%.3e, Q=+%.3e", eta_x, n_perp, F_1, F_2, R, Q);
	real_t N_par2 = sin(20 * M_PI / 180) * 1.8 / 5.4;
    	real_t tau = 6.6e-22 * n_e * T_e * pow(1 - N_par2, 1.5) / cos(14 * M_PI / 180);
    	//printf("\n SCENPLINT:   tau=%.3e\n", tau);
    	//printf("\n 1/n_c*mc2=%.3e\n", 1/(n_c*mc2));
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
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detao_dTe;
    
    if (N == 1) {
	detao_dTe = A * (1.0 / mc2) * n_e / n_c * sqrt(1.0 - theta) / ( 1.0 + n_par2 * ( 0.5 - theta ) );
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2) * theta;

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
    	
    	detao_dTe = G * (N - 1.0) / mc2 * pow((T_e / mc2), N - 2) * pow(n_perp, 2 * N - 3) * F_1 / F_2;
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
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detax_dTe;
    
    if (N == 1) {
	detax_dTe = A * 1.0 / mc2 * n_par2 * pow(2.0 - theta, 1.5) * ((1.0 + theta)*(1.0 + theta)) / theta;
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2) * theta;

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
    	
    	detax_dTe = G * (N - 1.0) / mc2 * pow((T_e / mc2), N - 2) * pow(n_perp, 2 * N - 3) * F_1 / F_2;
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
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detao_dne;
    
    if (N == 1) {
	detao_dne = A * (T_e / mc2) * 1 / n_c * sqrt(1.0 - theta) / ( 1.0 + n_par2 * ( 0.5 - theta ) );
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2) * theta;

    	real_t alpha = n_e / n_c;
    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	real_t dalpha_dne = 1 / n_c;
    	
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
    	
    	detao_dne = G * pow((T_e / mc2), N - 1) * pow(n_perp, 2 * N - 4) * ((2.0 * N - 3.0) * F_1 / F_2 * dnperp_dne + n_perp * (1.0 / F_2 * dF1_dne - F_1 / (F_2*F_2) * dF2_dne));
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
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detax_dne;
    
    if (N == 1) {
	detax_dne = 0.0;
    } else {
    	real_t G = A * N3 / tgamma(N) * pow(N2 / 2.0, N - 2) * theta;

    	real_t alpha = n_e / n_c;
    	real_t alpha2 = alpha*alpha;
    	real_t alpha3 = alpha2*alpha;
    	real_t dalpha_dne = 1 / n_c;
    	
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
    	
    	detax_dne = G * pow((T_e / mc2), N - 1) * pow(n_perp, 2 * N - 4) * ((2.0 * N - 3.0) * F_1 / F_2 * dnperp_dne + n_perp * (1.0 / F_2 * dF1_dne - F_1 / (F_2*F_2) * dF2_dne));
    }
    
    return detax_dne;
}


/*
void OpticalThickness::Initialize() {
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
}*/

/*

/**
 * Evaluates the optical thickness of the O mode 
 * /
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
    real_t omega = ((real_t) N) * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t eta_o;
    
    if (N == 1) {
	eta_o = A * (T_e / mc2) * n_e / n_c * sqrt(1.0 - theta) / ( 1.0 + pow(n_parallel, 2) * ( 0.5 - theta ) );
    } else {
    	real_t G = A * ((real_t) pow(N,3)) / ((real_t)tgamma(N)) * pow( ((real_t)pow(N,2)) / 2.0, N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	//printf("\nS=%.4e", S);
    	real_t P = 1.0 - n_e / n_c;
    	//printf("\nP=%.4e", P);
    	real_t D = - n_e / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	//printf("\nD=%.4e", D);
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpo2 = (- E + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	//printf("\nn_perpo2=%.4e", n_perpo2);
    	real_t n_perpo  = sqrt(n_perpo2);
    	real_t n_o2     = n_perpo2 + pow(n_parallel,2);
    	//printf("\nn_o2=%.4e", n_o2);
    	
    	printf("\nS=%.4e, D=%.4e, P=%.4e, n_o2=%.4e, n_perpo2=%.4e", S, D, P, n_o2, n_perpo2);
    	//real_t F_1o = pow( (S - D - n_o2) * (P - n_perpo2) , 2);
    	real_t SmP = - n_e / n_c * 1 / ( ((real_t)pow(N,2)) - 1.0 );
    	real_t fac1 = (SmP - pow(n_parallel,2)) / 2.0 - D + (P * pow(n_parallel,2) + pow(D,2) - sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t F_1o = pow( fac1 * (P - n_perpo2) , 2);
    	printf(", F_1o=%.4e", F_1o);
    	real_t F_2o = pow(D,2) * pow( P - n_perpo2 , 2) + pow(n_parallel,2) * P * pow( S - n_o2 , 2);
    	//printf("\nF_2o=%.4e", F_2o);
    	
    	eta_o = G * pow((T_e / mc2), N - 1) * pow(n_perpo, 2 * N - 3) * F_1o / F_2o;
    	//printf("\neta_o=%.4e", eta_o);
    }
    
    return eta_o;
}

/**
 * Evaluates the optical thickness of the X mode 
 * /
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
    real_t omega = ((real_t) N) * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t eta_x;
    
    if (N == 1) {
	eta_x = A * (T_e / mc2) * pow(n_parallel,2) * pow(2.0 - theta, 1.5) * pow(1 + theta, 2) / theta;
    } else {
    	real_t G = A * ((real_t) pow(N,3)) / ((real_t)tgamma(N)) * pow( ((real_t)pow(N,2)) / 2.0, N - 2) * theta;
    	//printf("G=%.4e\n", G);
    	
    	real_t S = 1.0 - n_e / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	//printf("\nS=%.4e", S);
    	real_t P = 1.0 - n_e / n_c;
    	//printf("\nP=%.4e", P);
    	real_t D = - n_e / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	//printf("\nD=%.4e", D);
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	//printf("E=%.4e\n", E);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	//printf("C=%.4e\n", C);
    	
    	real_t n_perpx2 = (- E - sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	//printf("\nn_perpx2=%.4e", n_perpx2);
    	real_t n_perpx  = sqrt(n_perpx2);
    	real_t n_x2     = n_perpx2 + pow(n_parallel,2);
    	//printf("\nn_x2=%.4e", n_x2);
    	
    	printf("\nS=%.4e, D=%.4e, P=%.4e, n_x2=%.4e, n_perpx2=%.4e", S, D, P, n_x2, n_perpx2);
    	//real_t F_1x = pow( (S - D - n_x2) * (P - n_perpx2) , 2);
    	real_t SmP = - n_e / n_c * 1 / ( ((real_t)pow(N,2)) - 1.0 );
    	real_t fac1 = (SmP - pow(n_parallel,2)) / 2.0 - D + (P * pow(n_parallel,2) + pow(D,2) + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t F_1x = pow( fac1 * (P - n_perpx2) , 2);
    	printf(", F_1x=%.4e", F_1x);
    	real_t F_2x = pow(D,2) * pow( P - n_perpx2 , 2) + pow(n_parallel,2) * P * pow( S - n_x2 , 2);
    	//printf("\nF_2x=%.4e", F_2x);
    	
    	eta_x = G * pow((T_e / mc2), N - 1) * pow(n_perpx, 2 * N - 3) * F_1x / F_2x;
    	//printf("\neta_x=%.4e", eta_x);
    }
    
    return eta_x;
}

/**
 * Evaluates the derivative of the optical thickness of the O mode with regards to the electron temperature
 * /
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
    real_t omega = ((real_t) N) * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detao_dTe;
    
    if (N == 1) {
	detao_dTe = A * 1 / mc2 * n_e / n_c * sqrt(1.0 - theta) / ( 1.0 + pow(n_parallel, 2) * ( 0.5 - theta ) );
    } else {
    	real_t G = A * ((real_t) pow(N,3)) / ((real_t)tgamma(N)) * pow( ((real_t)pow(N,2)) / 2.0, N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpo2 = (- E + sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpo  = sqrt(n_perpo2);
    	real_t n_o2     = n_perpo2 + pow(n_parallel,2);
    	
    	real_t F_1o = pow( (S - D - n_o2) * (P - n_perpo2) , 2);
    	real_t F_2o = pow(D,2) * pow( P - n_perpo2 , 2) + pow(n_parallel,2) * P * pow( S - n_o2 , 2);
    	
    	detao_dTe = G * (N - 1.0) * 1 / mc2 * pow((T_e / mc2), N - 2) * pow(n_perpo, 2 * N - 3) * F_1o / F_2o;
    }
    
    return detao_dTe;
}

/**
 * Evaluates the derivative of the optical thickness of the O mode with regards to the electron temperature
 * /
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
    real_t omega = ((real_t) N) * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detax_dTe;
    
    if (N == 1) {
	detax_dTe = A * 1 / mc2 * pow(n_parallel,2) * pow(2.0 - theta, 1.5) * pow(1 + theta, 2) / theta;
    } else {
    	real_t G = A * ((real_t) pow(N,3)) / ((real_t)tgamma(N)) * pow( ((real_t)pow(N,2)) / 2.0, N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t n_perpx2 = (- E - sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpx  = sqrt(n_perpx2);
    	real_t n_x2     = n_perpx2 + pow(n_parallel,2);
    	
    	real_t F_1x = pow( (S - D - n_x2) * (P - n_perpx2) , 2);
    	real_t F_2x = pow(D,2) * pow( P - n_perpx2 , 2) + pow(n_parallel,2) * P * pow( S - n_x2 , 2);
    	
    	detax_dTe = G * (N - 1.0) * 1 / mc2 * pow((T_e / mc2), N - 2) * pow(n_perpx, 2 * N - 3) * F_1x / F_2x;
    }
    
    return detax_dTe;
}

/**
 * Evaluates the derivative of the optical thickness of the O mode with regards to the electron density
 * /
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
    real_t omega = ((real_t) N) * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detao_dne;
    
    if (N == 1) {
	detao_dne = A * (T_e / mc2) * 1.0 / n_c * sqrt(1.0 - theta) / ( 1.0 + pow(n_parallel, 2) * ( 0.5 - theta ) );
    } else {
    	real_t G = A * ((real_t) pow(N,3)) / ((real_t)tgamma(N)) * pow( ((real_t)pow(N,2)) / 2.0, N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	
    	real_t dS_dne = - 1.0 / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	real_t dP_dne = - 1.0 / n_c;
    	real_t dD_dne = - 1.0 / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	
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
    	
    	detao_dne = G * pow((T_e / mc2), N - 1) * pow(n_perpo, 2 * N - 4) * ( (2.0 * N - 3.0) * F_1o / F_2o * dnperpo_dne + n_perpo * (1.0 / F_2o * dF1o_dne - F_1o / pow(F_2o,2) * dF2o_dne) );
    }
    
    return detao_dne;
}

/**
 * Evaluates the derivative of the optical thickness of the X mode with regards to the electron density
 * /
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
    real_t omega = ((real_t) N) * omega_ce;
    real_t lambda = 2.0 * M_PI * c / omega;
    real_t n_c = me * eps0 * pow(omega,2) / pow(ec,2);
    real_t n_parallel = cos(phi);
    
    real_t A = (M_PI*M_PI) * R0 / lambda;
    
    real_t detax_dne;
    
    if (N == 1) {
	detax_dne = 0.0;
    } else {
    	real_t G = A * ((real_t) pow(N,3)) / ((real_t)tgamma(N)) * pow( ((real_t)pow(N,2)) / 2.0, N - 2) * theta;
    	
    	real_t S = 1.0 - n_e / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	real_t P = 1.0 - n_e / n_c;
    	real_t D = - n_e / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	
    	real_t dS_dne = - 1.0 / n_c * ( ((real_t)pow(N,2)) / ( ((real_t)pow(N,2)) - 1.0 ) );
    	real_t dP_dne = - 1.0 / n_c;
    	real_t dD_dne = - 1.0 / n_c * ( N / ( ((real_t)pow(N,2)) - 1.0 ) );
    	
    	real_t E = - (S + P) * (S - pow(n_parallel,2)) + pow(D,2);
    	real_t C = P * ( pow(S - pow(n_parallel,2),2) - pow(D,2));
    	
    	real_t dE_dne = (pow(n_parallel,2) - 2.0 * S - P) * dS_dne + (pow(n_parallel,2) - S) * dP_dne + 2.0 * D * dD_dne;
    	real_t dC_dne = (S - pow(n_parallel,2)) * (2.0 * P * dS_dne + (S - pow(n_parallel,2)) * dP_dne) - D * (D * dP_dne + 2 * P *dD_dne);
    	
    	real_t n_perpx2 = (- E - sqrt( pow(E,2) - 4.0 * S * C)) / (2.0 * S);
    	real_t n_perpx  = sqrt(n_perpx2);
    	real_t n_x2     = n_perpx2 + pow(n_parallel,2);
    	
    	real_t dnperpx2_dne = 1.0 / sqrt(pow(E,2) - 4.0 * S * C) * ( (-2.0 * S * C + E * sqrt(pow(E,2) - 4.0 * S * C) + pow(E,2)) / (2.0 * pow(S,2)) * dS_dne - (E + sqrt(pow(E,2) - 4.0 * S * C)) / (2.0 * S) * dE_dne + dC_dne );
    	real_t dnperpx_dne  = 1.0 / (2.0 * n_perpx) * dnperpx2_dne;
    	
    	real_t F_1x = pow( (S - D - n_x2) * (P - n_perpx2) , 2);
    	real_t F_2x = pow(D,2) * pow( P - n_perpx2 , 2) + pow(n_parallel,2) * P * pow( S - n_x2 , 2);
    	
    	real_t dF1x_dne = 2.0 * (S - D - n_x2) * (P - n_perpx2) * ( (P - n_perpx2) * (dS_dne - dD_dne - dnperpx2_dne) + (S - D - n_x2) * (dP_dne - dnperpx2_dne) );
    	real_t dF2x_dne = pow(n_parallel,2) * (S - n_x2) * ( 2.0 * P * (dS_dne - dnperpx2_dne) + (S - n_x2) * dP_dne ) + 2 * D * (P - n_perpx2) * ( D * (dP_dne - dnperpx2_dne) + (P - n_perpx2) * dD_dne ); 
    	
    	detax_dne = G * pow((T_e / mc2), N - 1) * pow(n_perpx, 2 * N - 4) * ( (2.0 * N - 3.0) * F_1x / F_2x * dnperpx_dne + n_perpx * (1.0 / F_2x * dF1x_dne - F_1x / pow(F_2x,2) * dF2x_dne) );
    }
    
    return detax_dne;
}
/*
void OpticalThickness::Initialize() {
    this->id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
}* /

*/
