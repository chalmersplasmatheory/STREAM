#include "STREAM/Equations/OpticalThickness.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

// TODO : Add derivatives
 
/**
 * Constructor
 */
OpticalThickness::OpticalThickness(
	FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r
) {
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
    
    real_t ec   = DREAM::Constants::ec;
    real_t c    = DREAM::Constants::c;
    real_t me   = DREAM::Constants::me;
    real_t eps0 = DREAM::Constants::eps0;
    
    real_t R0 = radials->GetMajorRadius();
    real_t B = radials->GetMagneticField();
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2 * M_PI * c / omega;
    real_t n_c = m_e * eps0 * omega*omega / (ec*ec);
    real_t n_parallell = cos(phi);
    real_t alpha = n_e / n_c;
    real_t global_factor = (M_PI*M_PI * R0 / lambda) * (T_e / me);
    
    real_t eta_o;
    if (N == 1) {
	eta_o = global_factor * alpha * sqrt(1 - theta) * 1 / (1 + n_parallell*n_parallell * (0.5 - theta)); // (n_e / n_c) = alpha
    } else {
    	real_t S = 1 - alpha * (N*N / (N*N - 1));
    	real_t P = 1 - alpha;
    	real_t D = - alpha * (N / (N*N - 1));
    	
    	real_t A = S;
    	real_t B = - (S + P) * (S - n_parallell*n_parallell) + D*D;
    	real_t C = P * (pow(S - n_parallell*n_parallell,2) - D*D);
    	
    	real_t n_perpo = sqrt((- B + sqrt(B*B - 4 * A * C)) / (2 * A));
    	real_t n_o = sqrt(n_perpo*n_perpo + n_parallell*n_parallell);
        real_t factor1o = pow(n_perpo*n_perpo * N*N * T_e / (2 * me * c*c), N-2);
        real_t factor2o = pow((S - D - n_o*n_o) * (P - n_perpo*n_perpo), 2) / (D*D * pow(P - n_perpo*n_perpo, 2) + n_parallell*n_parallell * P * pow(S - n_perpo*n_perpo, 2));
        
        eta_o = global_factor * 1 / (c*c) * pow(N, 3) / tgamma(N) * theta * n_perpo * factor1o * factor2o;
    }
    this->eta_polo = eta_o / cos(theta);
    this->eta_polx = eta_x / cos(theta);
    
    return eta_polo;
}

/**
 * Evaluates the optical thickness of the X mode 
 */
real_t OpticalThickness::EvaluateOpticalThickness_o(len_t ir){ 
    real_t T_e = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t n_e = unknowns->GetUnknownData(id_ncold)[ir];
    
    real_t ec   = DREAM::Constants::ec;
    real_t c    = DREAM::Constants::c;
    real_t me   = DREAM::Constants::me;
    real_t eps0 = DREAM::Constants::eps0;
    
    real_t R0 = radials->GetMajorRadius();
    real_t B = radials->GetMagneticField();
    
    real_t omega_ce = ec * B / me;
    real_t omega = N * omega_ce;
    real_t lambda = 2 * M_PI * c / omega;
    real_t n_c = m_e * eps0 * omega*omega / (ec*ec);
    real_t n_parallell = cos(phi);
    real_t alpha = n_e / n_c;
    real_t global_factor = (M_PI*M_PI * R0 / lambda) * (T_e / me);
    
    real_t eta_o;
    if (N == 1) {
	eta_x = global_factor * n_parallell*n_parallell * pow(2 - theta, 1.5) * pow(1 + theta, 2) / theta;
    } else {
    	real_t S = 1 - alpha * (N*N / (N*N - 1));
    	real_t P = 1 - alpha;
    	real_t D = - alpha * (N / (N*N - 1));
    	
    	real_t A = S;
    	real_t B = - (S + P) * (S - n_parallell*n_parallell) + D*D;
    	real_t C = P * (pow(S - n_parallell*n_parallell,2) - D*D);
    	
    	real_t n_perpx = sqrt((- B - sqrt(B*B - 4 * A * C)) / (2 * A));
    	real_t n_x = sqrt(n_perpx*n_perpx + n_parallell*n_parallell);
        real_t factor1x = pow(n_perpx*n_perpx * N*N * T_e / (2 * me * c*c), N-2);
        real_t factor2x = pow((S - D - n_x*n_x) * (P - n_perpx*n_perpx), 2) / (D*D * pow(P - n_perpx*n_perpx, 2) + n_parallell*n_parallell * P * pow(S - n_perpx*n_perpx, 2));
        
        eta_x = global_factor * 1 / (c*c) * pow(N, 3) / tgamma(N) * theta * n_perpx * factor1o * factor2o;
    }
    this->eta_polx = eta_x / cos(theta);
    
    return eta_polo;
}
