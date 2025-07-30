/**
 * Class for calculating plasma volume and neutral volume
 */
 
#include "STREAM/Equations/PlasmaVolume.hpp"
#include <cmath> 
 
using namespace STREAM;
using namespace DREAM;

    PlasmaVolume::PlasmaVolume(FVM::Grid *g, real_t vessel_vol, FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, ADAS *adas, IonHandler *ions): grid(g), vessel_vol(vessel_vol), unknowns(u), radials(r), adas(adas), ions(ions){
        id_lambda_i = u->GetUnknownID(OptionConstants::UQTY_LAMBDA_I);
        id_W_i = u->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
        id_n_i = u->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
        id_T_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD); 
        id_n_cold = u->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD); 
    }
    
	real_t PlasmaVolume::GetPlasmaCrossSection() const {
		return grid->GetVpVol(0) * radials->GetMinorRadius() / (2*M_PI);
	}

    real_t PlasmaVolume::GetPlasmaVolume() const{
        return grid->GetVpVol(0) * radials->GetMinorRadius() * radials->GetMajorRadius(); 
    }
    
    real_t PlasmaVolume::GetNeutralVolume(const len_t iz){ 
        real_t a = radials->GetMinorRadius();       
        real_t lambda_i = unknowns->GetUnknownData(id_lambda_i)[iz];  

        //Eq. 13 in startup-appendix + added triangularity 
        if (lambda_i <= a){
            real_t R0 = radials->GetMajorRadius();
            real_t kappa = radials->GetElongation(); 
            //real_t delta = radials->GetTriangularity();  
            
            return 2.0 * M_PI*M_PI * R0 * kappa * (a*a - (a-lambda_i)*(a-lambda_i));
            /*return  2.0*M_PI*M_PI*R0*kappa*(a*a-(a-lambda_i)*(a-lambda_i)) +
                    2.0*(8.0-3.0*M_PI*M_PI) *kappa *delta*(a*a*a - (a-lambda_i)*(a-lambda_i)*(a-lambda_i))/3.0;*/
        } else {
            return GetPlasmaVolume();
        }
    }
    
    //gamma_n,i * vessel volume V, with gamma_n,i from eq. 12 in startup-appendix
    real_t PlasmaVolume::GetTotalNeutralVolume(const len_t iz){
        return vessel_vol - GetPlasmaVolume() + GetNeutralVolume(iz);    
    }
    
    /*real_t PlasmaVolume::GetNeutralVolume_dT(const len_t iz){
        real_t a = radials->GetMinorRadius(); 
        real_t lambda_i = unknowns->GetUnknownData(id_lambda_i)[iz];

        if (lambda_i <= a){
            real_t R0 = radials->GetMajorRadius();
            real_t kappa = radials->GetElongation(); 
            real_t delta = radials->GetTriangularity(); 
            real_t W_i = unknowns->GetUnknownData(id_W_i)[iz];
            real_t n_i = unknowns->GetUnknownData(id_n_i)[iz];
            real_t T_cold = unknowns->GetUnknownData(id_T_cold)[0];
            real_t n_cold = unknowns->GetUnknownData(id_n_cold)[0];
        
            len_t Z = ions->GetZ(iz);
            real_t I_i = adas->GetSCD(Z)->Eval(0, n_cold, T_cold); //Evaluate I_i^(0)
            real_t dIdT = adas->GetSCD(Z)->Eval_deriv_T(0, n_cold, T_cold); // Derivative w.r.t. T_cold

            real_t v_i = 0;
            if(n_i != 0) 
                v_i = sqrt(4.0 * W_i/(3.0 * n_i * ions->GetIonSpeciesMass(iz)));
            
            return  (4.0*M_PI*M_PI*R0*kappa*(a-lambda_i)+2.0*(8.0-3.0*M_PI*M_PI)*kappa*delta*(a-lambda_i)*(a-lambda_i)) * -v_i/(n_cold*I_i*I_i) * dIdT; 
        } else {
            return 0;
        }
    }
    
    real_t PlasmaVolume::GetNeutralVolume_dn(const len_t iz){              
        real_t a = radials->GetMinorRadius(); 
        real_t lambda_i = unknowns->GetUnknownData(id_lambda_i)[iz];
        
        if (lambda_i <= a){
            real_t R0 = radials->GetMajorRadius();
            real_t kappa = radials->GetElongation(); 
            real_t delta = radials->GetTriangularity(); 
            real_t W_i = unknowns->GetUnknownData(id_W_i)[iz];
            real_t n_i = unknowns->GetUnknownData(id_n_i)[iz];
            real_t T_cold = unknowns->GetUnknownData(id_T_cold)[0];
            real_t n_cold = unknowns->GetUnknownData(id_n_cold)[0];
            
            len_t Z = ions->GetZ(iz);
            real_t I_i = adas->GetSCD(Z)->Eval(0, n_cold, T_cold); //Evaluate I_i^(0)
            real_t dIdn = adas->GetSCD(Z)->Eval_deriv_n(0, n_cold, T_cold); // Derivative w.r.t. n_cold
            
            real_t v_i = 0;
            if(n_i != 0)
                v_i = sqrt(4.0 * W_i/(3.0 * n_i * ions->GetIonSpeciesMass(iz)));
        
            return  (4.0*M_PI*M_PI*R0*kappa*(a-lambda_i)+2.0*(8.0-3.0*M_PI*M_PI)*kappa*delta*(a-lambda_i)*(a-lambda_i)) * -v_i/(n_cold*I_i) * (1.0/n_cold + dIdn/I_i); 
        } else {
            return 0;
        }
    }*/

    real_t PlasmaVolume::GetNeutralVolume_dLambdai(const len_t iz) {
        real_t a = radials->GetMinorRadius();
        real_t lambda_i = unknowns->GetUnknownData(id_lambda_i)[iz];

        if (lambda_i <= a){
            real_t R0 = radials->GetMajorRadius();
            real_t kappa = radials->GetElongation(); 
            //real_t delta = radials->GetTriangularity(); 
            //real_t a_l = a-lambda_i;
	    
	    return 4.0 * M_PI*M_PI * R0 * kappa * (a-lambda_i);
            //return (4.0*M_PI*M_PI*R0*kappa*a_l + 2.0*kappa*delta*(8.0-3.0*M_PI*M_PI)*a_l*a_l);
        } else
            return 0;
    }

    real_t PlasmaVolume::GetTotalNeutralVolume_dLambdai(const len_t iz) {
        return GetNeutralVolume_dLambdai(iz);
    }
    
    /*real_t PlasmaVolume::GetTotalNeutralVolume_dT(const len_t iz){
        return GetNeutralVolume_dT(iz);    
    }
    
    real_t PlasmaVolume::GetTotalNeutralVolume_dn(const len_t iz){
        return GetNeutralVolume_dn(iz);    
    }*/
