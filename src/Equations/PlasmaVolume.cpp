/**
 * Class for calculating plasma volume and neutral volume
 */
 
#include "STREAM/Equations/PlasmaVolume.hpp"
#include <cmath> 
 
using namespace STREAM;
using namespace DREAM;

    PlasmaVolume::PlasmaVolume(FVM::Grid *g, len_t iz, real_t vessel_vol, FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r): grid(g), iz(iz), vessel_vol(vessel_vol), unknowns(u), radials(r){
        id_lambda_i = u->GetUnknownID(OptionConstants::UQTY_LAMBDA_I); //Needs to be added as an unknown quantity
    }
    
    real_t PlasmaVolume::GetPlasmaVolume() const{
        return grid->GetVpVol(0) * grid->GetRadialGrid()->GetDr(0); //ok?
    }
    
    real_t PlasmaVolume::GetNeutralVolume(const len_t iz){ 
        real_t a = radials->GetMinorRadius(); 
        real_t R0 = grid->GetR0();
        real_t kappa = radials->GetElongation(); 
        real_t delta = radials->GetTriangularity(); 
        real_t lambda_i = unknowns->GetUnknownData(id_lambda_i)[0]; //How do we get lambda with correct index iz? Don't we need to use input argument iz somehow
         
        //Eq. 13 in startup-appendix + added triangularity 
        if (lambda_i <= a){
            return  2*M_PI*M_PI*R0*kappa*(a*a-(a-lambda_i)*(a-lambda_i)) + 2*(8-3*M_PI*M_PI) *kappa *delta*(a*a*a - (a-lambda_i)*(a-lambda_i)*(a-lambda_i))/3;
        } else {
            return GetPlasmaVolume();
        }
        
        /**
        //Eq. 13 in startup-appendix
        if (lambda_i <= a){
            return  2*M_PI*M_PI*R0*kappa*(a*a-(a-lambda_i)*(a-lambda_i));
        } else {
            return GetPlasmaVolume();
        }
        */ 
    }
