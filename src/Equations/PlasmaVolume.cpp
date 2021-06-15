/**
 * Class for calculating plasma volume and neutral volume
 */
 
#include "STREAM/Equations/PlasmaVolume.hpp"
#include <cmath> 
 
using namespace STREAM;
using namespace DREAM;

    PlasmaVolume::PlasmaVolume(FVM::Grid *g, len_t iz, real_t vessel_vol, FVM::UnknownQuantityHandler *u): grid(g), iz(iz), vessel_vol(vessel_vol), unknowns(u){
        id_lambda_i = u->GetUnknownID(OptionConstants::UQTY_); //Needs to be added as an unknown quantity
    }
    
    real_t PlasmaVolume::GetPlasmaVolume() const{
        return g->GetVpVol(0) * g->GetRadialGrid()->GetDr(0); //ok?
    }
    
    real_t PlasmaVolume::GetNeutralVolume(const len_t iz, real_t t){
        real_t a = grid->GetMinorRadius()->Eval(t); //Evaluate this at t or not?
        real_t R0 = grid->GetR0()->Eval(t); //Evaluate this at t or not?
        real_t lambda_i = unknowns->GetUnknownData(id_lambda_i)[0]; //How do we get lambda with correct index iz? We need tou use input argument iz somehow
        
        //TODO: Need to get these from somewhere or evaluate in terms of other available geometric properties? Are b and R_edge available?
        real_t kappa //=b/a, where b is vertical radius?
        real_t delta //=(R0- R_edge)/a, where R_edge is the major radius at the top/bottom of the non-elliptical cross-sectional area?        
        
        //Eq. 13 in startup-appendix
        if (lambda_i <= a){
            return  2*M_PI*M_PI*R0*kappa*(a*a-(a-lambda_i)*(a-lambda_i));
        } else {
            return GetPlasmaVolume();
        }
         
         /**
         
         //Attempt to include delta
        if (lambda_i <= a){
            return  2*M_PI*M_PI*R0*kappa*(a*a-(a-lambda_i)*(a-lambda_i)) + 2*(8-3*M_PI*M_PI) *kappa *delta*(a*a*a - (a-lambda_i)*(a-lambda_i)*(a-lambda_i))/3;
        } else {
            return GetPlasmaVolume();
        }
        
        */
         
    }
