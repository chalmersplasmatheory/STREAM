/**
 * Implementation of an equation term for the mean free path of
 * neutrals of species with index iz
 */

#include "STREAM/Equations/MeanFreePathTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <cmath>

using namespace STREAM;
using namespace DREAM;
using namespace std; 

    MeanFreePathTerm::MeanFreePathTerm(FVM::Grid *g, len_t iz, FVM::UnknownQuantityHandler *u, ADAS *adas, IonHandler *ions): FVM::EvaluableEquationTerm(g), iz(iz), unknowns(u), adas(adas), ions(ions){
        id_W_i = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);
        id_n_i = u->GetUnknownID(OptionConstants::UQTY_NI_DENS);
        id_T_cold = u->GetUnknownID(OptionConstants::UQTY_T_COLD); 
        id_n_cold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
        
        len_t Z = ions->GetZ(iz);
        this->interp = adas->GetSCD(Z);
    }

    /**
     * Evaluate this mean-free path term, i.e. calculate
     *
     *   lambda_i = v_i / (n_i*I_i)
     * 
     * vec: Vector containing 'lambda_i' on return.
     */
    void MeanFreePathTerm::EvaluableTransform(real_t *vec) {
        this->SetVectorElements(vec, nullptr);
    }
    
    //Rebuild the equation term 
    void MeanFreePathTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*){
        real_t W_i = unknowns->GetUnknownData(id_W_i)[this->iz];
        real_t n_i = unknowns->GetUnknownData(id_n_i)[this->iz];
        real_t T_cold = unknowns->GetUnknownData(id_T_cold)[0];
        real_t n_cold = unknowns->GetUnknownData(id_n_cold)[0]; 
                
        real_t I_i = interp->Eval(0, n_cold, T_cold);    //Evaluate I_i^(0)
        
        real_t v_i = 0;
        this->lambda_i = 0;
        if(n_i != 0){
            v_i = sqrt(4*W_i/(3*n_i*ions->GetIonSpeciesMass(iz))); 
            this->lambda_i = v_i/(n_cold*I_i); 
        }
    }
    

    bool MeanFreePathTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*){
        if (derivId != id_W_i && derivId != id_n_i && derivId != id_T_cold && derivId != id_n_cold)             
            return false;
    
        real_t W_i = unknowns->GetUnknownData(id_W_i)[this->iz];
        real_t n_i = unknowns->GetUnknownData(id_n_i)[this->iz];
        real_t T_cold = unknowns->GetUnknownData(id_T_cold)[0];
        real_t n_cold = unknowns->GetUnknownData(id_n_cold)[0];
        
        real_t I_i = interp->Eval(0, n_cold, T_cold); //Evaluate I_i^(0)
        real_t dIdn = interp->Eval_deriv_n(0, n_cold, T_cold); // Derivative w.r.t. n_cold
        real_t dIdT = interp->Eval_deriv_T(0, n_cold, T_cold); // Derivative w.r.t. T_cold

        real_t v_i = 0;
        if(n_i != 0){
            v_i = sqrt(4*W_i/(3*n_i*ions->GetIonSpeciesMass(iz))); 
        }

        if (derivId == id_W_i){
            if(W_i !=0)
                jac->SetElement(iz, iz, v_i/(2*n_cold*I_i*W_i)); 
        } else if (derivId == id_n_i){
            if(n_i != 0)
                jac->SetElement(iz, iz, -v_i/(2*n_cold*I_i*n_i)); 
        } else if (derivId == id_T_cold){
            jac->SetElement(iz, 0, -v_i/(n_cold*I_i*I_i) * dIdT);
        } else if (derivId == id_n_cold){
            jac->SetElement(iz, 0, -v_i/(n_cold*I_i) * (1/n_cold + dIdn/I_i)); 
        }
        
        return true; 
    }
    
    //This is the linear implicit solver
    void MeanFreePathTerm::SetMatrixElements(FVM::Matrix*, real_t *s){
        s[iz] = lambda_i; 
    }
    
    //This is part of the non-linear solver (F)
    void MeanFreePathTerm::SetVectorElements(real_t *vec, const real_t*){
        SetMatrixElements(nullptr, vec);
    }


