/**
 * Implementation of an equation term for the mean free path of
 * neutrals of species with index iz
 */

#include "STREAM/Equations/MeanFreePathTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

using namespace STREAM;
using namespace DREAM;

    MeanFreePathTerm::MeanFreePathTerm(FVM::Grid *g, len_t iz, FVM::UnknownQuantityHandler *u, ADAS *adas): FVM::EvaluableEquationTerm(g), iz(iz), unknowns(u), adas(adas) interp(adas->getSCD(iz)){
    id_Wi = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    id_Ni
    id_Te
    id_ne
    

    }
    
    bool MeanFreePathTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(id_Wi==derivId)
    jac->set(iz,iz,) //dlambdai/dWi
    
    ne = unknowns->GetUnknownData(id_ne)[0];
    Te = unknowns->GetUnknownData(id_Te)[0];
    Eval_deriv_n(0, ne, Te); // Derivative evaluation
    Eval_deriv_T(0, ne, Te); // Derivative evaulation
    Eval((0, ne, Te); //Evaluate I_i^(0)
    }
    
    
    
    


