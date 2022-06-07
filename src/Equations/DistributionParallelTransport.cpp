
#include "DREAM/Constants.hpp"
#include "STREAM/Equations/DistributionParallelTransport.hpp"

using namespace STREAM;
using namespace DREAM;

/**
 * Constructor
 */
DistributionParallelTransport::DistributionParallelTransport(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u,
    ConnectionLength *CL, FVM::Grid *operandGrid, real_t pcutoff
    ) : DiagonalComplexTerm(g, u, operandGrid), CL(CL), pcutoff(pcutoff) {
    
    this->DiagonalTerm::SetName("DistributionParallelTransport");
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DistributionParallelTransport::SetMatrixElements(FVM::Matrix *mat, real_t*) { 
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    for(len_t j = 0; j < n2; j++) {
        for(len_t i = 0; i < n1; i++) {
            mat->SetElement(i, j, weights[j*n1+i]); 
        }
    }
}

/**
 * Set function vector for this term.
 */
void DistributionParallelTransport::SetVectorElements(real_t *vec, const real_t *x) { 
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    for(len_t j = 0; j < n2; j++) {
        for(len_t i = 0; i < n1; i++) {
            vec[j*n1+i] += weights[j*n1+i] * x[j*n1+i]; 
        }
    }
}

/**
 * Set of weights for this diagonal term
 */
void DistributionParallelTransport::SetWeights(){ 
    if (id_Ip == 0) this->Initialize();
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    real_t xi, p; 
    real_t frac = 0;
    for(len_t j = 0; j < n2; j++) {
        for(len_t i = 0; i < n1; i++) {
            p  = grid->GetMomentumGrid(0)->GetP(i, j);
            if (p > pcutoff) {
	        xi = grid->GetMomentumGrid(0)->GetXi0(i, j); 
	    
	        real_t vpar = std::abs(p * xi / sqrt(1 + p*p)); 
	        real_t Lfinv = CL->EvaluateInverseConnectionLength(0); // ?? Should do arbitrary ir?
	    
	        weights[j*n1+i] =- vpar * Lfinv; 
	        
	        frac += 1.0 / (n2 * n1);
    
	    } else {
	        weights[j*n1+i] = 0;
	    }
        }
    }
    //printf("p_cutoff=%.2e, frac=%.2e\n", pcutoff, frac);
}

/**
 * Set of derivatives of weights for this diagonal term
 */
void DistributionParallelTransport::SetDiffWeights(len_t derivId, len_t){
    ResetDiffWeights();
    
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    real_t xi, p;
    for(len_t j = 0; j < n2; j++) {
        for(len_t i = 0; i < n1; i++) {
            p  = grid->GetMomentumGrid(0)->GetP(i, j);
            if (p > pcutoff) {
                xi = grid->GetMomentumGrid(0)->GetXi0(i, j);
            
                real_t vpar = std::abs(p * xi / sqrt(1 + p*p)); 
                real_t dLfinv; 
            
                if(derivId == id_Ip) {
                    dLfinv = CL->EvaluateInverseConnectionLength_dIp(0); // ?? Should do arbitrary ir?
                } else if(derivId == id_Iwall) {
	            dLfinv = CL->EvaluateInverseConnectionLength_dIwall(0); // ?? Should do arbitrary ir?
                } 
                diffWeights[j*n1+i] =- vpar * dLfinv; 
            } else {
                diffWeights[j*n1+i] = 0;
            }
        }
    }
}


/**
 * Set all diffweights to 0.
 */
void DistributionParallelTransport::ResetDiffWeights(){
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    for(len_t j = 0; j < n2; j++) {
        for(len_t i = 0; i < n1; i++) {
            diffWeights[j*n1+i] = 0;
        }
    }
}

/**
 * Allocate differentiation coefficients.
 */
void DistributionParallelTransport::AllocateDiffWeights() {
    DeallocateDiffWeights();
    
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    
    diffWeights = new real_t[n2*n1];
    ResetDiffWeights();
}

/**
 * Deallocate differentiation coefficients.
 */
void DistributionParallelTransport::DeallocateDiffWeights() {
    if(diffWeights != nullptr)
        delete [] diffWeights;
}

void DistributionParallelTransport::Initialize() {
    id_Ip    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    id_Iwall = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_WALL);
}
