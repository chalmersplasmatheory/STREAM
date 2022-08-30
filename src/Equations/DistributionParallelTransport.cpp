
#include "DREAM/Constants.hpp"
#include "STREAM/Equations/DistributionParallelTransport.hpp"

using namespace STREAM;
using namespace DREAM;

/**
 * Constructor
 */
DistributionParallelTransport::DistributionParallelTransport(
    FVM::Grid *g, FVM::UnknownQuantityHandler *u,
    ConnectionLength *CL, RunawayElectronConfinementTime *REC, FVM::Grid *operandGrid, real_t pcutoff
    ) : DiagonalComplexTerm(g, u, operandGrid), CL(CL), REC(REC), pcutoff(pcutoff) {
    
    this->DiagonalTerm::SetName("DistributionParallelTransport");
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DistributionParallelTransport::SetMatrixElements(FVM::Matrix *mat, real_t*) { 
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    const len_t Nr = this->grid->GetNr();
    for (len_t ir = 0; ir < Nr; ir++) {
        for(len_t j = 0; j < n2; j++) {
            for(len_t i = 0; i < n1; i++) {
                mat->SetElement(i, j, weights[ir*n2*n1+j*n1+i]); 
            }
        }
    }
}

/**
 * Set function vector for this term.
 */
void DistributionParallelTransport::SetVectorElements(real_t *vec, const real_t *x) { 
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    const len_t Nr = this->grid->GetNr();
    for (len_t ir = 0; ir < Nr; ir++) {
        for(len_t j = 0; j < n2; j++) {
            for(len_t i = 0; i < n1; i++) {
                vec[ir*n2*n1+j*n1+i] += weights[ir*n2*n1+j*n1+i] * x[ir*n2*n1+j*n1+i]; 
            }
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
    const len_t Nr = this->grid->GetNr();
    real_t xi, p; 
    real_t c = 299792458.0;
    
    for (len_t ir = 0; ir < Nr; ir++) {
        for(len_t j = 0; j < n2; j++) {
            for(len_t i = 0; i < n1; i++) {
                p  = grid->GetMomentumGrid(0)->GetP(i, j);
                if (p > pcutoff) {
	            xi = grid->GetMomentumGrid(0)->GetXi0(i, j); 
	    
	            real_t vpar = std::abs(c * p * xi / sqrt(1 + p*p)); 
                    real_t Lfinv = CL->EvaluateInverseConnectionLength(0); // ?? Should do arbitrary ir?
                    real_t tauinv = REC->EvaluateInverse(0); // ?? Should do arbitrary ir?
	            weights[ir*n2*n1+j*n1+i] =- vpar * tauinv / c; 
                } else {
	            weights[ir*n2*n1+j*n1+i] = 0;
	        }
	    }
        }
    }
}

/**
 * Set of derivatives of weights for this diagonal term
 */
void DistributionParallelTransport::SetDiffWeights(len_t derivId, len_t){
    ResetDiffWeights();
    
    len_t n2 = grid->GetMomentumGrid(0)->GetNp2();
    len_t n1 = grid->GetMomentumGrid(0)->GetNp1();
    const len_t Nr = this->grid->GetNr();
    real_t xi, p; 
    real_t c = 299792458.0;
    
    for (len_t ir = 0; ir < Nr; ir++) {
        for(len_t j = 0; j < n2; j++) {
            for(len_t i = 0; i < n1; i++) {
                p  = grid->GetMomentumGrid(0)->GetP(i, j);
                if (p > pcutoff) {
                    xi = grid->GetMomentumGrid(0)->GetXi0(i, j);
            
                    real_t vpar = std::abs(c * p * xi / sqrt(1 + p*p)); 
                    real_t dLfinv, dtauinv; 
            
                    if(derivId == id_Ip) {
                        dLfinv = CL->EvaluateInverseConnectionLength_dIp(0); // ?? Should do arbitrary ir?
                        dtauinv = REC->Evaluate_dIp(0); // ?? Should do arbitrary ir?
                    } else if(derivId == id_Iwall) {
	                dLfinv = CL->EvaluateInverseConnectionLength_dIwall(0); // ?? Should do arbitrary ir?
	                dtauinv = REC->Evaluate_dIp(0); // ?? Should do arbitrary ir?
                    } else if(derivId = id_Efield) {
                        dtauinv = REC->Evaluate_dE(0);
                    }
                    diffWeights[ir*n2*n1+j*n1+i] =- vpar * dtauinv / c; 
                } else {
                   diffWeights[ir*n2*n1+j*n1+i] = 0;
               }
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
    const len_t Nr = this->grid->GetNr();
    for (len_t ir = 0; ir < Nr; ir++) {
        for(len_t j = 0; j < n2; j++) {
            for(len_t i = 0; i < n1; i++) {
                diffWeights[ir*n2*n1+j*n1+i] = 0;
            }
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
    id_Efield = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
}
