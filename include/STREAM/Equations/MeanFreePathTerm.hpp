#ifndef _STREAM_EQUATION_MEAN_FREE_PATH_TERM_HPP
#define _STREAM_EQUATION_MEAN_FREE_PATH_TERM_HPP


#include "FVM/Equation/EvaluableEquationTerm.hpp" 
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/IonHandler.hpp"

namespace STREAM{ 
	class MeanFreePathTerm : public DREAM::FVM::EvaluableEquationTerm{
	private:
		len_t iz;
		DREAM::FVM::UnknownQuantityHandler *unknowns;
		DREAM::ADAS *adas;
		DREAM::ADASRateInterpolator *interp;
		DREAM::IonHandler *ions;
		real_t lambda_i;
		len_t id_W_i, id_n_i, id_T_cold, id_n_cold;
	public:
	    //Constructor
		MeanFreePathTerm(DREAM::FVM::Grid *g, len_t iz, DREAM::FVM::UnknownQuantityHandler *u, DREAM::ADAS *adas, DREAM::IonHandler *ions);

        virtual void EvaluableTransform(real_t*) override;
		
		//Methods implemented here
		virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 4; }
        
        //Methods implemented in .cpp-file
        virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*) override;
		virtual bool SetJacobianBlock(const len_t, const len_t derivId, DREAM::FVM::Matrix *jac, const real_t*) override;
		virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t *s) override;
		virtual void SetVectorElements(real_t *vec, const real_t*) override;
	};
}

#endif/*_STREAM_EQUATION_MEAN_FREE_PATH_TERM_HPP*/
