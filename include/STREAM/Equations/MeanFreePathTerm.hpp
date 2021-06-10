#ifndef _STREAM_EQUATION_MEAN_FREE_PATH_TERM_HPP
#define _STREAM_EQUATION_MEAN_FREE_PATH_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp" 
#include "FVM/UnknownQuantityHandler.hpp"

namespace STREAM{
	class MeanFreePathTerm : public DREAM::FVM::EvaluableEquationTerm{
	private:
		len_t iz;
		DREAM::FVM::UnknownQuantityHandler *unknowns;
		DREAM::ADAS *adas;
		DREAM::ADASRateInterpolator *interp;
	public:
		MeanFreePathTerm(DREAM::FVM::Grid *g, len_t iz, DREAM::FVM::UnknownQuantityHandler *u, DREAM::ADAS *adas); //Constructor declaration, TODO: unknown id here in cpp 
		virtual void Rebuild(const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*) override;
		virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        	virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 3; }
		virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*) override;
		virtual void SetMatrixElements(FVM::Matrix *mat, real_t*) override;
		virtual void SetVectorElements(real_t *vec, const real_t *xi) override;
	};
}

#endif/*_STREAM_EQUATION_MEAN_FREE_PATH_TERM_HPP*/
