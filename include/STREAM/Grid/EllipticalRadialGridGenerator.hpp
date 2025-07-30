#ifndef _STREAM_GRID_ELLIPTICAL_RADIAL_GRID_GENERATOR_HPP
#define _STREAM_GRID_ELLIPTICAL_RADIAL_GRID_GENERATOR_HPP

#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "FVM/Interpolator1D.hpp"

namespace STREAM { 
    class EllipticalRadialGridGenerator : public DREAM::FVM::RadialGridGenerator {
    protected:
        DREAM::FVM::Interpolator1D *a, *B0, *kappa, *delta;

        // Radial grid points considered
        real_t r, r_f[2], dr, dr_f=0, R0;

        // Most recently evaluated minor radius, magnetic field and elongation...
        real_t currA=0, currB0=0, currKappa=0, currTriang=0;

		len_t ntheta_out = 80;

    public:
        EllipticalRadialGridGenerator(
            DREAM::FVM::Interpolator1D*, DREAM::FVM::Interpolator1D*,
            DREAM::FVM::Interpolator1D*, DREAM::FVM::Interpolator1D*, 
            real_t
        );
        ~EllipticalRadialGridGenerator();

        virtual bool NeedsRebuild(const real_t) const override;
        virtual bool Rebuild(const real_t, DREAM::FVM::RadialGrid*) override;

        virtual real_t JacobianAtTheta(const len_t, const real_t) override;
        virtual real_t ROverR0AtTheta(const len_t, const real_t) override;
        virtual real_t NablaR2AtTheta(const len_t, const real_t) override;

        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t) override;
        virtual real_t ROverR0AtTheta_f(const len_t, const real_t) override;
        virtual real_t NablaR2AtTheta_f(const len_t, const real_t) override;

        virtual void EvaluateGeometricQuantities(
            const len_t, const real_t,
            real_t&, real_t&, real_t&, real_t&
        ) override;
        virtual void EvaluateGeometricQuantities_fr(
            const len_t, const real_t,
            real_t&, real_t&, real_t&, real_t&
        ) override;

        void GetRThetaFromCartesian(
            real_t*, real_t*, real_t, real_t, real_t, real_t, real_t
        );
        void GetGradRCartesian(
            real_t*, real_t, real_t
        );
        real_t FindClosestApproach(
            real_t, real_t, real_t, real_t, real_t, real_t
        );
        
        real_t GetMinorRadius(){ return currA; }
        real_t GetElongation(){ return currKappa; }
        real_t GetTriangularity(){return currTriang; }
        real_t GetMagneticField(){ return currB0; }
        real_t GetMajorRadius(){ return R0; }
        
        void GetRThetaPhiFromCartesian(real_t*, real_t*, real_t*, real_t , real_t , real_t, real_t, real_t ) override {};
        void GetGradRCartesian(real_t* ,real_t, real_t, real_t) override {};

		virtual const real_t GetZ0() override { return 0; }
		virtual const len_t GetNPsi() override { return this->GetNr(); }
		virtual const len_t GetNTheta() override { return this->ntheta_out; }
		virtual const real_t *GetFluxSurfaceRMinusR0() override;
		virtual const real_t *GetFluxSurfaceRMinusR0_f() override;
		virtual const real_t *GetFluxSurfaceZMinusZ0() override;
		virtual const real_t *GetFluxSurfaceZMinusZ0_f() override;
		virtual const real_t *GetPoloidalAngle() override;
    };
}

#endif/*_STREAM_GRID_ELLIPTICAL_RADIAL_GRID_GENERATOR_HPP*/
