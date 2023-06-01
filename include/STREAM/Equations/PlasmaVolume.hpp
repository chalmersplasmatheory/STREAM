#ifndef _PLASMA_VOLUME_HPP
#define _PLASMA_VOLUME_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/IonHandler.hpp"


namespace STREAM{
    class PlasmaVolume{
    private:
        DREAM::FVM::Grid *grid; 
        real_t vessel_vol;
        DREAM::FVM::UnknownQuantityHandler *unknowns;
        EllipticalRadialGridGenerator *radials;
        DREAM::ADAS *adas;
		DREAM::IonHandler *ions;
        len_t id_lambda_i, id_W_i, id_n_i, id_T_cold, id_n_cold;
    public:
        PlasmaVolume(DREAM::FVM::Grid *g, real_t vessel_vol, DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r, DREAM::ADAS *adas, DREAM::IonHandler *ions); 
		real_t GetPlasmaCrossSection() const;
        real_t GetPlasmaVolume() const; 
        real_t GetNeutralVolume(const len_t iz);
        real_t GetTotalNeutralVolume(const len_t iz);
        //real_t GetNeutralVolume_dT(const len_t iz);
        //real_t GetNeutralVolume_dn(const len_t iz);
        real_t GetNeutralVolume_dLambdai(const len_t iz);
        //real_t GetTotalNeutralVolume_dT(const len_t iz);
        //real_t GetTotalNeutralVolume_dn(const len_t iz);
        real_t GetTotalNeutralVolume_dLambdai(const len_t iz);

        real_t GetVesselVolume() { return vessel_vol; }
    };
}

#endif/*_PLASMA_VOLUME_HPP*/
