#ifndef _PLASMA_VOLUME_HPP
#define _PLASMA_VOLUME_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"


namespace STREAM{
    class PlasmaVolume{ //: public { should it be a subclass?
    private:
        DREAM::FVM::Grid *grid; //Is this correct when it is not a subclass?
        len_t iz;
        real_t vessel_vol;
        DREAM::FVM::UnknownQuantityHandler *unknowns;
        EllipticalRadialGridGenerator *radials;
        len_t id_lambda_i;
    public:
        PlasmaVolume(DREAM::FVM::Grid *g, len_t iz, real_t vessel_vol, DREAM::FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r); 
        real_t GetPlasmaVolume() const; 
        real_t GetNeutralVolume(const len_t iz);
    };
}

#endif/*_PLASMA_VOLUME_HPP*/
