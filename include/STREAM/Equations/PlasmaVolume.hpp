#ifndef _PLASMA_VOLUME_HPP
#define _PLASMA_VOLUME_HPP

#include "FVM/Grid/Grid.hpp"

namespace STREAM{
    class PlasmaVolume{ //: public { should it be a subclass?
    private:
        len_t iz;
        real_t vessel_vol;
        DREAM::FVM::UnknownQuantityHandler *unknowns;
        DREAM::FVM::Grid *grid; //Is this correct when it is not a subclass?
    public:
        PlasmaVolume(DREAM::FVM::Grid *g, len_t iz, real_t vessel_vol, DREAM::FVM::UnknownQuantityHandler *u); 
        real_t GetPlasmaVolume() const; 
        real_t GetNeutralVolume(const len_t iz);
    }
}

#endif/*_PLASMA_VOLUME_HPP*/
