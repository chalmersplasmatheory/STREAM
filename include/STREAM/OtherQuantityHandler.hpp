#ifndef _STREAM_OTHER_QUANTITY_HANDLER_HPP
#define _STREAM_OTHER_QUANTITY_HANDLER_HPP

#include "DREAM/OtherQuantityHandler.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"

namespace STREAM {
    class OtherQuantityHandler : public DREAM::OtherQuantityHandler {
    private:
        ConfinementTime *confinementTime;
        NeutralInflux *neutralInflux;
        PlasmaVolume *plasmaVolume;

    public:
        OtherQuantityHandler(
            ConfinementTime*, NeutralInflux*, PlasmaVolume*,
            DREAM::CollisionQuantityHandler*, DREAM::CollisionQuantityHandler*,
            DREAM::PostProcessor*, DREAM::RunawayFluid*, DREAM::FVM::UnknownQuantityHandler*,
            std::vector<DREAM::UnknownQuantityEquation*>*, DREAM::IonHandler*,
            DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            struct eqn_terms*
        );
        virtual ~OtherQuantityHandler();

        void DefineQuantitiesSTREAM();
    };
}

#endif/*_STREAM_OTHER_QUANTITY_HANDLER_HPP*/
