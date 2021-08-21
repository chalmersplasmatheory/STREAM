#ifndef _STREAM_OTHER_QUANTITY_HANDLER_HPP
#define _STREAM_OTHER_QUANTITY_HANDLER_HPP

#include <vector>
#include "DREAM/OtherQuantityHandler.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/ChargeExchangeTerm.hpp"
#include "STREAM/Equations/ElectronHeatTransport.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp"
#include "STREAM/Equations/IonRateEquation.hpp"
#include "STREAM/Equations/IonTransport.hpp"
#include "STREAM/Equations/IonHeatTransport.hpp"

namespace STREAM {
    class OtherQuantityHandler : public DREAM::OtherQuantityHandler {
    public:
        struct eqn_terms {
            IonTransport **iontransport=nullptr;
            IonHeatTransport **Wi_iontransport=nullptr;
            ChargeExchangeTerm **Wi_chargeexchange=nullptr;
            ElectronHeatTransport *Tcold_transport=nullptr;
        };
    private:
        ConfinementTime *confinementTime;
        NeutralInflux *neutralInflux;
        PlasmaVolume *plasmaVolume;

        std::vector<IonRateEquation*> ionRateEquations;
        struct eqn_terms *stream_terms;

        len_t id_ni;

    public:
        OtherQuantityHandler(
            ConfinementTime*, NeutralInflux*, PlasmaVolume*,
            std::vector<IonRateEquation*>, struct eqn_terms*,
            DREAM::CollisionQuantityHandler*, DREAM::CollisionQuantityHandler*,
            DREAM::PostProcessor*, DREAM::RunawayFluid*, DREAM::FVM::UnknownQuantityHandler*,
            std::vector<DREAM::UnknownQuantityEquation*>*, DREAM::IonHandler*,
            DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        virtual ~OtherQuantityHandler();

        void DefineQuantitiesSTREAM();
    };
}

#endif/*_STREAM_OTHER_QUANTITY_HANDLER_HPP*/
