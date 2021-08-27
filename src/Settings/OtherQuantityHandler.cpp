/**
 * Construct an 'OtherQuantityHandler' object.
 */

#include <iostream>
#include <string>
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/PostProcessor.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "STREAM/OtherQuantityHandler.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


#define MODULENAME "other"

using namespace STREAM;
using namespace std;

/**
 * Construct an 'OtherQuantityHandler' object.
 */
void SimulationGenerator::ConstructOtherQuantityHandler(
    EquationSystem *eqsys, DREAM::Settings *s,
    struct OtherQuantityHandler::eqn_terms *stream_terms,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    OtherQuantityHandler *oqh = new OtherQuantityHandler(
        eqsys->GetConfinementTime(), eqsys->GetNeutralInflux(),
        eqsys->GetPlasmaVolume(), eqsys->GetRunawayElectronConfinementTime(),
        eqsys->GetIonRateEquations(), stream_terms,
        eqsys->GetHotTailCollisionHandler(), eqsys->GetRunawayCollisionHandler(),
        eqsys->GetPostProcessor(), eqsys->GetREFluid(), eqsys->GetUnknownHandler(),
        eqsys->GetEquations(), eqsys->GetIonHandler(), eqsys->GetFluidGrid(),
        eqsys->GetHotTailGrid(), eqsys->GetRunawayGrid(), eqsys->GetScalarGrid(),
        oqty_terms
    );

    const vector<string> other = s->GetStringList(MODULENAME "/include");

    if (other.size() == 1 && other[0] == "all")
        oqh->RegisterAllQuantities();
    else {
        for (auto it = other.begin(); it != other.end(); it++)
            oqh->RegisterQuantity(*it);
    }

    eqsys->SetOtherQuantityHandler(oqh);
}

