#include "STREAM/EquationSystem.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "STREAM/Equations/ElectronTransportDiffusion.hpp"

using namespace STREAM;

#define MODULENAME "eqsys/n_re"

void SimulationGenerator::DefineOptions_n_re(
    DREAM::Settings *s
) {
    DREAM::SimulationGenerator::DefineOptions_n_re(s);
}

void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, DREAM::Settings *s,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    DREAM::SimulationGenerator::ConstructEquation_n_re(eqsys, s, oqty_terms, nullptr);
    
    len_t id_n_re = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_N_RE);
    
    DREAM::FVM::Operator *Op_nRE = eqsys->GetEquation(id_n_re)->GetOperatorUnsafe(id_n_re);
    Op_nRE->AddTerm(new ElectronTransportDiffusion(eqsys->GetFluidGrid(), eqsys->GetEllipticalRadialGridGenerator(), eqsys->GetConfinementTime(), eqsys->GetUnknownHandler()));
}
