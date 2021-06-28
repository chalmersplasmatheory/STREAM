#include "STREAM/EquationSystem.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace STREAM;

#define MODULENAME "eqsys/n_re"

void SimulationGenerator::DefineOptions_n_re(
    Settings *s
) {
    DREAM::SimulationGenerator::DefineOptions_n_re(s);
}

void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, Settings *s,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    DREAM::SimulationGenerator::ConstructEquation_n_re(eqsys, s, oqty_terms);
    
    len_t id_n_re = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    
    const DREAM::FVM::Operator *Op_nRE = eqsys->GetEquation(id_n_re)->GetOperator(id_n_re);
    
    
}
