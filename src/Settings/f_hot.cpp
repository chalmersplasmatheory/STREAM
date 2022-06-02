
#include "DREAM/Settings/Settings.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


using namespace STREAM;


/**
 * Construct the equation for the hot electron distribution function.
 */
void SimulationGenerator::ConstructEquation_f_hot(
	EquationSystem *eqsys, DREAM::Settings *s,
	struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
	ConstructEquation_f_hot(eqsys, s, oqty_terms);
}


