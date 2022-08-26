
#include "DREAM/Settings/Settings.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/Equations/DistributionParallelTransport.hpp"


using namespace STREAM;


/**
 * Construct the equation for the hot electron distribution function.
 */
void SimulationGenerator::ConstructEquation_f_re(
	EquationSystem *eqsys, DREAM::Settings *s,
	struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms,
	struct STREAM::OtherQuantityHandler::eqn_terms *stream_terms
) {
	DREAM::SimulationGenerator::ConstructEquation_f_re(eqsys, s, oqty_terms, nullptr);
	
	len_t id_f = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_F_RE);
	
	DREAM::FVM::Operator *eqn = eqsys->GetEquation(id_f)->GetOperatorUnsafe(id_f);
	
        real_t p_cutoff   = s->GetReal("radialgrid/p_cutoff");
	DistributionParallelTransport *DPT = new DistributionParallelTransport(eqsys->GetHotTailGrid(), 
	                                           eqsys->GetUnknownHandler(), eqsys->GetConnectionLength(), 
	                                           eqsys->GetHotTailGrid(), p_cutoff);
   	//eqn->AddTerm(DPT);
   	//stream_terms->DPT=DPT;
}


