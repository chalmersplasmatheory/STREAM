/**
 * Common routines for adding transport terms to equations.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


using namespace STREAM;
using namespace std;


/**
 * Define options for the transport model.
 */
void SimulationGenerator::DefineOptions_Transport(
    const string& mod, DREAM::Settings *s, bool kinetic, const string& subname
) {
    DREAM::SimulationGenerator::DefineOptions_Transport(mod, s, kinetic, subname);

    s->DefineSetting(mod + "/" + subname + "/type", "Type of transport model to use.", (int_t)OptionConstants::EQTERM_TRANSPORT_NONE);
}

/**
 * Construct transport terms.
 *
 * oprtr:        Operator to add the transport term to.
 * mod:          Name of module to load settings from.
 * grid:         Grid on which the operator will be defined.
 * momtype:      Type of momentum grid.
 * unknowns:     Unknown quantity handler.
 * s:            Object to load settings from.
 * kinetic:      If 'true', the term is assumed to be applied to a kinetic
 *               grid and the transport coefficient is expected to be 4D
 *               (time + radius + p1 + p2).
 * heat:         Indicates that the quantity to which this operator is
 *               applied represents heat (i.e. a temperature) and that
 *               operators for heat transport should be used where available.
 * advective_bc: If not 'nullptr', sets this pointer to the newly created
 *               advective B.C. (if any, otherwise nullptr). This can be used
 *               to gain access to the B.C. in, for example, the
 *               OtherQuantityHandler.
 * diffusive_bc: If not 'nullptr', sets this pointer to the newly created
 *               diffusive B.C. (if any, otherwise nullptr). This can be used
 *               to gain access to the B.C. in, for example, the
 *               OtherQuantityHandler.
 * oqty_terms:   List of other quantity terms to save (only for SvenssonTransport).
 * subname:      Name of section in the settings module which the transport
 *               settings are stored.
 * 
 * returns: true if non-zero transport, otherwise false
 */
bool SimulationGenerator::ConstructTransportTerm(
    DREAM::FVM::Operator *oprtr, const string& mod, DREAM::FVM::Grid *grid,
    enum DREAM::OptionConstants::momentumgrid_type momtype,
    DREAM::EquationSystem *eqsys,
    DREAM::Settings *s, bool kinetic, bool heat,
    DREAM::TransportAdvectiveBC **advective_bc,
    DREAM::TransportDiffusiveBC **diffusive_bc,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms,
    const string& subname
) {
    enum OptionConstants::eqterm_transport_type type =
        (enum OptionConstants::eqterm_transport_type)s->GetInteger(mod + "/type");

    if (type == OptionConstants::EQTERM_TRANSPORT_DYON) {
        // TODO
        
        //return true;
    } else if (type != OptionConstants::EQTERM_TRANSPORT_NONE)
        return DREAM::SimulationGenerator::ConstructTransportTerm(
            oprtr, mod, grid, momtype, eqsys, s, kinetic,
            heat, advective_bc, diffusive_bc, oqty_terms,
            subname
        );

    return false;
}

