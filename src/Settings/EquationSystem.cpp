
#include "DREAM/EquationSystem.hpp"
#include "DREAM/OtherQuantityHandler.hpp"


using namespace STREAM;


#define EQUATIONSYSTEM "eqsys"


/**
 * Construct an EquationSystem object.
 */
DREAM::EquationSystem *SimulationGenerator::ConstructEquationSystem(
    DREAM::Settings *s, DREAM::FVM::Grid *scalarGrid, DREAM::FVM::Grid *fluidGrid,
    enum DREAM::OptionConstants::momentumgrid_type ht_type, DREAM::FVM::Grid *hottailGrid,
    enum DREAM::OptionConstants::momentumgrid_type re_type, DREAM::FVM::Grid *runawayGrid,
    DREAM::ADAS *adas, DREAM::AMJUEL *amjuel, DREAM::NIST *nist
) {
    DREAM::EquationSystem *eqsys = new DREAM::EquationSystem(
        scalarGrid, fluidGrid,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, nullptr,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, nullptr,
    );

    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms =
        new DREAM::OtherQuantityHandler::eqn_terms;

    // Timing information
    eqsys->SetTiming(s->GetBool("/output/timingstdout"), s->GetBool("/output/timingfile"));

    // Initialize from previous simulation output?
    const real_t t0 = ConstructInitializer(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(
        eqsys, s, scalarGrid, fluidGrid, hottailGrid, runawayGrid,
        adas, amjuel, nist, oqty_terms
    );

    // Construct equations according to settings
    ConstructEquations(eqsys, s, adas, amjuel, nist, oqty_terms);

    // Figure out which unknowns must be part of the matrix,
    // and set initial values for those quantities which don't
    // yet have an initial value.
    eqsys->ProcessSystem(t0);

    // (these must be initialized AFTER calling 'ProcessSystem()' on
    // the equation system, since we need to know which unknowns are
    // "non-trivial", i.e. need to show up in the solver matrices,
    // in order to build them)

    // Construct the time stepper
    DREAM::SimulationGenerator::ConstructTimeStepper(eqsys, s);

    // Construct solver (must be done after processing equation system,
    // since we need to know which unknowns are "non-trivial",
    // i.e. need to show up in the solver matrices)
    DREAM::SimulationGenerator::ConstructSolver(eqsys, s);

    return eqsys;
}

/**
 * Construct equations to solve.
 */
void SimulationGenerator::ConstructEquations(
    DREAM::EquationSystem *eqsys, DREAM::Settings *s, DREAM::ADAS *adas,
    DREAM::AMJUEL *amjuel, DREAM::NIST *nist,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_tems
) {
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    DREAM::FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    DREAM::FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    enum DREAM::OptionConstants::momentumgrid_type ht_type = eqsys->GetHotTailGridType();
    enum DREAM::OptionConstants::momentumgrid_type re_type = eqsys->GetRunawayGridType();

    // TODO
    ConstructEquation_Ions(eqsys, s, adas, amjuel);

    IonHandler *ionHandler = eqsys->GetIonHandler();

    // Construct collision quantity handlers
    if (hottailGrid != nullptr)
        eqsys->SetHotTailCollisionHandler(
            DREAM::SimulationGenerator::ConstructCollisionQuantityHandler(
                ht_type, hottailGrid, unknowns, ionHandler, s
            )
        );
    if (runawayGrid != nullptr)
        eqsys->SetRunawayCollisionHandler(
            DREAM::SimulationGenerator::ConstructCollisionQuantityHandler(
                re_type, runawayGrid, unknowns, ionHandler, s
            )
        );

    DREAM::SimulationGenerator::ConstructRunawayFluid(
        fluidGrid, unknowns, ionHandler,
        re_type, eqsys, s
    );

    // Build post-processor
    real_t pThreshold = 0.0;
    FVM::MomentQuantity::pThresholdMode pMode =
        FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL;
    enum DREAM::OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum DREAM::OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");

    if (eqsys->HasHotTailGrid() && collfreq_mode == DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) {
        pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        pMode = (DREAM::FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
    }

    PostProcessor *postProcessor = new PostProcessor(
        fluidGrid, unknowns, pThreshold, pMode
    );
    eqsys->SetPostProcessor(postProcessor);

    // Hot electron quantities
    if (eqsys->HasHotTailGrid()) {
        DREAM::SimulationGenerator::ConstructEquation_f_hot(eqsys, s, oqty_terms);
    }

    // Runaway electron quantities
    if (eqsys->HasRunawayGrid()) {
        DREAM::SimulationGenerator::ConstructEquation_f_re(eqsys, s, oqty_terms);
    }

    // Standard equations
    DREAM::SimulationGenerator::ConstructEquation_E_field(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_hot(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_tot(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_ohm(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_re(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_n_cold(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_n_hot(eqsys, s);
    ConstructEquation_T_cold(eqsys, s);   // TODO

    enum DREAM::OptionConstants::uqty_T_i_eqn typeTi =
        (DREAM::OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if (typeTi != DREAM::OptionConstants::UQTY_T_I_INCLUDE)
        throw DREAM::SettingsException(
            "T_i not included: STREAM requires the ion temperatures to be evolved."
        );

    DREAM::SimulationGenerator::ConstructEquation_Ion_Ni(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_T_i(eqsys, s);

    // NOTE: The runaway number density may depend explicitly on
    // the hot-tail equation and must therefore be constructed
    // AFTER the call to 'ConstructEquation_f_hot()'.
    DREAM::SimulationGenerator::ConstructEquation_n_re(eqsys, s, oqty_terms);

    DREAM::SimulationGenerator::ConstructEquation_psi_p(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_psi_edge(eqsys, s);

    // Helper quantities
    DREAM::SimulationGenerator::ConstructEquation_n_tot(eqsys, s);
    enum DREAM::OptionConstants::eqterm_hottail_mode ht_mode = (enum OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    enum DREAM::OptionConstants::uqty_f_hot_dist_mode ht_dist_mode = (enum OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");
    if (ht_mode != DREAM::OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED &&
        ht_dist_mode == DREAM::OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL)
        DREAM::SimulationGenerator::ConstructEquation_tau_coll(eqsys);
}

/**
 * Construct the unknowns of the STREAM equation system.
 */
void SimulationGenerator::ConstructUnknowns(
    DREAM::EquationSystem *eqsys, DREAM::Settings *s, DREAM::ADAS *adas,
    DREAM::AMJUEL *amjuel, DREAM::NIST *nist,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_tems
) {
    DREAM::SimulationGenerator::ConstructUnknowns(
        eqsys, s,
        eqsys->GetScalarGrid(), eqsys->GetFluidGrid(),
        eqsys->GetHotTailGrid(), eqsys->GetRunawayGrid(),
        adas, nist, amjuel, oqty_terms
    );
}

