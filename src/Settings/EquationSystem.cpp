
#include "STREAM/EquationSystem.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/EqsysInitializer.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"

using namespace STREAM;


#define EQUATIONSYSTEM "eqsys"



/**
 * Construct an EquationSystem object.
 */
EquationSystem *SimulationGenerator::ConstructEquationSystem(
    DREAM::Settings *s, DREAM::FVM::Grid *scalarGrid, DREAM::FVM::Grid *fluidGrid,
    DREAM::ADAS *adas, DREAM::AMJUEL *amjuel, DREAM::NIST *nist,
    EllipticalRadialGridGenerator *ergg
) {
    EquationSystem *eqsys = new EquationSystem(
        scalarGrid, fluidGrid,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, nullptr,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, nullptr, 
        nullptr, nullptr, nullptr
    );
    eqsys->SetEllipticalRadialGridGenerator(ergg);

    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms =
        new DREAM::OtherQuantityHandler::eqn_terms;

    // Timing information
    eqsys->SetTiming(s->GetBool("/output/timingstdout"), s->GetBool("/output/timingfile"));

    // Initialize from previous simulation output?
    const real_t t0 = DREAM::SimulationGenerator::ConstructInitializer(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(
        eqsys, s
    );

    // Construct equations according to settings
    ConstructEquations(eqsys, s, adas, amjuel, nist, oqty_terms);

    // Construct the "other" quantity handler
    ConstructOtherQuantityHandler(eqsys, s, oqty_terms);

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
    EquationSystem *eqsys, DREAM::Settings *s, DREAM::ADAS *adas,
    DREAM::AMJUEL *amjuel, DREAM::NIST *nist,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    DREAM::FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    DREAM::FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();    

    enum DREAM::OptionConstants::momentumgrid_type ht_type = eqsys->GetHotTailGridType();
    enum DREAM::OptionConstants::momentumgrid_type re_type = eqsys->GetRunawayGridType();
    
    // TODO
    ConstructEquation_Ions(eqsys, s, adas, amjuel);
    
    DREAM::IonHandler *ionHandler = eqsys->GetIonHandler();

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
    enum DREAM::FVM::MomentQuantity::pThresholdMode pMode =
        DREAM::FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL;
    enum DREAM::OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum DREAM::OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");

    if (eqsys->HasHotTailGrid() && collfreq_mode == DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) {
        pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        pMode = (enum DREAM::FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
    }

    DREAM::PostProcessor *postProcessor = new DREAM::PostProcessor(
        fluidGrid, unknowns, pThreshold, pMode
    );
    eqsys->SetPostProcessor(postProcessor); 

    // Hot electron quantities
    if (eqsys->HasHotTailGrid()) {
        DREAM::SimulationGenerator::ConstructEquation_f_hot(eqsys, s, oqty_terms);
    }

    // Runaway electron quantities
    if (eqsys->HasRunawayGrid()) {
        DREAM::SimulationGenerator::ConstructEquation_f_re(eqsys, s, oqty_terms, nullptr);
    }

    // Standard equations
    DREAM::SimulationGenerator::ConstructEquation_E_field(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_hot(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_tot(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_ohm(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_re(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_n_cold(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_n_hot(eqsys, s);
    ConstructEquation_T_cold(eqsys, s, adas, amjuel, nist, oqty_terms);   // TODO
    ConstructEquation_lambda_i(eqsys, s, adas); 

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
    ConstructEquation_n_re(eqsys, s, oqty_terms);

    DREAM::SimulationGenerator::ConstructEquation_psi_p(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_psi_edge(eqsys, s);

    ResetPoloidalFluxInitialization(eqsys, s);

    // Helper quantities
    DREAM::SimulationGenerator::ConstructEquation_n_tot(eqsys, s);
    enum DREAM::OptionConstants::eqterm_hottail_mode ht_mode =
        (enum DREAM::OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    enum DREAM::OptionConstants::uqty_f_hot_dist_mode ht_dist_mode =
        (enum DREAM::OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");
    if (ht_mode != DREAM::OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED &&
        ht_dist_mode == DREAM::OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL)
        DREAM::SimulationGenerator::ConstructEquation_tau_coll(eqsys);
        
    eqsys->GetConfinementTime()->Initialize();
}

/**
 * Construct the unknowns of the STREAM equation system.
 */
void SimulationGenerator::ConstructUnknowns(
    EquationSystem *eqsys, DREAM::Settings *s
) {
    DREAM::SimulationGenerator::ConstructUnknowns(
        eqsys, s,
        eqsys->GetScalarGrid(), eqsys->GetFluidGrid(),
        eqsys->GetHotTailGrid(), eqsys->GetRunawayGrid()
    );
    
    len_t nIonSpecies = DREAM::SimulationGenerator::GetNumberOfIonSpecies(s);
    eqsys->SetUnknown(OptionConstants::UQTY_LAMBDA_I, OptionConstants::UQTY_LAMBDA_I_DESC, eqsys->GetFluidGrid(), nIonSpecies);
}

/**
 * In STREAM, since we always have nr=1, we can workaround the inaccuracy
 * in the initialization of psi_p in DREAM by analytically solving for psi_p
 * in the differential form of Ampère's law. In this method we remove the
 * usual initialization rule from DREAM and replace it with a specialized
 * version for STREAM.
 */
void SimulationGenerator::ResetPoloidalFluxInitialization(
    EquationSystem *eqsys, DREAM::Settings *s
) {
    const len_t id_psi_p = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_POL_FLUX);
    const len_t id_psi_edge = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_PSI_EDGE);
    const len_t id_j_tot = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_J_TOT);
    
    DREAM::FVM::RadialGrid *rGrid = eqsys->GetFluidGrid()->GetRadialGrid();
    const real_t a = rGrid->GetMinorRadius();
    const real_t b = s->GetReal("radialgrid/wall_radius");
    const real_t M_inductance = DREAM::PlasmaEdgeToWallInductanceTerm::GetInductance(a, b);

    eqsys->initializer->RemoveRule(id_psi_p);

    std::function<void(DREAM::FVM::UnknownQuantityHandler*, real_t*)> initfunc_PsiP =
        [rGrid,M_inductance,id_psi_edge,id_j_tot](DREAM::FVM::UnknownQuantityHandler *u, real_t *psi_p_init)
    {
        const real_t j_tot = u->GetUnknownData(id_j_tot)[0];
        const real_t psi_edge = u->GetUnknownData(id_psi_edge)[0];
        
        const real_t a = rGrid->GetMinorRadius();
        const real_t r1 = rGrid->GetR(0);
        const real_t dr1 = rGrid->GetDr(0);
        const real_t BdotGradPhiOverB =
            rGrid->GetFSA_1OverR2(0) * rGrid->GetBTorG(0) / rGrid->GetBmin(0);
        const real_t Vp_1 = rGrid->GetVpVol(0);
        const real_t Vp_32 = rGrid->GetVpVol_f(1);

        psi_p_init[0] = psi_edge -
            (2*M_PI*DREAM::Constants::mu0 * BdotGradPhiOverB * (a-r1)*dr1*Vp_1) /
            (Vp_32 * rGrid->GetFSA_NablaR2OverR2_f(1))
            * j_tot;
    };

    eqsys->initializer->AddRule(
        id_psi_p,
        DREAM::EqsysInitializer::INITRULE_EVAL_FUNCTION,
        initfunc_PsiP,
        // Dependencies
        id_j_tot,
        id_psi_edge
    );
}

