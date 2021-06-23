#include "STREAM/EquationSystem.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"

using namespace STREAM;


#define EQUATIONSYSTEM "eqsys"

void SimulationGenerator::DefineOptions_wall(DREAM::Settings *s) {
    s->DefineSetting(
        "radialgrid/wall/c1",
        "Coefficients for deuterium recycling",
        (real_t)1.1
    );
    s->DefineSetting(
        "radialgrid/wall/c2",
        "Coefficients for deuterium recycling",
        (real_t)0.09
    );
    s->DefineSetting(
        "radialgrid/wall/c3",
        "Coefficients for deuterium recycling",
        (real_t)0.1
    );
    
    s->DefineSetting(
        "radialgrid/wall/vessel_volume", 
        "The vacuum vessel volume",
        (real_t)0
    );
}


/**
 * Construct an EquationSystem object.
 */
EquationSystem *SimulationGenerator::ConstructEquationSystem(
    DREAM::Settings *s, DREAM::FVM::Grid *scalarGrid, DREAM::FVM::Grid *fluidGrid,
    DREAM::ADAS *adas, DREAM::AMJUEL *amjuel, DREAM::NIST *nist
) {
    EquationSystem *eqsys = new EquationSystem(
        scalarGrid, fluidGrid,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, nullptr,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, nullptr, 
        nullptr, nullptr
    );

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
    
    DREAM::IonHandler *ionHandler = eqsys->GetIonHandler();
    
    // Confinement time 
    EllipticalRadialGridGenerator *r = eqsys->GetEllipticalRadialGridGenerator(); 
    
    real_t l_MK2 = s->GetReal("radialgrid/wall_radius");
    ConfinementTime *confinementTime = new ConfinementTime(
        unknowns, r, l_MK2
    );
    eqsys->SetConfinementTime(confinementTime);
    
    //Plasma Volume
    real_t vessel_volume = s->GetReal("radialgrid/wall/vessel_volume");
    if (vessel_volume == 0){
        throw DREAM::SettingsException(
            "Vessel volume is unspecified" //Is this an ok exception? 
        );
    }
    PlasmaVolume *volumes = new PlasmaVolume(
        fluidGrid, vessel_volume, unknowns, r, adas, ionHandler
    );
    eqsys->SetPlasmaVolume(volumes);
    
    /*
    // Neutral influx 
    SputteredRecycledCoefficient *SRC = eqsys->GetSputteredRecycledCoefficient(); //Korrekt?

    real_t c1 = s->GetReal("radialgrid/wall/c1"); 
    real_t c2 = s->GetReal("radialgrid/wall/c2"); 
    real_t c3 = s->GetReal("radialgrid/wall/c3"); 
    NeutralInflux *neutralInflux = new NeutralInflux(
        ionHandler, SRC, confinementTime, volumes, c1, c2, c3 
    );
    eqsys->SetNeutralInflux(neutralInflux); 
    */
    
    // TODO
    ConstructEquation_Ions(eqsys, s, adas, amjuel);

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
    DREAM::SimulationGenerator::ConstructEquation_n_re(eqsys, s, oqty_terms);

    DREAM::SimulationGenerator::ConstructEquation_psi_p(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_psi_edge(eqsys, s);

    // Helper quantities
    DREAM::SimulationGenerator::ConstructEquation_n_tot(eqsys, s);
    enum DREAM::OptionConstants::eqterm_hottail_mode ht_mode =
        (enum DREAM::OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    enum DREAM::OptionConstants::uqty_f_hot_dist_mode ht_dist_mode =
        (enum DREAM::OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");
    if (ht_mode != DREAM::OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED &&
        ht_dist_mode == DREAM::OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL)
        DREAM::SimulationGenerator::ConstructEquation_tau_coll(eqsys);
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
    eqsys->SetUnknown(OptionConstants::UQTY_LAMBDA_I, OptionConstants::UQTY_LAMBDA_I_DESC, eqsys->GetFluidGrid());
}

