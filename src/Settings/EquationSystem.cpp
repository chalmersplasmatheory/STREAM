
#include "STREAM/EquationSystem.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/EqsysInitializer.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Equations/OpticalThickness.hpp"

using namespace STREAM;


#define EQUATIONSYSTEM "eqsys"



/**
 * Construct an EquationSystem object.
 */
EquationSystem *SimulationGenerator::ConstructEquationSystem(
    DREAM::Settings *s, DREAM::FVM::Grid *scalarGrid, DREAM::FVM::Grid *fluidGrid,
	DREAM::FVM::Grid *hottailGrid, DREAM::FVM::Grid *runawayGrid, DREAM::ADAS *adas, DREAM::AMJUEL *amjuel,
	DREAM::NIST *nist, EllipticalRadialGridGenerator *ergg
) {
    EquationSystem *eqsys = new EquationSystem(
        scalarGrid, fluidGrid,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, hottailGrid,
        DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI, runawayGrid, 
        s, ergg
    );
    eqsys->SetEllipticalRadialGridGenerator(ergg);

    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms =
        new DREAM::OtherQuantityHandler::eqn_terms;
    struct OtherQuantityHandler::eqn_terms *stream_terms =
        new OtherQuantityHandler::eqn_terms;

    // Timing information
    eqsys->SetTiming(s->GetBool("/output/timingstdout"), s->GetBool("/output/timingfile"));

    // Initialize from previous simulation output?
    const real_t t0 = DREAM::SimulationGenerator::ConstructInitializer(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(
        eqsys, s
    );

    // Construct equations according to settings
    ConstructEquations(eqsys, s, adas, amjuel, nist, oqty_terms, stream_terms);

    // Construct the "other" quantity handler
    ConstructOtherQuantityHandler(eqsys, s, stream_terms, oqty_terms);

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
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms,
    struct OtherQuantityHandler::eqn_terms *stream_terms
) {
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    DREAM::FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    DREAM::FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();    

    enum DREAM::OptionConstants::momentumgrid_type ht_type = eqsys->GetHotTailGridType();
    enum DREAM::OptionConstants::momentumgrid_type re_type = eqsys->GetRunawayGridType();
    
    // TODO
    ConstructEquation_Ions(eqsys, s, adas, amjuel, stream_terms);
    
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
    
    real_t I_ref = s->GetReal("radialgrid/Iref");
    RunawayElectronConfinementTime *rect = 
        new RunawayElectronConfinementTime(
            eqsys->GetUnknownHandler(), eqsys->GetEllipticalRadialGridGenerator(),
            eqsys->GetConnectionLength(), I_ref
        );
    eqsys->SetRunawayElectronConfinementTime(rect);
    
    // Hot electron quantities
    if (eqsys->HasHotTailGrid()) {
        ConstructEquation_f_hot(eqsys, s, oqty_terms, stream_terms);
    }

    // Runaway electron quantities
    if (eqsys->HasRunawayGrid()) {
        ConstructEquation_f_re(eqsys, s, oqty_terms, stream_terms);
    }

    // Standard equations
    ConstructEquation_E_field(eqsys, s, oqty_terms);
    DREAM::SimulationGenerator::ConstructEquation_j_hot(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_tot(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_ohm(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_j_re(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_n_cold(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_n_hot(eqsys, s);
    ConstructEquation_T_cold(eqsys, s, adas, amjuel, nist, oqty_terms, stream_terms);
    ConstructEquation_lambda_i(eqsys, s, adas); 

    enum DREAM::OptionConstants::uqty_T_i_eqn typeTi =
        (DREAM::OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if (typeTi != DREAM::OptionConstants::UQTY_T_I_INCLUDE)
        throw DREAM::SettingsException(
            "T_i not included: STREAM requires the ion temperatures to be evolved."
        );
    
    DREAM::SimulationGenerator::ConstructEquation_Ion_Ni(eqsys, s);
    ConstructEquation_T_i(eqsys, s, adas, stream_terms);

    // NOTE: The runaway number density may depend explicitly on
    // the hot-tail equation and must therefore be constructed
    // AFTER the call to 'ConstructEquation_f_hot()'.
    ConstructEquation_n_re(eqsys, s, oqty_terms);

    // Helper quantities
    DREAM::SimulationGenerator::ConstructEquation_n_tot(eqsys, s);
    enum DREAM::OptionConstants::eqterm_hottail_mode ht_mode =
        (enum DREAM::OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    enum DREAM::OptionConstants::uqty_f_hot_dist_mode ht_dist_mode =
        (enum DREAM::OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");
    if (ht_mode != DREAM::OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED &&
        ht_dist_mode == DREAM::OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL)
        DREAM::SimulationGenerator::ConstructEquation_tau_coll(eqsys);
      
    eqsys->GetConnectionLength()->Initialize();
}

/**
 * Construct the unknowns of the STREAM equation system.
 */
void SimulationGenerator::ConstructUnknowns(
    EquationSystem *eqsys, DREAM::Settings *s
) {
    DREAM::FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    DREAM::FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();

    #define DEFU_HOT(NAME) eqsys->SetUnknown( \
        DREAM::OptionConstants::UQTY_ ## NAME, \
        DREAM::OptionConstants::UQTY_ ## NAME ## _DESC, \
        hottailGrid)
    #define DEFU_RE(NAME) eqsys->SetUnknown( \
        DREAM::OptionConstants::UQTY_ ## NAME, \
        DREAM::OptionConstants::UQTY_ ## NAME ## _DESC, \
        runawayGrid)
    #define DEFU_FLD(NAME) eqsys->SetUnknown( \
        DREAM::OptionConstants::UQTY_ ## NAME, \
        DREAM::OptionConstants::UQTY_ ## NAME ## _DESC, \
        fluidGrid)
    #define DEFU_FLD_N(NAME,NMULT) eqsys->SetUnknown( \
        DREAM::OptionConstants::UQTY_ ## NAME, \
        DREAM::OptionConstants::UQTY_ ## NAME ## _DESC, \
        fluidGrid, (NMULT))
    #define DEFU_SCL(NAME) eqsys->SetUnknown( \
        DREAM::OptionConstants::UQTY_ ## NAME, \
        DREAM::OptionConstants::UQTY_ ## NAME ## _DESC, \
        scalarGrid)
    #define DEFU_SCL_N(NAME,NMULT) eqsys->SetUnknown( \
        DREAM::OptionConstants::UQTY_ ## NAME, \
        DREAM::OptionConstants::UQTY_ ## NAME ## _DESC, \
        scalarGrid,(NMULT))

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        DEFU_HOT(F_HOT);
    }

    // Runaway quantities
    if (runawayGrid != nullptr) {
        DEFU_RE(F_RE);
    }

    // Fluid quantities
    len_t nIonChargeStates = DREAM::SimulationGenerator::GetNumberOfIonChargeStates(s);
    DEFU_FLD_N(ION_SPECIES, nIonChargeStates);
    DEFU_FLD(N_HOT);
    DEFU_FLD(N_COLD);
    DEFU_FLD(N_RE);
    DEFU_FLD(J_OHM);
    DEFU_FLD(J_HOT);
    DEFU_FLD(J_RE);
    DEFU_FLD(J_TOT);
    DEFU_FLD(T_COLD);
    DEFU_FLD(W_COLD);
    DEFU_FLD(E_FIELD);
    DEFU_SCL(I_P);

    enum DREAM::OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode = (enum DREAM::OptionConstants::eqterm_spi_ablation_mode)s->GetInteger("eqsys/spi/ablation");
    if(spi_ablation_mode!=DREAM::OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT){
        len_t nShard;
        s->GetRealArray("eqsys/spi/init/rp", 1, &nShard);
        DEFU_SCL_N(Y_P,nShard);
        DEFU_SCL_N(X_P,3*nShard);
        DEFU_SCL_N(V_P,3*nShard);
        
        if (hottailGrid != nullptr){
        	DEFU_FLD(Q_HOT);
        	DEFU_FLD(W_HOT);
    	}
    	if(spi_ablation_mode==DREAM::OptionConstants::EQTERM_SPI_ABLATION_MODE_NGPS){
    		DEFU_FLD_N(ION_SPECIES_ABL, nIonChargeStates);
    		DEFU_FLD(N_ABL);
    		DEFU_FLD(T_ABL);
    		DEFU_FLD(W_ABL);
		}
    }

    len_t nIonSpecies = DREAM::SimulationGenerator::GetNumberOfIonSpecies(s);
    if( (DREAM::OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi") == DREAM::OptionConstants::UQTY_T_I_INCLUDE ){
        DEFU_FLD_N(WI_ENER, nIonSpecies);
        DEFU_FLD_N(NI_DENS, nIonSpecies);
    }
    
    // Fluid helper quantities
    DEFU_FLD(N_TOT);
    if (hottailGrid != nullptr){
        DEFU_FLD(S_PARTICLE);
    }
    DREAM::OptionConstants::eqterm_hottail_mode hottail_mode = (enum DREAM::OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    DREAM::OptionConstants::uqty_f_hot_dist_mode ht_dist_mode = (enum DREAM::OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");    
    if(hottail_mode != DREAM::OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED && ht_dist_mode == DREAM::OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL){
        DEFU_FLD(TAU_COLL);
    }

    // Mean-free path
    eqsys->SetUnknown(
        OptionConstants::UQTY_LAMBDA_I,
        OptionConstants::UQTY_LAMBDA_I_DESC,
        fluidGrid,
        nIonSpecies
    );
}

