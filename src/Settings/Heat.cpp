/**
 * Configuration of the heat dynamics in STREAM.
 */

#include <string>
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesIdentityTerm.hpp"
#include "STREAM/Equations/IonHeatTransport.hpp"
#include "STREAM/Settings/OptionConstants.hpp"


using namespace std;
using namespace STREAM;

#define MODULENAME "eqsys/T_cold"


/**
 * Define options for the electron temperature module.
 */
void SimulationGenerator::DefineOptions_T_cold(DREAM::Settings *s) {
    s->DefineSetting(
        MODULENAME "/type",
        "Type of equation to use for determining the electron temperature evolution",
        (int_t)DREAM::OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED
    );
    
    // Prescribed data
    DREAM::SimulationGenerator::DefineDataRT(MODULENAME, s, "data");
    // Prescribed initial profile (when evolving T_cold self-consistently)
    DREAM::SimulationGenerator::DefineDataR(MODULENAME, s, "init");

    // TODO transport settings
}

/**
 * Construct the equation for the temperature.
 */
void SimulationGenerator::ConstructEquation_T_cold(
    EquationSystem *eqsys, DREAM::Settings *s,
    DREAM::ADAS *adas, DREAM::AMJUEL *amjuel, DREAM::NIST *nist,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    enum DREAM::OptionConstants::uqty_T_cold_eqn type =
        (enum DREAM::OptionConstants::uqty_T_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case DREAM::OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            DREAM::SimulationGenerator::ConstructEquation_T_cold_prescribed(eqsys, s);
            break;

        case DREAM::OptionConstants::UQTY_T_COLD_SELF_CONSISTENT:
            ConstructEquation_T_cold_selfconsistent(eqsys, s, adas, amjuel, nist, oqty_terms);
            break;

        default:
            throw DREAM::SettingsException(
                "Unrecognized equation type '%s': %d.",
                DREAM::OptionConstants::UQTY_T_COLD, type
            );
    }
}

/**
 * Construct the equation for a self-consistent temperature evolution.
 */
void SimulationGenerator::ConstructEquation_T_cold_selfconsistent(
    EquationSystem *eqsys, DREAM::Settings *s,
    DREAM::ADAS *adas, DREAM::AMJUEL *amjuel, DREAM::NIST *nist,
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::IonHandler *ionHandler = eqsys->GetIonHandler();
    DREAM::FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    len_t id_T_cold  = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    len_t id_W_cold  = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_W_COLD);
    len_t id_n_cold  = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    len_t id_E_field = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);

    DREAM::FVM::Operator *op_W_cold  = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *op_E_field = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *op_n_cold  = new DREAM::FVM::Operator(fluidGrid);

    // Add dW/dt term
    op_W_cold->AddTerm(new DREAM::FVM::TransientTerm(fluidGrid, id_W_cold));

    // Add ohmic heating term (j*E)
    oqty_terms->T_cold_ohmic = new DREAM::OhmicHeatingTerm(fluidGrid, unknowns);
    op_E_field->AddTerm(oqty_terms->T_cold_ohmic);

    // Gather opacity settings
    len_t ni_types;
    const int_t *iopacity_modes = s->GetIntegerArray(MODULENAME "/opacity_modes", 1, &ni_types);
    enum DREAM::OptionConstants::ion_opacity_mode *opacity_mode =
        new enum DREAM::OptionConstants::ion_opacity_mode[ni_types];
    for (len_t i = 0; i < ni_types; i++)
        opacity_mode[i] = (enum DREAM::OptionConstants::ion_opacity_mode)iopacity_modes[i];

    // Add radiation loss term (TODO TODO TODO)
    bool withRecombinationRadiation = s->GetBool(MODULENAME "/recombination");
    oqty_terms->T_cold_radiation = new DREAM::RadiatedPowerTerm(
        fluidGrid, unknowns, ionHandler, adas, nist,
        amjuel, opacity_mode, withRecombinationRadiation
    );
    op_n_cold->AddTerm(oqty_terms->T_cold_radiation);

    // Add transport (TODO TODO TODO)
    bool hasTransport = false;      // TODO Load from settings...
    //op_W_cold->AddTerm(...)

    eqsys->SetOperator(id_T_cold, id_E_field, op_E_field);
    eqsys->SetOperator(id_T_cold, id_n_cold, op_n_cold);
    string desc = "dWc/dt = j_ohm*E - radiation";

    if (hasTransport) {
        oqty_terms->T_cold_transport = op_W_cold->GetAdvectionDiffusion();
        eqsys->SetOperator(id_T_cold, id_W_cold, op_W_cold);
        desc += " - transport";
    }

    // Energy transfer from runaways to cold electrons.
    // If the kinetic runaway grid is enabled and we do not resolve the cold
    // electrons kinetically...
    enum DREAM::OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum DREAM::OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if (eqsys->HasRunawayGrid() && collfreq_mode != DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) {
        len_t id_f_re = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_F_RE);

        oqty_terms->T_cold_fre_coll = new DREAM::CollisionalEnergyTransferKineticTerm(
            fluidGrid, eqsys->GetRunawayGrid(),
            id_T_cold, id_f_re, eqsys->GetRunawayCollisionHandler(), eqsys->GetUnknownHandler(),
            eqsys->GetRunawayGridType(), -1.0
        );

        DREAM::FVM::Operator *op_f_re = new DREAM::FVM::Operator(fluidGrid);
        op_f_re->AddTerm(oqty_terms->T_cold_nre_coll);
        eqsys->SetOperator(id_T_cold, id_f_re, op_f_re);

        desc += " + int(W*nu_E*f_re)";
    // ...otherwise, add contribution from fluid runaways.
    } else {
        len_t id_n_re = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_RE);
        oqty_terms->T_cold_nre_coll = new DREAM::CollisionalEnergyTransferREFluidTerm(
            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid()->GetLnLambda(), -1.0
        );

        DREAM::FVM::Operator *op_n_re = new DREAM::FVM::Operator(fluidGrid);
        op_n_re->AddTerm(oqty_terms->T_cold_nre_coll);
        eqsys->SetOperator(id_T_cold, id_n_re, op_n_re);

        desc += " + e*c*Ec*n_re";
    }

    // Collisional energy transfer with ions
    enum DREAM::OptionConstants::uqty_T_i_eqn Ti_type =
        (enum DREAM::OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if (Ti_type == DREAM::OptionConstants::UQTY_T_I_INCLUDE) {
        DREAM::CoulombLogarithm *lnLambda = eqsys->GetREFluid()->GetLnLambda();
        const len_t nZ = ionHandler->GetNZ();
        const len_t id_Wi = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);

        oqty_terms->T_cold_ion_coll = new DREAM::FVM::Operator(fluidGrid);
        for (len_t iz = 0; iz < nZ; iz++) {
            oqty_terms->T_cold_ion_coll->AddTerm(
                new DREAM::MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    0, false,
                    iz, true,
                    unknowns, lnLambda, ionHandler, -1.0
                )
            );
        }

        eqsys->SetOperator(id_T_cold, id_Wi, oqty_terms->T_cold_ion_coll);
        desc += " + sum_i Q_ei";
    }

    eqsys->SetOperator(id_T_cold, id_W_cold, op_W_cold, desc);

    // Initialize T_cold
    // (if the input temperature profile is not explicitly set, then 'SetInitialValue()' is
    // called with a null-pointer which results in T=0 at t=0)
    real_t *Tcold_init = DREAM::SimulationGenerator::LoadDataR(MODULENAME, fluidGrid->GetRadialGrid(), s, "init");
    eqsys->SetInitialValue(id_T_cold, Tcold_init);
    delete [] Tcold_init;

    DREAM::SimulationGenerator::ConstructEquation_W_cold(eqsys, s);
}

void SimulationGenerator::ConstructEquation_T_i_selfconsistent(EquationSystem *eqsys, Settings* /*s*/){
    const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER); 
    const len_t id_Wcold = eqsys->GetUnknownID(OptionConstants::UQTY_W_COLD);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    IonHandler *ionHandler = eqsys->GetIonHandler();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();    
    const len_t nZ = ionHandler->GetNZ();


    FVM::Operator *Op_Wij = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_Wie = new FVM::Operator(fluidGrid);

    CoulombLogarithm *lnLambda = eqsys->GetREFluid()->GetLnLambda();
    for(len_t iz=0; iz<nZ; iz++){
        Op_Wij->AddTerm( 
            new IonSpeciesTransientTerm(fluidGrid, iz, id_Wi, -1.0)
        );
        for(len_t jz=0; jz<nZ; jz++){
            if(jz==iz) // the term is trivial =0 for self collisions and can be skipped
                continue;
            Op_Wij->AddTerm(
                new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    jz, true,
                    unknowns, lnLambda, ionHandler)
            );
        }
        Op_Wie->AddTerm(
            new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    0, false,
                    unknowns, lnLambda, ionHandler)
        );
        Op_Wij->AddTerm(new IonHeatTransport(eqsys->GetFluidGrid(), ions, iz, eqsys->GetConfinementTime(), eqsys->GetUnknownHandler()));
    }
    eqsys->SetOperator(id_Wi, id_Wi, Op_Wij, "dW_i/dt = sum_j Q_ij + Q_ie");
    eqsys->SetOperator(id_Wi, id_Wcold, Op_Wie);
}

