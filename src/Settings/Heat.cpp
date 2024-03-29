/**
 * Configuration of the heat dynamics in STREAM.
 */

#include <string>
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesIdentityTerm.hpp"
#include "STREAM/Equations/IonHeatTransport.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Fluid/IonSpeciesTransientTerm.hpp"
#include "STREAM/Equations/RadiatedPowerTerm.hpp"
#include "STREAM/Equations/ElectronHeatTransport.hpp"
#include "STREAM/Equations/ElectronHeatTransportDiffusion.hpp"
#include "STREAM/Equations/ChargeExchangeTerm.hpp"

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
    s->DefineSetting(MODULENAME "/recombination", "Whether to include recombination radiation (true) or ionization energy loss (false)", (bool)false);
    
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
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms,
    struct OtherQuantityHandler::eqn_terms *stream_terms
) {
    enum DREAM::OptionConstants::uqty_T_cold_eqn type =
        (enum DREAM::OptionConstants::uqty_T_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case DREAM::OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            DREAM::SimulationGenerator::ConstructEquation_T_cold_prescribed(eqsys, s);
            break;

        case DREAM::OptionConstants::UQTY_T_COLD_SELF_CONSISTENT:
            ConstructEquation_T_cold_selfconsistent(eqsys, s, adas, amjuel, nist, oqty_terms, stream_terms);
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
    struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms,
    struct OtherQuantityHandler::eqn_terms *stream_terms
) {
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::IonHandler *ionHandler = eqsys->GetIonHandler();
    DREAM::FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    EllipticalRadialGridGenerator *radials = eqsys->GetEllipticalRadialGridGenerator();

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
    len_t ni_types = ionHandler->GetNZ();
    enum DREAM::OptionConstants::ion_opacity_mode *opacity_mode =
        new enum DREAM::OptionConstants::ion_opacity_mode[ni_types];
    for (len_t i = 0; i < ni_types; i++)
        opacity_mode[i] = DREAM::OptionConstants::OPACITY_MODE_TRANSPARENT;
    
    PlasmaVolume *pv = eqsys->GetPlasmaVolume();    
    
    // Add radiation loss term 
    bool withRecombinationRadiation = s->GetBool(MODULENAME "/recombination");
    oqty_terms->T_cold_radiation = new STREAM::RadiatedPowerTerm(
        fluidGrid, unknowns, ionHandler, adas, nist,
        amjuel, opacity_mode, withRecombinationRadiation, pv
    );
    op_n_cold->AddTerm(oqty_terms->T_cold_radiation);

    // Add transport
    ElectronHeatTransport *ht = new ElectronHeatTransport(
        eqsys->GetFluidGrid(), eqsys->GetConfinementTime(),
        eqsys->GetEllipticalRadialGridGenerator(),
        eqsys->GetUnknownHandler()
    );
    stream_terms->Tcold_transport = ht;
    op_W_cold->AddTerm(ht);
    
    string desc = "dWc/dt = j_ohm*E - radiation - transport";
    
    // Add ECH
    real_t P_inj = s->GetReal("radialgrid/P_inj"); 
    real_t f_o   = s->GetReal("radialgrid/f_o"); 
    real_t f_x   = s->GetReal("radialgrid/f_x"); 
    real_t theta = s->GetReal("radialgrid/theta"); 
    real_t phi   = s->GetReal("radialgrid/phi");
    len_t N      = s->GetInteger("radialgrid/N");
    
    
    if (P_inj != 0.0) {
        OpticalThickness *opticalThickness = new OpticalThickness(unknowns, radials, N, theta, phi);
        eqsys->SetOpticalThickness(opticalThickness);
    
        ElectronCyclotronHeating *ech = new ElectronCyclotronHeating(
            eqsys->GetFluidGrid(), eqsys->GetEllipticalRadialGridGenerator(),
            eqsys->GetUnknownHandler(), eqsys->GetOpticalThickness(), eqsys->GetPlasmaVolume(), P_inj, f_o, f_x, theta
        );
        stream_terms->Tcold_ECH = ech;
        op_W_cold->AddTerm(ech);
        
        desc += " + ECH";
    }

    eqsys->SetOperator(id_T_cold, id_E_field, op_E_field);
    eqsys->SetOperator(id_T_cold, id_n_cold, op_n_cold);

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

void SimulationGenerator::ConstructEquation_T_i(
    EquationSystem *eqsys, DREAM::Settings *s, DREAM::ADAS *adas,
    STREAM::OtherQuantityHandler::eqn_terms *stream_terms
) {
    /**
     * if the electron heat W_cold is evolved self-consistently,
     * also evolve the ion heat W_i. Otherwise set it to constant.
     */
    enum DREAM::OptionConstants::uqty_T_cold_eqn TcoldType = (enum DREAM::OptionConstants::uqty_T_cold_eqn)s->GetInteger("eqsys/T_cold/type");
    if(TcoldType==DREAM::OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED)
        DREAM::SimulationGenerator::ConstructEquation_T_i_trivial(eqsys, s);
    else if (TcoldType == DREAM::OptionConstants::UQTY_T_COLD_SELF_CONSISTENT)
        ConstructEquation_T_i_selfconsistent(eqsys, s, adas, stream_terms);
    else 
        throw DREAM::SettingsException(
            "T_i: Unrecognized equation type for '%s': %d.",
            DREAM::OptionConstants::UQTY_T_COLD, TcoldType
        );

    // Initialize heat from ion densities and input ion temperatures
    /*real_t *Ti_init = DREAM::SimulationGenerator::LoadDataIonR("eqsys/n_i", eqsys->GetFluidGrid()->GetRadialGrid(), s, eqsys->GetIonHandler()->GetNZ(), "initialTi");
    const real_t *Ni_init = eqsys->GetUnknownHandler()->GetUnknownInitialData(eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS));
    for(len_t it=0; it<eqsys->GetIonHandler()->GetNZ()*eqsys->GetFluidGrid()->GetNr(); it++)
        Ti_init[it] *= 1.5*DREAM::Constants::ec*Ni_init[it];
    eqsys->SetInitialValue(eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER), Ti_init);*/

	real_t *Ti_init = DREAM::SimulationGenerator::LoadDataIonR("eqsys/n_i", eqsys->GetFluidGrid()->GetRadialGrid(), s, eqsys->GetIonHandler()->GetNZ(), "initialTi");
	len_t id_Ni = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);

	std::function<void(DREAM::FVM::UnknownQuantityHandler*, real_t*)> initfunc_Ti =
		[Ti_init,id_Ni,eqsys](DREAM::FVM::UnknownQuantityHandler *u, real_t *Ti) {
		const real_t *Ni_init = u->GetUnknownInitialData(id_Ni);
		for(len_t it=0; it<eqsys->GetIonHandler()->GetNZ()*eqsys->GetFluidGrid()->GetNr(); it++)
			Ti[it] = Ti_init[it] * 1.5*DREAM::Constants::ec*Ni_init[it];
	};
	len_t id_Wi = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
	eqsys->initializer->AddRule(
		id_Wi,
		DREAM::EqsysInitializer::INITRULE_EVAL_FUNCTION,
		initfunc_Ti,
		id_Ni
	);
}

void SimulationGenerator::ConstructEquation_T_i_selfconsistent(
    EquationSystem *eqsys, DREAM::Settings* , DREAM::ADAS *adas,
    struct OtherQuantityHandler::eqn_terms *stream_terms
){
    const len_t id_Wi = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER); 
    const len_t id_Wcold = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_W_COLD);
    const len_t id_ni = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES);

    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    DREAM::IonHandler *ionHandler = eqsys->GetIonHandler();
    DREAM::FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();    
    const len_t nZ = ionHandler->GetNZ();


    DREAM::FVM::Operator *Op_Wij = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_Wie = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_ni = new DREAM::FVM::Operator(fluidGrid);

    stream_terms->Wi_chargeexchange = new ChargeExchangeTerm*[nZ];
    stream_terms->Wi_e_coll = new DREAM::MaxwellianCollisionalEnergyTransferTerm*[nZ];
    stream_terms->Wi_iontransport = new IonHeatTransport*[nZ];

    // Locate deuterium
    len_t D_index;
    const vector<string>& ionNames = ionHandler->GetNameList();
    for (D_index = 0; D_index < ionNames.size(); D_index++)
        if (ionNames[D_index] == "D" || ionNames[D_index] == "H")
            break;

    if (D_index == ionNames.size())
        throw DREAM::SettingsException(
            "T_i: Expected to find one deuterium/hydrogen ion species named 'D'/'H'."
        );

    D_index = ionHandler->GetIndex(D_index, 0);

    DREAM::CoulombLogarithm *lnLambda = eqsys->GetREFluid()->GetLnLambda();
    for(len_t iz=0; iz<nZ; iz++){
        Op_Wij->AddTerm( 
            new DREAM::IonSpeciesTransientTerm(fluidGrid, iz, id_Wi, -1.0)
        );

        // i-i collisions
        for(len_t jz=0; jz<nZ; jz++){
            if(jz==iz) // the term is trivial =0 for self collisions and can be skipped
                continue;
            Op_Wij->AddTerm(
                new DREAM::MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    iz, true,
                    jz, true,
                    unknowns, lnLambda, ionHandler)
            );
        }

        // i-e collisions
        stream_terms->Wi_e_coll[iz] = new DREAM::MaxwellianCollisionalEnergyTransferTerm(
            fluidGrid,
            iz, true,
            0, false,
            unknowns, lnLambda, ionHandler
        );
        Op_Wie->AddTerm(stream_terms->Wi_e_coll[iz]);

        // Heat transport
        stream_terms->Wi_iontransport[iz] = new IonHeatTransport(eqsys->GetFluidGrid(), eqsys->GetIonHandler(), iz, eqsys->GetConfinementTime(), eqsys->GetUnknownHandler(), eqsys->GetEllipticalRadialGridGenerator());
        Op_ni->AddTerm(stream_terms->Wi_iontransport[iz]);

        // Charge exchange
        if (iz == D_index){
            stream_terms->Wi_chargeexchange[iz] = new ChargeExchangeTerm(eqsys->GetFluidGrid(), eqsys->GetUnknownHandler(), eqsys->GetIonHandler(), iz, adas, eqsys->GetPlasmaVolume(), eqsys->GetEllipticalRadialGridGenerator(), fluidGrid, D_index);
            Op_ni->AddTerm(stream_terms->Wi_chargeexchange[iz]);
        } else
            stream_terms->Wi_chargeexchange[iz] = nullptr;
    }
    eqsys->SetOperator(id_Wi, id_Wi, Op_Wij, "dW_i/dt = sum_j Q_ij + Q_ie - Q_CX - Q_transport");
    eqsys->SetOperator(id_Wi, id_Wcold, Op_Wie);
    eqsys->SetOperator(id_Wi, id_ni, Op_ni/*, "dW_i/dt = V_n,i/V_p * [ 3/2 * n_i^(0) * ( T_i - 0.026 ) * ( sum R_j,cx^(1) n_j^(1) ) ]"*/);
}

