/**
 * Implementation of a STREAM extension to the OtherQuantityHandler.
 */

#include "STREAM/OtherQuantityHandler.hpp"


using namespace STREAM;


/**
 * Constructor. 
 */
OtherQuantityHandler::OtherQuantityHandler(
    ConnectionLength *connectionLength, ConfinementTime *confinementTime, NeutralInflux *neutralInflux,
    PlasmaVolume *plasmaVolume, RunawayElectronConfinementTime *rect, 
    OpticalThickness *opticalThickness,
    std::vector<IonRateEquation*> ionRateEquations, struct eqn_terms *stream_terms,
    // Carried over from DREAM...
    DREAM::CollisionQuantityHandler *cqtyHottail, DREAM::CollisionQuantityHandler *cqtyRunaway,
    DREAM::PostProcessor *postProcessor, DREAM::RunawayFluid *REFluid, DREAM::FVM::UnknownQuantityHandler *unknowns,
    std::vector<DREAM::UnknownQuantityEquation*> *unknown_equations, DREAM::IonHandler *ions,
    DREAM::FVM::Grid *fluidGrid, DREAM::FVM::Grid *hottailGrid, DREAM::FVM::Grid *runawayGrid,
    DREAM::FVM::Grid *scalarGrid, struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) : DREAM::OtherQuantityHandler(cqtyHottail, cqtyRunaway, postProcessor, REFluid,
        unknowns, unknown_equations, ions, nullptr, fluidGrid, hottailGrid, runawayGrid,
        scalarGrid, oqty_terms),
    connectionLength(connectionLength), confinementTime(confinementTime), neutralInflux(neutralInflux), plasmaVolume(plasmaVolume),
    reConfinementTime(rect), opticalThickness(opticalThickness), ionRateEquations(ionRateEquations), stream_terms(stream_terms) {

    this->id_ni = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES);

    this->DefineQuantitiesSTREAM();
}

/**
 * Destructor.
 */
OtherQuantityHandler::~OtherQuantityHandler() {
}

/**
 * Define STREAM-specific other quantities.
 */
void OtherQuantityHandler::DefineQuantitiesSTREAM() {
    // XXX here we assume that all momentum grids are the same
    const len_t nr = this->fluidGrid->GetNr();
    const len_t nr_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetNr());
    const len_t n1_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp2());

    /*const len_t nr_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetNr());
    const len_t n1_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp2());*/

    if(hottailGrid != nullptr)
        kineticVectorHot = new real_t[hottailGrid->GetNCells()];
    if(runawayGrid != nullptr)
        kineticVectorRE  = new real_t[runawayGrid->GetNCells()];

    // HELPER MACROS (to make definitions more compact)
    // Define on scalar grid
    #define DEF_SC(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), scalarGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr](const real_t t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_SC_MUL(NAME, MUL, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), scalarGrid, (MUL), DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr](const real_t t, DREAM::FVM::QuantityData *qd) {FUNC}));

    // Define on fluid grid
    #define DEF_FL(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), fluidGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_FL_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), fluidGrid, 1, DREAM::FVM::FLUXGRIDTYPE_RADIAL, [this,nr](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_FL_MUL(NAME, MUL, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), fluidGrid, (MUL), DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));

    // Define on hot-tail grid
    #define DEF_HT(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_ht,n1_ht,n2_ht](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_HT_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_RADIAL, [this,nr_ht,n1_ht,n2_ht](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_HT_F1(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P1, [this,nr_ht,n1_ht,n2_ht](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_HT_F2(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P2, [this,nr_ht,n1_ht,n2_ht](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));

    // Define on runaway grid
    #define DEF_RE(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_re,n1_re,n2_re](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_RE_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_RADIAL, [this,nr_re,n1_re,n2_re](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_RE_F1(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P1, [this,nr_re,n1_re,n2_re](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_RE_F2(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P2, [this,nr_re,n1_re,n2_re](const real_t, DREAM::FVM::QuantityData *qd) {FUNC}));

    const len_t nIons = this->ions->GetNZ();
    const len_t nChargeStates = this->ions->GetNzs();

    if (this->stream_terms->iontransport != nullptr) {
        DEF_SC_MUL("stream/ni_iontransport", nChargeStates, "Ion particle transport rate for each species",
            const len_t nZ = this->ions->GetNZ();
            const len_t nZs = this->ions->GetNzs();
            real_t *v = qd->StoreEmpty();

            for (len_t i = 0; i < nZs; i++)
                v[i] = 0;

            for (len_t iZ = 0; iZ < nZ; iZ++) {
                if (this->stream_terms->iontransport[iZ] != nullptr)
                    this->stream_terms->iontransport[iZ]->SetVectorElements(v, nullptr);
            }
        );
    }

    if (this->stream_terms->Tcold_transport != nullptr) {
        DEF_SC("stream/Tcold_transport", "Electron heat transport (from STREAM::ElectronHeatTransport)",
            real_t *v = qd->StoreEmpty();
            real_t *Wcold = this->unknowns->GetUnknownData(id_Wcold);
            v[0] = 0;
            this->stream_terms->Tcold_transport->SetVectorElements(v, Wcold);
        );
    }
    
    if (this->stream_terms->Tcold_ECH != nullptr) {
        DEF_SC("fluid/Tcold_ECH", "Electron Cyclotron Heating (from STREAM::ElectronCyclotronHeating)",
            real_t *v = qd->StoreEmpty();
            v[0] = 0;
            this->stream_terms->Tcold_ECH->SetVectorElements(v, nullptr);
        );
    }

    if (this->stream_terms->Wi_e_coll != nullptr) {
        DEF_SC_MUL("stream/Wi_e_coll", nIons, "Ion-electron collision heat loss rate",
            const len_t nZ = this->ions->GetNZ();
            const real_t *nions = this->unknowns->GetUnknownData(this->id_ni);
            real_t *v = qd->StoreEmpty();

            for (len_t i = 0; i < nZ; i++)
                v[i] = 0;

            for (len_t iZ = 0; iZ < nZ; iZ++) {
                if (this->stream_terms->Wi_e_coll[iZ] != nullptr)
                    this->stream_terms->Wi_e_coll[iZ]->SetVectorElements(v, nions);
            }
        );
    }

    if (this->stream_terms->Wi_iontransport != nullptr) {
        DEF_SC_MUL("stream/Wi_iontransport", nIons, "Ion heat transport rate for each species",
            const len_t nZ = this->ions->GetNZ();
            const real_t *nions = this->unknowns->GetUnknownData(this->id_ni);
            real_t *v = qd->StoreEmpty();

            for (len_t i = 0; i < nZ; i++)
                v[i] = 0;

            for (len_t iZ = 0; iZ < nZ; iZ++) {
                if (this->stream_terms->Wi_iontransport[iZ] != nullptr)
                    this->stream_terms->Wi_iontransport[iZ]->SetVectorElements(v, nions);
            }
        );
    }

    if (this->stream_terms->Wi_chargeexchange != nullptr) {
        DEF_SC_MUL("stream/Wi_chargeexchange", nIons, "Ion energy loss due to charge-exchange",
            const len_t nZ = this->ions->GetNZ();
            const real_t *nions = this->unknowns->GetUnknownData(this->id_ni);
            real_t *v = qd->StoreEmpty();

            for (len_t i = 0; i < nZ; i++)
                v[i] = 0;

            for (len_t iZ = 0; iZ < nZ; iZ++) { 
                if (this->stream_terms->Wi_chargeexchange[iZ] != nullptr)
                    this->stream_terms->Wi_chargeexchange[iZ]->SetVectorElements(v, nions);
            }
        );
    }

    DEF_SC_MUL("stream/neutralinflux", nIons, "Influx rate of neutral particles of each species",
        const len_t nZ = this->ions->GetNZ();
        real_t *v = qd->StoreEmpty();
        for (len_t iz = 0; iz < nZ; iz++)
            v[iz] = this->neutralInflux->EvaluateNeutralInflux(iz, t);
    );

    DEF_FL("stream/tau_D", "Deuterium confinement time",
        real_t *v = qd->StoreEmpty();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = 1/this->confinementTime->EvaluateConfinementTime(ir);
    );
    DEF_FL("stream/tau_D_par", "Parallel deuterium confinement time",
        real_t *v = qd->StoreEmpty();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = 1/this->confinementTime->EvaluateParallelConfinementTime(ir);
    );
    DEF_FL("stream/tau_D_perp", "Perpendicular deuterium confinement time",
        real_t *v = qd->StoreEmpty();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = 1/this->confinementTime->EvaluatePerpendicularConfinementTime(ir);
    );

    DEF_FL("stream/tau_RE", "Runaway electron confinement time [s]",
        real_t *v = qd->StoreEmpty();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = 1/this->reConfinementTime->EvaluateInverse(ir);
    );
    DEF_FL("stream/tau_RE1", "Runaway electron confinement time (component dominating at low I_p) [s]",
        real_t *v = qd->StoreEmpty();
        const real_t I_p = this->unknowns->GetUnknownData(id_Ip)[0];
        const real_t I_ref = this->reConfinementTime->GetIref();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = exp(I_p/I_ref) * this->reConfinementTime->EvaluateRunawayElectronConfinementTime1(ir);
    );
    DEF_FL("stream/tau_RE2", "Runaway electron confinement time (component dominating at high I_p) [s]",
        real_t *v = qd->StoreEmpty();
        const real_t I_p = this->unknowns->GetUnknownData(id_Ip)[0];
        const real_t I_ref = this->reConfinementTime->GetIref();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = 1/(1-exp(-I_p/I_ref)) * this->reConfinementTime->EvaluateRunawayElectronConfinementTime2(ir);
    );
    
    
    DEF_FL("stream/Lf", "Effective distance travelled by a particle before colliding with the wall [m]",
        real_t *v = qd->StoreEmpty();
        for (len_t ir = 0; ir < nr; ir++)
            v[ir] = this->connectionLength->EvaluateConnectionLength(ir);
    );
    
	if (this->opticalThickness) {
		DEF_FL("stream/eta_o", "Optical thickness for O mode",
			real_t *v = qd->StoreEmpty();
			for (len_t ir = 0; ir < nr; ir++)
				v[ir] = this->opticalThickness->EvaluateOpticalThickness_o(ir);
		);
		DEF_FL("stream/eta_x", "Optical thickness for X mode",
			real_t *v = qd->StoreEmpty();
			for (len_t ir = 0; ir < nr; ir++)
				v[ir] = this->opticalThickness->EvaluateOpticalThickness_x(ir);
		);
	}

	DEF_SC("stream/a", "Plasma minor radius [m]",
		real_t a = this->fluidGrid->GetRadialGrid()->GetDr(0);
		qd->Store(&a);
	);
	DEF_SC("stream/A", "Plasma poloidal cross-section [m^2]",
		real_t A = this->plasmaVolume->GetPlasmaCrossSection();
		qd->Store(&A);
	);
    DEF_SC_MUL("stream/V_n", nIons, "Plasma volume occupied by neutrals [m^3]",
        const len_t nZ = this->ions->GetNZ();
        real_t *v = qd->StoreEmpty();
        for (len_t i = 0; i < nZ; i++)
            v[i] = this->plasmaVolume->GetNeutralVolume(i);
    );
    DEF_SC_MUL("stream/V_n_tot", nIons, "Total volume occupied by neutrals (including outside plasma) [m^3]",
        const len_t nZ = this->ions->GetNZ();
        real_t *v = qd->StoreEmpty();
        for (len_t i = 0; i < nZ; i++)
            v[i] = this->plasmaVolume->GetTotalNeutralVolume(i);
    );
    DEF_SC("stream/V_p", "Plasma volume [m^3]", 
        real_t v = this->plasmaVolume->GetPlasmaVolume();
        qd->Store(&v);
    );
    DEF_SC("stream/V", "Tokamak vessel volume [m^3]",
        real_t v = this->plasmaVolume->GetVesselVolume();
        qd->Store(&v);
    );
    
    if (this->stream_terms->DPT != nullptr)
		DEF_HT("hottail/parallel_transport", "Parallel transport",
			real_t *v = qd->StoreEmpty();

			DistributionParallelTransport *DPT = stream_terms->DPT;
			real_t *f = this->unknowns->GetUnknownData(id_f_hot);
			DPT->SetVectorElements(v, f);
		);

    // Diagnostics for ion rate equations
    DEF_FL_MUL("stream/ionrateequation_posIonization", nChargeStates, "Positive ionization term in ion rate equation",
        real_t *v = qd->StoreEmpty();
        len_t offset = 0;
        for (len_t iz = 0; iz < this->ionRateEquations.size(); iz++) {
            IonRateEquation *ire = this->ionRateEquations[iz];
            len_t Z = ire->GetZ();

            real_t **t = ire->GetPositiveIonizationTerm();
            for (len_t Z0 = 0; Z0 <= Z; Z0++)
                for (len_t ir = 0; ir < nr; ir++)
                    v[offset+Z0*nr+ir] = t[Z0][ir];
            offset+=(Z+1)*nr;
        }
    );

    // Diagnostics for ion rate equations
    DEF_FL_MUL("stream/ionrateequation_negIonization", nChargeStates, "Negative ionization term in ion rate equation",
        real_t *v = qd->StoreEmpty();
        len_t offset = 0;
        for (len_t iz = 0; iz < this->ionRateEquations.size(); iz++) {
            IonRateEquation *ire = this->ionRateEquations[iz];
            len_t Z = ire->GetZ();

            real_t **t = ire->GetNegativeIonizationTerm();
            for (len_t Z0 = 0; Z0 <= Z; Z0++)
                for (len_t ir = 0; ir < nr; ir++)
                    v[offset+Z0*nr+ir] = t[Z0][ir];
            offset+=(Z+1)*nr;
        }
    );

    // Diagnostics for ion rate equations
    DEF_FL_MUL("stream/ionrateequation_posRecombination", nChargeStates, "Positive recombination term in ion rate equation",
        real_t *v = qd->StoreEmpty();
        len_t offset = 0;
        for (len_t iz = 0; iz < this->ionRateEquations.size(); iz++) {
            IonRateEquation *ire = this->ionRateEquations[iz];
            len_t Z = ire->GetZ();

            real_t **t = ire->GetPositiveRecombinationTerm();
            for (len_t Z0 = 0; Z0 <= Z; Z0++)
                for (len_t ir = 0; ir < nr; ir++)
                    v[offset+Z0*nr+ir] = t[Z0][ir];
            offset+=(Z+1)*nr;
        }
    );

    // Diagnostics for ion rate equations
    DEF_FL_MUL("stream/ionrateequation_negRecombination", nChargeStates, "Negative recombination term in ion rate equation",
        real_t *v = qd->StoreEmpty();
        len_t offset = 0;
        for (len_t iz = 0; iz < this->ionRateEquations.size(); iz++) {
            IonRateEquation *ire = this->ionRateEquations[iz];
            len_t Z = ire->GetZ();

            real_t **t = ire->GetNegativeRecombinationTerm();
            for (len_t Z0 = 0; Z0 <= Z; Z0++)
                for (len_t ir = 0; ir < nr; ir++)
                    v[offset+Z0*nr+ir] = t[Z0][ir];
            offset+=(Z+1)*nr;
        }
    );

    // Diagnostics for ion rate equations
    DEF_FL_MUL("stream/ionrateequation_posChargeExchange", nChargeStates, "Positive charge-exchange term in ion rate equation",
        real_t *v = qd->StoreEmpty();
        len_t offset = 0;
        for (len_t iz = 0; iz < this->ionRateEquations.size(); iz++) {
            IonRateEquation *ire = this->ionRateEquations[iz];
            len_t Z = ire->GetZ();

            real_t **t = ire->GetPositiveChargeExchangeTerm();
            for (len_t Z0 = 0; Z0 <= Z; Z0++)
                for (len_t ir = 0; ir < nr; ir++)
                    v[offset+Z0*nr+ir] = t[Z0][ir];
            offset+=(Z+1)*nr;
        }
    );

    // Diagnostics for ion rate equations
    DEF_FL_MUL("stream/ionrateequation_negChargeExchange", nChargeStates, "Negative charge-exchange term in ion rate equation",
        real_t *v = qd->StoreEmpty();
        len_t offset = 0;
        for (len_t iz = 0; iz < this->ionRateEquations.size(); iz++) {
            IonRateEquation *ire = this->ionRateEquations[iz];
            len_t Z = ire->GetZ();

            real_t **t = ire->GetNegativeChargeExchangeTerm();
            for (len_t Z0 = 0; Z0 <= Z; Z0++)
                for (len_t ir = 0; ir < nr; ir++)
                    v[offset+Z0*nr+ir] = t[Z0][ir];
            offset+=(Z+1)*nr;
        }
    );


    for (auto qty : all_quantities) {
        if (qty->GetName().substr(0, 6) == "stream")
            this->groups["stream"].push_back(qty->GetName());
    }
}

