/**
 * Implementation of a STREAM extension to the OtherQuantityHandler.
 */

#include "STREAM/OtherQuantityHandler.hpp"


using namespace STREAM;


/**
 * Constructor.
 */
OtherQuantityHandler::OtherQuantityHandler(
    PlasmaVolume *plasmaVolume,
    // Carried over from DREAM...
    DREAM::CollisionQuantityHandler *cqtyHottail, DREAM::CollisionQuantityHandler *cqtyRunaway,
    DREAM::PostProcessor *postProcessor, DREAM::RunawayFluid *REFluid, DREAM::FVM::UnknownQuantityHandler *unknowns,
    std::vector<DREAM::UnknownQuantityEquation*> *unknown_equations, DREAM::IonHandler *ions,
    DREAM::FVM::Grid *fluidGrid, DREAM::FVM::Grid *hottailGrid, DREAM::FVM::Grid *runawayGrid,
    DREAM::FVM::Grid *scalarGrid, struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) : DREAM::OtherQuantityHandler(cqtyHottail, cqtyRunaway, postProcessor, REFluid,
        unknowns, unknown_equations, ions, fluidGrid, hottailGrid, runawayGrid,
        scalarGrid, oqty_terms),
    plasmaVolume(plasmaVolume) {

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
    /*const len_t nr_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetNr());
    const len_t n1_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_ht = (this->hottailGrid==nullptr ? 0 : this->hottailGrid->GetMomentumGrid(0)->GetNp2());

    const len_t nr_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetNr());
    const len_t n1_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp1());
    const len_t n2_re = (this->runawayGrid==nullptr ? 0 : this->runawayGrid->GetMomentumGrid(0)->GetNp2());*/

    if(hottailGrid != nullptr)
        kineticVectorHot = new real_t[hottailGrid->GetNCells()];
    if(runawayGrid != nullptr)
        kineticVectorRE  = new real_t[runawayGrid->GetNCells()];

    // HELPER MACROS (to make definitions more compact)
    // Define on scalar grid
    #define DEF_SC(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), scalarGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_SC_MUL(NAME, MUL, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), scalarGrid, (MUL), DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](DREAM::FVM::QuantityData *qd) {FUNC}));

    // Define on fluid grid
    #define DEF_FL(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), fluidGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_FL_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), fluidGrid, 1, DREAM::FVM::FLUXGRIDTYPE_RADIAL, [this](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_FL_MUL(NAME, MUL, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), fluidGrid, (MUL), DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this](DREAM::FVM::QuantityData *qd) {FUNC}));

    // Define on hot-tail grid
    #define DEF_HT(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_ht,n1_ht,n2_ht](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_HT_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_RADIAL, [this,nr_ht,n1_ht,n2_ht](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_HT_F1(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P1, [this,nr_ht,n1_ht,n2_ht](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_HT_F2(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), hottailGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P2, [this,nr_ht,n1_ht,n2_ht](DREAM::FVM::QuantityData *qd) {FUNC}));

    // Define on runaway grid
    #define DEF_RE(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION, [this,nr_re,n1_re,n2_re](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_RE_FR(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_RADIAL, [this,nr_re,n1_re,n2_re](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_RE_F1(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P1, [this,nr_re,n1_re,n2_re](DREAM::FVM::QuantityData *qd) {FUNC}));
    #define DEF_RE_F2(NAME, DESC, FUNC) \
        this->all_quantities.push_back(new DREAM::OtherQuantity((NAME), (DESC), runawayGrid, 1, DREAM::FVM::FLUXGRIDTYPE_P2, [this,nr_re,n1_re,n2_re](DREAM::FVM::QuantityData *qd) {FUNC}));

    const len_t nIons = this->ions->GetNZ();
    DEF_SC_MUL("stream/V_n", nIons, "Plasma volume occupied by neutrals",
        const len_t nZ = this->ions->GetNZ();
        real_t *v = qd->StoreEmpty();
        for (len_t i = 0; i < nZ; i++)
            v[i] = this->plasmaVolume->GetNeutralVolume(i);
    );
    DEF_SC_MUL("stream/V_n_tot", nIons, "Total volume occupied by neutrals (including outside plasma)",
        const len_t nZ = this->ions->GetNZ();
        real_t *v = qd->StoreEmpty();
        for (len_t i = 0; i < nZ; i++)
            v[i] = this->plasmaVolume->GetTotalNeutralVolume(i);
    );
    DEF_SC("stream/V_p", "Plasma volume", 
        real_t v = this->plasmaVolume->GetPlasmaVolume();
        qd->Store(&v);
    );
    DEF_SC("stream/V", "Tokamak vessel volume",
        real_t v = this->plasmaVolume->GetVesselVolume();
        qd->Store(&v);
    );

    for (auto qty : all_quantities) {
        if (qty->GetName().substr(0, 6) == "stream")
            this->groups["stream"].push_back(qty->GetName());
    }
}

