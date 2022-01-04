/**
 * Construction of the DYON circuit equations.
 */

#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "STREAM/Settings/OptionConstants.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


using namespace STREAM;


namespace STREAM {
    class VloopTerm : public DREAM::FVM::DiagonalLinearTerm {
    public:
        VloopTerm(DREAM::FVM::Grid* g) : DREAM::FVM::DiagonalLinearTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            DREAM::FVM::RadialGrid *rg = grid->GetRadialGrid();
            for (len_t ir = 0; ir < nr; ir++){
                // E_field is defined as E = <E.B>/sqrt(<B^2>), and Vloop is
                // defined as Vloop = 2*pi * <E.B>/G
                // To go from E_field to Vloop we must therefore multiply E by
                // a factor 2*pi * sqrt(<B^2>) / G
                // (the Bmin comes from the fact that FSA_B2 is normalized to Bmin)
                real_t w = 2*M_PI * (sqrt(rg->GetFSA_B2(ir)) * rg->GetBmin(ir)) / rg->GetBTorG(ir);
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}

#define MODULENAME "eqsys/E_field"

void SimulationGenerator::DefineOptions_ElectricField(DREAM::Settings *s) {
    DREAM::SimulationGenerator::DefineOptions_ElectricField(s);

    // Circuit equation options
    s->DefineSetting(MODULENAME "/circuit/Lp", "Plasma self-inductance", (real_t)0);
    s->DefineSetting(MODULENAME "/circuit/Lwall", "Wall self-inductance", (real_t)0);
    s->DefineSetting(MODULENAME "/circuit/M", "Plasma-wall mutual inductance", (real_t)0);
    s->DefineSetting(MODULENAME "/circuit/Rwall", "Wall resistance", (real_t)0);
    s->DefineSetting(MODULENAME "/circuit/Iwall0", "Wall current at t=0", (real_t)0);

    DREAM::SimulationGenerator::DefineDataT(MODULENAME "/circuit", s, "Vloop");
}

void SimulationGenerator::ConstructEquation_E_field(
    EquationSystem *eqsys, DREAM::Settings *s,  struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    enum OptionConstants::eqterm_E_field_eqn type =
        (enum OptionConstants::eqterm_E_field_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED:
            throw DREAM::SettingsException(
                "STREAM does not allow for the electric field evolution to be prescribed."
            );

        case OptionConstants::UQTY_E_FIELD_EQN_SELFCONSISTENT:
            ConstructEquation_E_field_selfconsistent(eqsys, s, oqty_terms);
            break;

        case OptionConstants::UQTY_E_FIELD_EQN_CIRCUIT:
            ConstructEquation_E_field_circuit(eqsys, s);
            break;

        default:
            throw DREAM::SettingsException(
                "Unrecognized equation type for '%s': %d.",
                DREAM::OptionConstants::UQTY_E_FIELD, type
            );
    }
}

/**
 * Construct the usual self-consistent electric field equation.
 */
void SimulationGenerator::ConstructEquation_E_field_selfconsistent(
    EquationSystem *eqsys, DREAM::Settings *s, struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms
) {
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    // Define poloidal fluxes
    eqsys->SetUnknown(
        DREAM::OptionConstants::UQTY_POL_FLUX,
        DREAM::OptionConstants::UQTY_POL_FLUX_DESC,
        fluidGrid
    );
    eqsys->SetUnknown(
        DREAM::OptionConstants::UQTY_PSI_EDGE,
        DREAM::OptionConstants::UQTY_PSI_EDGE_DESC,
        fluidGrid
    );
    eqsys->SetUnknown(
        DREAM::OptionConstants::UQTY_PSI_WALL,
        DREAM::OptionConstants::UQTY_PSI_WALL_DESC,
        fluidGrid
    );

    // Define equation for E
    DREAM::SimulationGenerator::ConstructEquation_E_field_selfconsistent(
        eqsys, s, oqty_terms
    );

    // Define equation for poloidal flux
    DREAM::SimulationGenerator::ConstructEquation_psi_p(eqsys, s);
    DREAM::SimulationGenerator::ConstructEquation_psi_edge(eqsys, s);

    ResetPoloidalFluxInitialization(eqsys, s);
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

/**
 * Construct the STREAM analogue of the DYON circuit equations.
 */
void SimulationGenerator::ConstructEquation_E_field_circuit(
    EquationSystem *eqsys, DREAM::Settings *s
) {
    DREAM::FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    DREAM::FVM::Grid *fluidGrid = eqsys->GetFluidGrid();


    // Define wall current as an unknown
    eqsys->SetUnknown(
        DREAM::OptionConstants::UQTY_I_WALL,
        DREAM::OptionConstants::UQTY_I_WALL_DESC,
        scalarGrid
    );
    // Define external loop voltage as an unknown
    eqsys->SetUnknown(
        OptionConstants::UQTY_VLOOP,
        OptionConstants::UQTY_VLOOP_DESC,
        scalarGrid
    );
    
    const len_t id_E_field = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
    const len_t id_V_loop  = eqsys->GetUnknownID(OptionConstants::UQTY_VLOOP);
    const len_t id_I_p     = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    const len_t id_I_w     = eqsys->GetUnknownID(DREAM::OptionConstants::UQTY_I_WALL);

    // Get circuit parameters
    const real_t Lp    = s->GetReal(MODULENAME "/circuit/Lp");
    const real_t Lwall = s->GetReal(MODULENAME "/circuit/Lwall");
    const real_t M     = s->GetReal(MODULENAME "/circuit/M");
    const real_t Rwall = s->GetReal(MODULENAME "/circuit/Rwall");

    const real_t R0    = eqsys->GetEllipticalRadialGridGenerator()->GetMajorRadius();

    if (Lp <= 0)
        throw DREAM::SettingsException(
            "%s: parameter '" MODULENAME "/circuit/Lp' has an invalid value. Must be > 0.",
            DREAM::OptionConstants::UQTY_E_FIELD
        );
    else if (Lwall <= 0)
        throw DREAM::SettingsException(
            "%s: parameter '" MODULENAME "/circuit/Lwall' has an invalid value. Must be > 0.",
            DREAM::OptionConstants::UQTY_E_FIELD
        );
    else if (M <= 0)
        throw DREAM::SettingsException(
            "%s: parameter '" MODULENAME "/circuit/M' has an invalid value. Must be > 0.",
            DREAM::OptionConstants::UQTY_E_FIELD
        );
    else if (Rwall <= 0)
        throw DREAM::SettingsException(
            "%s: parameter '" MODULENAME "/circuit/Rwall' has an invalid value. Must be > 0.",
            DREAM::OptionConstants::UQTY_E_FIELD
        );

    // Load loop voltage
    DREAM::FVM::Interpolator1D *Vloop =
        DREAM::SimulationGenerator::LoadDataT(MODULENAME "/circuit", s, "Vloop");

    ///////////////////
    // LOOP VOLTAGE
    ///////////////////
    DREAM::FVM::Operator *Op1 = new DREAM::FVM::Operator(scalarGrid);
    Op1->AddTerm(new DREAM::FVM::PrescribedParameter(scalarGrid, Vloop));
    eqsys->SetOperator(id_V_loop, id_V_loop, Op1, "Prescribed");

    ///////////////////
    // ELECTRIC FIELD
    ///////////////////
    DREAM::FVM::Operator *Op_E_Efield = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_E_Vloop  = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_E_Ip     = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_E_Iw     = new DREAM::FVM::Operator(fluidGrid);

    Op_E_Efield->AddTerm(new VloopTerm(fluidGrid));
    Op_E_Vloop->AddTerm(new DREAM::FVM::IdentityTerm(fluidGrid, -1/R0));
    Op_E_Ip->AddTerm(new DREAM::FVM::TransientTerm(fluidGrid, id_I_p, Lp/R0));
    Op_E_Iw->AddTerm(new DREAM::FVM::TransientTerm(fluidGrid, id_I_w, M/R0));

    eqsys->SetOperator(id_E_field, id_E_field, Op_E_Efield, "Vloop = 2*pi*R*E + Lp*dIp/dt + M*dIw/dt");
    eqsys->SetOperator(id_E_field, id_V_loop, Op_E_Vloop);
    eqsys->SetOperator(id_E_field, id_I_p, Op_E_Ip);
    eqsys->SetOperator(id_E_field, id_I_w, Op_E_Iw);

    ///////////////////
    // WALL CURRENT
    ///////////////////
    DREAM::FVM::Operator *Op_Iw_Vloop = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_Iw_Iw = new DREAM::FVM::Operator(fluidGrid);
    DREAM::FVM::Operator *Op_Iw_Ip = new DREAM::FVM::Operator(fluidGrid);

    Op_Iw_Vloop->AddTerm(new DREAM::FVM::IdentityTerm(fluidGrid, -1.0/R0));
    Op_Iw_Iw->AddTerm(new DREAM::FVM::IdentityTerm(fluidGrid, Rwall/R0));
    Op_Iw_Iw->AddTerm(new DREAM::FVM::TransientTerm(fluidGrid, id_I_w, Lwall/R0));
    Op_Iw_Ip->AddTerm(new DREAM::FVM::TransientTerm(fluidGrid, id_I_p, M/R0));

    eqsys->SetOperator(id_I_w, id_V_loop, Op_Iw_Vloop, "Vloop = Rw*Iw + Lw*dIw/dt + M*dIp/dt");
    eqsys->SetOperator(id_I_w, id_I_w, Op_Iw_Iw);
    eqsys->SetOperator(id_I_w, id_I_p, Op_Iw_Ip);

    // Set initial values
    real_t *Efield_init = DREAM::SimulationGenerator::LoadDataR(
        MODULENAME, fluidGrid->GetRadialGrid(), s, "init"
    );
    eqsys->SetInitialValue(id_E_field, Efield_init);
    delete [] Efield_init;

    // I_w(t=0) = 0
    const real_t Iwall0 = s->GetReal(MODULENAME "/circuit/Iwall0");
    eqsys->SetInitialValue(id_I_w, &Iwall0);

    // Vloop(t=0) = prescribed value
    eqsys->initializer->AddRule(
        id_V_loop,
        DREAM::EqsysInitializer::INITRULE_EVAL_EQUATION
    );
}

