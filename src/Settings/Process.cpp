/**
 * Process settings and build a Simulation object.
 */

#include "DREAM/ADAS.hpp"
#include "DREAM/AMJUEL.hpp"
#include "DREAM/NIST.hpp"
#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"

using namespace STREAM;

/**
 * Construct a new settings object with all available options
 * defined.
 */
DREAM::Settings *SimulationGenerator::CreateSettings() {
    DREAM::Settings *s = new DREAM::Settings();
    DefineOptions(s);

    return s;
}

/**
 * Process the given settings and construct a simulation object.
 *
 * s: Settings specifying how to construct the simulation.
 */
DREAM::Simulation *SimulationGenerator::ProcessSettings(DREAM::Settings *s) {
    const real_t t0 = 0;

    // Construct grids
    enum DREAM::OptionConstants::momentumgrid_type ht_type, re_type;
    DREAM::FVM::Grid *scalarGrid  = DREAM::SimulationGenerator::ConstructScalarGrid();

    // Fluid grid
    EllipticalRadialGridGenerator *ergg = ConstructRadialGrid_Elliptical(s);
    auto rg = new DREAM::FVM::RadialGrid(ergg);
    DREAM::FVM::Grid *fluidGrid   = new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));

    DREAM::FVM::Grid *hottailGrid = DREAM::SimulationGenerator::ConstructHotTailGrid(s, fluidGrid->GetRadialGrid(), &ht_type);

    scalarGrid->Rebuild(t0);
    fluidGrid->Rebuild(t0);
    if (hottailGrid)
        hottailGrid->Rebuild(t0);

    // The runaway grid depends on the hot-tail grid (if it exists)
    DREAM::FVM::Grid *runawayGrid = DREAM::SimulationGenerator::ConstructRunawayGrid(s, fluidGrid->GetRadialGrid(), hottailGrid, &re_type);
    if (runawayGrid)
        runawayGrid->Rebuild(t0);

    // Load databases
    DREAM::ADAS *adas     = DREAM::SimulationGenerator::LoadADAS(s);
    DREAM::AMJUEL *amjuel = DREAM::SimulationGenerator::LoadAMJUEL(s);
    DREAM::NIST *nist     = DREAM::SimulationGenerator::LoadNIST(s);

    // Construct equation system
    EquationSystem *eqsys = ConstructEquationSystem(
        s, scalarGrid, fluidGrid, hottailGrid, runawayGrid, adas, amjuel, nist, ergg
    );

    // Construct simulation object
    DREAM::Simulation *sim = new DREAM::Simulation();
    sim->SetADAS(adas);
    sim->SetNIST(nist);
    sim->SetAMJUEL(amjuel);
    sim->SetEquationSystem(eqsys);

    // Load data from specified output file
    DREAM::SimulationGenerator::LoadOutput(s, sim);

    return sim;
}

