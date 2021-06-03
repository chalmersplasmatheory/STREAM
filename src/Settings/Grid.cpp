/**
 * Constructor of the main grid in STREAM.
 */

#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


#define MODULENAME "radialgrid"

using namespace STREAM;

DREAM::FVM::Grid *SimulationGenerator::ConstructRadialGrid(DREAM::Settings *s) {
    // TODO add a radial grid with elongated flux surfaces and variable minor radius
    DREAM::FVM::RadialGrid *rg = ConstructRadialGrid_Cylindrical(s);

    return new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));
}

/**
 * Construction of basic cylindrical radial grid.
 */
DREAM::FVM::RadialGrid *SimulationGenerator::ConstructRadialGrid_Cylindrical(
    DREAM::Settings *s
) {
    real_t B0 = s->GetReal(MODULENAME "/B0");   // magnetic field strength
    real_t a  = s->GetReal(MODULENAME "/a");    // plasma minor radius (constant)

    const real_t r0 = 0;    // innermost radius simulated
    const len_t nr = 1;     // number of radial grid points

    return new DREAM::FVM::RadialGrid(
        new DREAM::FVM::CylindricalRadialGridGenerator(
            nr, B0, r0, a
        )
    );
}

