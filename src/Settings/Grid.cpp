/**
 * Constructor of the main grid in STREAM.
 */

#include "FVM/Grid/CylindricalRadialGridGenerator.hpp"
#include "FVM/Grid/EmptyMomentumGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Settings/SimulationGenerator.hpp"


#define MODULENAME "radialgrid"

using namespace STREAM;


void SimulationGenerator::DefineOptions_Grid(DREAM::Settings *s) {
    s->DefineSetting(MODULENAME "/wall_radius",  "Tokamak wall radius", (real_t)0.5);
    s->DefineSetting(MODULENAME "/R0",  "Tokamak major radius", (real_t)0);
    
    DREAM::SimulationGenerator::DefineDataT(MODULENAME, s, "a");
    DREAM::SimulationGenerator::DefineDataT(MODULENAME, s, "B0");
    DREAM::SimulationGenerator::DefineDataT(MODULENAME, s, "kappa");
    DREAM::SimulationGenerator::DefineDataT(MODULENAME, s, "delta");
    
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

    s->DefineSetting(
        "radialgrid/Iref",
        "Reference plasma current at which flux surfaces form [A]",
        (real_t)100e3
    );

    s->DefineSetting(
        "radialgrid/Bv",
        "Stray vertical magnetic field [T]",
        (real_t)1e-3
    );
    
    s->DefineSetting(
        "radialgrid/uncertaintyFactor",
        "Runaway confinement uncertainty factor",
        (real_t)1.0
    );
    
    DREAM::SimulationGenerator::DefineOptions_f_ripple(MODULENAME, s);
}

/**
 * Main routine for creating a new main grid in STREAM.
 *
 * s: Settings object containing a specification of the grid.
 */
/*DREAM::FVM::Grid *SimulationGenerator::ConstructRadialGrid(DREAM::Settings *s) {
    //DREAM::FVM::RadialGrid *rg = ConstructRadialGrid_Cylindrical(s);
    DREAM::FVM::RadialGrid *rg = ConstructRadialGrid_Elliptical(s);

    return new DREAM::FVM::Grid(rg, new DREAM::FVM::EmptyMomentumGrid(rg));
}*/

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

/**
 * Construction of basic elliptical radial grid.
 */
EllipticalRadialGridGenerator *SimulationGenerator::ConstructRadialGrid_Elliptical(
    DREAM::Settings *s
) {
    // Plasma minor radius
    DREAM::FVM::Interpolator1D *a
        = DREAM::SimulationGenerator::LoadDataT(MODULENAME, s, "a");
    // Magnetic field strength on-axis
    DREAM::FVM::Interpolator1D *B0
        = DREAM::SimulationGenerator::LoadDataT(MODULENAME, s, "B0");
    // Plasma elongation
    DREAM::FVM::Interpolator1D *kappa
        = DREAM::SimulationGenerator::LoadDataT(MODULENAME, s, "kappa");
    // Plasma triangularity
    DREAM::FVM::Interpolator1D *delta
        = DREAM::SimulationGenerator::LoadDataT(MODULENAME, s, "delta");
    
    real_t R0 = s->GetReal(MODULENAME "/R0");
    
    if(R0 == 0)
        throw DREAM::SettingsException(
            "Tokamak major radius has not been specified."
        );

    return new EllipticalRadialGridGenerator(
        a, B0, kappa, delta, R0
    );
}

