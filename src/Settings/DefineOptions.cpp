/**
 * This file contains the routine 'DefineOptions', which
 * defines all available options in STREAM.
 */

#include "STREAM/Settings/SimulationGenerator.hpp"


using namespace STREAM;


/**
 * Define options for the given settings object.
 *
 * s: Settings object to define available options for.
 */
void SimulationGenerator::DefineOptions(DREAM::Settings *s) {
    // STREAM-specific
    DefineOptions_ElectricField(s);
    DefineOptions_Grid(s);
    DefineOptions_T_cold(s);
    DefineOptions_Ions(s);

    // Same as for DREAM
    DREAM::SimulationGenerator::DefineOptions_ADAS(s);
    DREAM::SimulationGenerator::DefineOptions_CollisionQuantityHandler(s);
    DREAM::SimulationGenerator::DefineOptions_EquationSystem(s);
    DREAM::SimulationGenerator::DefineOptions_Initializer(s);
    DREAM::SimulationGenerator::DefineOptions_HotTailGrid(s);
    DREAM::SimulationGenerator::DefineOptions_f_hot(s);
    DREAM::SimulationGenerator::DefineOptions_f_re(s);
    DREAM::SimulationGenerator::DefineOptions_j_ohm(s);
    DREAM::SimulationGenerator::DefineOptions_j_tot(s);
    DREAM::SimulationGenerator::DefineOptions_n_re(s);
    DREAM::SimulationGenerator::DefineOptions_Output(s);
    DREAM::SimulationGenerator::DefineOptions_RunawayGrid(s);
    DREAM::SimulationGenerator::DefineOptions_Solver(s);
    DREAM::SimulationGenerator::DefineOptions_TimeStepper(s);
    DREAM::SimulationGenerator::DefineOptions_OtherQuantities(s);
    DREAM::SimulationGenerator::DefineOptions_SPI(s);
}

