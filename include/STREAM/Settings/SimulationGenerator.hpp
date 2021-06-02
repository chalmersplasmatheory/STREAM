#ifndef _STREAM_SIMULATION_GENERATOR_HPP
#define _STREAM_SIMULATION_GENERATOR_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/AMJUEL.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

namespace STREAM {
    class SimulationGenerator {
    public:
        static DREAM::Simulation *ProcessSettings(DREAM::Settings*);

        // Equation system
        static DREAM::EquationSystem *ConstructEquationSystem(
            DREAM::Settings*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*,
            enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*
        );
        static void ConstructEquations(
            DREAM::EquationSystem*, DREAM::Settings*, DREAM::ADAS*,
            DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        static void ConstructUnknowns(
            DREAM::EquationSystem*, DREAM::Settings*, DREAM::ADAS*,
            DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );

        // Ion density equation
        static void ConstructEquation_Ions(
            DREAM::EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*
        );

        // Temperature equation
        static void ConstructEquation_T_cold(
            DREAM::EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::NIST*,
            struct OtherQuantityHandler::eqn_terms*
        );
        static void ConstructEquation_T_cold_selfconsistent(
            DREAM::EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::NIST*,
            struct OtherQuantityHandler::eqn_terms*
        );

        // STREAM main grid
        static DREAM::FVM::Grid *ConstructRadialGrid();
        static DREAM::FVM::RadialGrid *ConstructRadialGrid_Cylindrical(Settings*);
    };
}

#endif/*_STREAM_SIMULATION_GENERATOR_HPP*/
