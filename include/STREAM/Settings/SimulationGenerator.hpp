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
        static DREAM::Settings *CreateSettings();
        static void DefineOptions(DREAM::Settings*);
        static DREAM::Simulation *ProcessSettings(DREAM::Settings*);

        // Define options
        static void DefineOptions_Grid(DREAM::Settings*);
        static void DefineOptions_T_cold(DREAM::Settings*);
        static void DefineOptions_Ions(DREAM::Settings*);
        static void DefineOptions_Transport(const std::string&, DREAM::Settings*, bool, const std::string& subname="transport");

        // Equation system
        static DREAM::EquationSystem *ConstructEquationSystem(
            DREAM::Settings*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*
        );
        static void ConstructEquations(
            DREAM::EquationSystem*, DREAM::Settings*, DREAM::ADAS*,
            DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        static void ConstructUnknowns(
            DREAM::EquationSystem*, DREAM::Settings*
        );

        // Ion density equation
        static void ConstructEquation_Ions(
            DREAM::EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*
        );

        // Temperature equation
        static void ConstructEquation_T_cold(
            DREAM::EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        static void ConstructEquation_T_cold_selfconsistent(
            DREAM::EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        
        //Mean free path equation
        static void ConstructEquation_lambda_i(
            DREAM::EquationSystem*, DREAM::Settings*, DREAM::ADAS*
        );
        
        // General transport interface
        static bool ConstructTransportTerm(
            DREAM::FVM::Operator*, const std::string&, DREAM::FVM::Grid*,
            enum DREAM::OptionConstants::momentumgrid_type, DREAM::EquationSystem*,
            DREAM::Settings*, bool, bool, DREAM::TransportAdvectiveBC** abc=nullptr,
            DREAM::TransportDiffusiveBC** dbc=nullptr,
            struct DREAM::OtherQuantityHandler::eqn_terms *oqty_terms=nullptr,
            const std::string& subname="transport"
        );
        // STREAM main grid
        static DREAM::FVM::Grid *ConstructRadialGrid(DREAM::Settings*);
        static DREAM::FVM::RadialGrid *ConstructRadialGrid_Cylindrical(DREAM::Settings*);
        static DREAM::FVM::RadialGrid *ConstructRadialGrid_Elliptical(DREAM::Settings*);
    };
}

#endif/*_STREAM_SIMULATION_GENERATOR_HPP*/
