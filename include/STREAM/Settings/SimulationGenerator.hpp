#ifndef _STREAM_SIMULATION_GENERATOR_HPP
#define _STREAM_SIMULATION_GENERATOR_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/AMJUEL.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/EquationSystem.hpp"

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
        static void DefineOptions_wall(DREAM::Settings *s);

        // Equation system
        static EquationSystem *ConstructEquationSystem(
            DREAM::Settings*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*
        );
        static void ConstructEquations(
            EquationSystem*, DREAM::Settings*, DREAM::ADAS*,
            DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        static void ConstructUnknowns(
            EquationSystem*, DREAM::Settings*
        );

        // Ion density equation
        static void ConstructEquation_Ions(
            EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*
        );

        // Temperature equation
        static void ConstructEquation_T_cold(
            EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        static void ConstructEquation_T_cold_selfconsistent(
            EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        
        //Mean free path equation
        static void ConstructEquation_lambda_i(
            EquationSystem*, DREAM::Settings*, DREAM::ADAS*
        );
        
        //Ion heat transport
        static void ConstructEquation_ion_heat_transport(
            EquationSystem*, DREAM::Settings*
        );
        
        //Ion heat transport
        static void ConstructEquation_ion_transport(
            EquationSystem*, DREAM::Settings*
        );
        
        // General transport interface
        static bool ConstructTransportTerm(
            DREAM::FVM::Operator*, const std::string&, DREAM::FVM::Grid*,
            enum DREAM::OptionConstants::momentumgrid_type, EquationSystem*,
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
