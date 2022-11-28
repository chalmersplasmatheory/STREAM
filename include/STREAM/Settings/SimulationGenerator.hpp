#ifndef _STREAM_SIMULATION_GENERATOR_HPP
#define _STREAM_SIMULATION_GENERATOR_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/AMJUEL.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/EquationSystem.hpp"
#include "STREAM/OtherQuantityHandler.hpp"

namespace STREAM { 
    class SimulationGenerator {
    public:
        static DREAM::Settings *CreateSettings();
        static void DefineOptions(DREAM::Settings*);
        static DREAM::Simulation *ProcessSettings(DREAM::Settings*);

        // Define options
        static void DefineOptions_ElectricField(DREAM::Settings*);
        static void DefineOptions_Grid(DREAM::Settings*);
        static void DefineOptions_T_cold(DREAM::Settings*);
        static void DefineOptions_Ions(DREAM::Settings*);
        static void DefineOptions_Transport(const std::string&, DREAM::Settings*, bool, const std::string& subname="transport");

        // Equation system
        static EquationSystem *ConstructEquationSystem(
            DREAM::Settings*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
			DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            EllipticalRadialGridGenerator*
        );
        static void ConstructEquations(
            EquationSystem*, DREAM::Settings*, DREAM::ADAS*,
            DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*,
            struct OtherQuantityHandler::eqn_terms*
        );
        static void ConstructOtherQuantityHandler(
            EquationSystem*, DREAM::Settings*,
            struct OtherQuantityHandler::eqn_terms*,
            struct DREAM::OtherQuantityHandler::eqn_terms*
        );
        static void ConstructUnknowns(
            EquationSystem*, DREAM::Settings*
        );

	// Distribution function
	static void ConstructEquation_f_hot(EquationSystem*, DREAM::Settings*, struct DREAM::OtherQuantityHandler::eqn_terms*, struct STREAM::OtherQuantityHandler::eqn_terms*);
	static void ConstructEquation_f_re(EquationSystem*, DREAM::Settings*, struct DREAM::OtherQuantityHandler::eqn_terms*, struct STREAM::OtherQuantityHandler::eqn_terms*);
	
        // Electric field and related equations
        static void ConstructEquation_E_field(EquationSystem*, DREAM::Settings*, struct DREAM::OtherQuantityHandler::eqn_terms*);
        static void ConstructEquation_E_field_selfconsistent(EquationSystem*, DREAM::Settings*, struct DREAM::OtherQuantityHandler::eqn_terms*);
        static void ConstructEquation_E_field_circuit(EquationSystem*, DREAM::Settings*);

        // Ion density equation
        static void ConstructEquation_Ions(
            EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*,
            struct OtherQuantityHandler::eqn_terms*
        );

        // Temperature equation
        static void ConstructEquation_T_cold(
            EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*,
            struct OtherQuantityHandler::eqn_terms*
        );
        static void ConstructEquation_T_cold_selfconsistent(
            EquationSystem*, DREAM::Settings*,
            DREAM::ADAS*, DREAM::AMJUEL*, DREAM::NIST*,
            struct DREAM::OtherQuantityHandler::eqn_terms*,
            struct OtherQuantityHandler::eqn_terms*
        );
        static void ConstructEquation_T_i(
            EquationSystem *eqsys, DREAM::Settings*, DREAM::ADAS*,
            struct OtherQuantityHandler::eqn_terms*
        );
        static void ConstructEquation_T_i_selfconsistent(
            EquationSystem *eqsys, DREAM::Settings* /*s*/, DREAM::ADAS *,
            struct OtherQuantityHandler::eqn_terms*
        );
        
        //Mean free path equation
        static void ConstructEquation_lambda_i(
            EquationSystem*, DREAM::Settings*, DREAM::ADAS*
        );

        static void ResetPoloidalFluxInitialization(EquationSystem*, DREAM::Settings*);
        
        
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
        static DREAM::FVM::RadialGrid *ConstructRadialGrid_Cylindrical(DREAM::Settings*);
        static EllipticalRadialGridGenerator *ConstructRadialGrid_Elliptical(DREAM::Settings*);
        
        static void DefineOptions_n_re(DREAM::Settings*);
        static void ConstructEquation_n_re(EquationSystem*, DREAM::Settings*, struct DREAM::OtherQuantityHandler::eqn_terms*);
        
        // Data loading routines
        static void DefineDataT_3D(const std::string&, DREAM::Settings*, const std::string& name="data");
        static DREAM::FVM::Interpolator1D ***LoadDataT_3D(const std::string&, DREAM::Settings*, const std::string& name="data");
    };
}

#endif/*_STREAM_SIMULATION_GENERATOR_HPP*/
