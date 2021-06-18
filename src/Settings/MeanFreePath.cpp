#include "STREAM/Settings/SimulationGenerator.hpp"
#include "STREAM/Equations/MeanFreePathTerm.hpp"
#include "STREAM/Settings/OptionConstants.hpp"

using namespace STREAM;

    void SimulationGenerator::ConstructEquation_lambda_i(
        DREAM::EquationSystem *eqsys, DREAM::Settings *s, DREAM::ADAS *adas){
        DREAM::FVM::Operator *op_lambda_i = new DREAM::FVM::Operator(eqsys->GetFluidGrid());
        DREAM::FVM::Operator *op_W_i = new DREAM::FVM::Operator(eqsys->GetFluidGrid());
        
        DREAM::IonHandler *ions = eqsys->GetIonHandler();
        
        for (len_t iz = 0; iz<ions->GetNZ(); iz++){
            op_lambda_i->AddTerm(new IonIdentityTerm(eqsys->GetFluidGrid(), iz));
            op_W_i->AddTerm(new MeanFreePathTerm(eqsys->GetFluidGrid(), iz, eqsys->GetUnknownHandler(), adas, ions));
        }
        
        eqsys->SetOperator(OptionConstants::UQTY_LAMBDA_I, OptionConstants::UQTY_LAMBDA_I, op_lambda_i, "lambda_i = v_i/(n_e * I_i^(0))");
        eqsys->SetOperator(OptionConstants::UQTY_LAMBDA_I, DREAM::OptionConstants::UQTY_WI_ENER, op_W_i);
        
        eqsys->initializer->AddRule(OptionConstants::UQTY_LAMBDA_I, DREAM::EqsysInitializer::INITRULE_EVAL_EQUATION);
    }
