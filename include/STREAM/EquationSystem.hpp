#ifndef _STREAM_EQUATION_SYSTEM_HPP
#define _STREAM_EQUATION_SYSTEM_HPP

namespace STREAM {
    class EquationSystem : DREAM::EquationSystem {
        public: 
            ConfinementTime *CT;
            NeutralInflux *NI;
            
            EquationSystem(DREAM::FVM::Grid*, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, 
            ConfinementTime*, NeutralInflux*);
    }
}
