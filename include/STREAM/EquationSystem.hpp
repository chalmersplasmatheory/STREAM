#ifndef _STREAM_EQUATION_SYSTEM_HPP
#define _STREAM_EQUATION_SYSTEM_HPP

#include "DREAM/EquationSystem.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp" 
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"

namespace STREAM {
    class EquationSystem : public DREAM::EquationSystem {
        private: 
            EllipticalRadialGridGenerator *r=nullptr;
            SputteredRecycledCoefficient *SRC=nullptr;
            PlasmaVolume *PV=nullptr;

        public: 
            ConfinementTime *CT;
            NeutralInflux *NI; 
            
            EquationSystem(DREAM::FVM::Grid*, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, 
            ConfinementTime*, NeutralInflux*);
            
            void SetConfinementTime(ConfinementTime *CT);
            void SetNeutralInflux(NeutralInflux *NI);
            
            EllipticalRadialGridGenerator *GetEllipticalRadialGridGenerator();
            SputteredRecycledCoefficient *GetSputteredRecycledCoefficient();
            PlasmaVolume *GetPlasmaVolume();
            ConfinementTime *GetConfinementTime();
            
            void SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r);
            void SetSputteredRecycledCoefficient(SputteredRecycledCoefficient *SRC);
            void SetPlasmaVolume(PlasmaVolume *PV);
    };
}
#endif /*_STREAM_EQUATION_SYSTEM_HPP*/
