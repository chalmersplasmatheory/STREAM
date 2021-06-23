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
        public: 
            EllipticalRadialGridGenerator *r=nullptr;
            //PlasmaVolume *PV=nullptr;  Antar ska vara som ConfinementTime istället för EllipticalRadialGrid?

        public: 
            ConfinementTime *CT;
            PlasmaVolume *PV
            
            EquationSystem(DREAM::FVM::Grid*, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, 
            ConfinementTime*, PlasmaVolume*);
            
            void SetConfinementTime(ConfinementTime *CT);
            //void SetNeutralInflux(NeutralInflux *NI);
            void SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r);
            //void SetSputteredRecycledCoefficient(SputteredRecycledCoefficient *SRC);
            void SetPlasmaVolume(PlasmaVolume *PV);
            
            EllipticalRadialGridGenerator *GetEllipticalRadialGridGenerator();
            //SputteredRecycledCoefficient *GetSputteredRecycledCoefficient();
            PlasmaVolume *GetPlasmaVolume();
            ConfinementTime *GetConfinementTime();
            //NeutralInflux *GetNeutralInflux();
            
    };
}
#endif /*_STREAM_EQUATION_SYSTEM_HPP*/
