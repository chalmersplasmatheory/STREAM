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
            ConfinementTime *CT;
            PlasmaVolume *PV;
            EllipticalRadialGridGenerator *r;
            
            EquationSystem(DREAM::FVM::Grid*, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, 
            ConfinementTime*, PlasmaVolume*, EllipticalRadialGridGenerator*);
            
            void SetConfinementTime(ConfinementTime *CT);
            void SetPlasmaVolume(PlasmaVolume *PV);
            void SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r);
            
            PlasmaVolume *GetPlasmaVolume();
            ConfinementTime *GetConfinementTime();
            EllipticalRadialGridGenerator *GetEllipticalRadialGridGenerator();
            
    };
}
#endif /*_STREAM_EQUATION_SYSTEM_HPP*/
