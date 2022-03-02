#ifndef _STREAM_EQUATION_SYSTEM_HPP
#define _STREAM_EQUATION_SYSTEM_HPP

#include <vector>
#include "DREAM/EquationSystem.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp" 
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Equations/IonRateEquation.hpp"
#include "STREAM/Equations/RunawayElectronConfinementTime.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"

namespace STREAM { 
    class EquationSystem : public DREAM::EquationSystem {
        public: 
            ConfinementTime *CT;
            NeutralInflux *NI;
            PlasmaVolume *PV;
            EllipticalRadialGridGenerator *r;
            std::vector<IonRateEquation*> ire;
            RunawayElectronConfinementTime *rect;
            
            EquationSystem(
                DREAM::FVM::Grid*, DREAM::FVM::Grid*,
                enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*,
                enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, 
                DREAM::Settings*, EllipticalRadialGridGenerator*
            );
            
            void SetConfinementTime(ConfinementTime *CT);
            void SetPlasmaVolume(PlasmaVolume *PV);
            void SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r);
            void SetNeutralInflux(NeutralInflux *NI);
            void AddIonRateEquation(IonRateEquation*);
            void SetRunawayElectronConfinementTime(RunawayElectronConfinementTime*);
            
            PlasmaVolume *GetPlasmaVolume();
            ConfinementTime *GetConfinementTime();
            EllipticalRadialGridGenerator *GetEllipticalRadialGridGenerator();
            NeutralInflux *GetNeutralInflux();
            std::vector<IonRateEquation*> GetIonRateEquations();
            RunawayElectronConfinementTime *GetRunawayElectronConfinementTime();
            
    };
}
#endif /*_STREAM_EQUATION_SYSTEM_HPP*/
