#ifndef _STREAM_EQUATION_SYSTEM_HPP
#define _STREAM_EQUATION_SYSTEM_HPP

#include <vector>
#include "DREAM/EquationSystem.hpp"
#include "STREAM/Equations/ConnectionLength.hpp"
#include "STREAM/Equations/ConfinementTime.hpp"
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/PlasmaVolume.hpp" 
#include "STREAM/Equations/NeutralInflux.hpp"
#include "STREAM/Equations/IonRateEquation.hpp"
#include "STREAM/Equations/RunawayElectronConfinementTime.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "STREAM/Equations/OpticalThickness.hpp"

namespace STREAM { 
    class EquationSystem : public DREAM::EquationSystem {
        public: 
            ConnectionLength *CL=nullptr;
            ConfinementTime *CT=nullptr;
            NeutralInflux *NI=nullptr;
            PlasmaVolume *PV=nullptr;
            EllipticalRadialGridGenerator *r=nullptr;
            std::vector<IonRateEquation*> ire;
            RunawayElectronConfinementTime *rect=nullptr;
            OpticalThickness *OT=nullptr;
            
            EquationSystem(
                DREAM::FVM::Grid*, DREAM::FVM::Grid*,
                enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*,
                enum DREAM::OptionConstants::momentumgrid_type, DREAM::FVM::Grid*, 
                DREAM::Settings*, EllipticalRadialGridGenerator*
            );
            
            void SetConnectionLength(ConnectionLength *CL);
            void SetConfinementTime(ConfinementTime *CT);
            void SetPlasmaVolume(PlasmaVolume *PV);
            void SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r);
            void SetNeutralInflux(NeutralInflux *NI);
            void AddIonRateEquation(IonRateEquation*);
            void SetRunawayElectronConfinementTime(RunawayElectronConfinementTime*);
            void SetOpticalThickness(OpticalThickness *OT);
            
            PlasmaVolume *GetPlasmaVolume();
            ConnectionLength *GetConnectionLength();
            ConfinementTime *GetConfinementTime();
            EllipticalRadialGridGenerator *GetEllipticalRadialGridGenerator();
            NeutralInflux *GetNeutralInflux();
            std::vector<IonRateEquation*> GetIonRateEquations();
            RunawayElectronConfinementTime *GetRunawayElectronConfinementTime();
            OpticalThickness *GetOpticalThickness();
            
    };
}
#endif /*_STREAM_EQUATION_SYSTEM_HPP*/
