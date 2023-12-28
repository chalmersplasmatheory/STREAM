#include "STREAM/EquationSystem.hpp"

using namespace STREAM;
using namespace std;

/**
 * Constructor. 
 */
EquationSystem::EquationSystem(
    DREAM::FVM::Grid *emptygrid, DREAM::FVM::Grid *rgrid,
    enum DREAM::OptionConstants::momentumgrid_type ht_type, DREAM::FVM::Grid *hottailGrid,
    enum DREAM::OptionConstants::momentumgrid_type re_type, DREAM::FVM::Grid *runawayGrid,
    DREAM::Settings *s, EllipticalRadialGridGenerator *r
) : DREAM::EquationSystem(emptygrid, rgrid, ht_type, hottailGrid, re_type, runawayGrid, s), r(r)
    {/*s->DisplaySettings();*/}
    
EllipticalRadialGridGenerator* EquationSystem::GetEllipticalRadialGridGenerator()
    {return this->r;}
PlasmaVolume* EquationSystem::GetPlasmaVolume()
    {return this->PV;}
ConnectionLength* EquationSystem::GetConnectionLength()
    {return this->CL;}
ConfinementTime* EquationSystem::GetConfinementTime()
    {return this->CT;}
NeutralInflux *EquationSystem::GetNeutralInflux()
    {return this->NI;}
vector<IonRateEquation*> EquationSystem::GetIonRateEquations()
    { return ire; }
RunawayElectronConfinementTime *EquationSystem::GetRunawayElectronConfinementTime()
    { return rect; }
OpticalThickness* EquationSystem::GetOpticalThickness()
    {return this->OT;}
    
void EquationSystem::SetConnectionLength(ConnectionLength *CL)
    { this->CL = CL; }
    
void EquationSystem::SetConfinementTime(ConfinementTime *CT)
    { this->CT = CT; }

void EquationSystem::SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r){
    this->r=r;
}

void EquationSystem::SetNeutralInflux(NeutralInflux *NI) {
    this->NI = NI;
}

void EquationSystem::SetPlasmaVolume(PlasmaVolume *PV){
    this->PV=PV;
}
    
void EquationSystem::AddIonRateEquation(IonRateEquation *IRE) {
    this->ire.push_back(IRE);
}

void EquationSystem::SetRunawayElectronConfinementTime(RunawayElectronConfinementTime *r) {
    this->rect = r;
}

void EquationSystem::SetOpticalThickness(OpticalThickness *OT)
    { this->OT = OT; }
