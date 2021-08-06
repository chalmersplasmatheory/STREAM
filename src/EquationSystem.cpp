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
    ConfinementTime *CT, PlasmaVolume *PV, EllipticalRadialGridGenerator *r
) : DREAM::EquationSystem(emptygrid, rgrid, ht_type, hottailGrid, re_type, runawayGrid),
    CT(CT), PV(PV), r(r) {}
    
EllipticalRadialGridGenerator* EquationSystem::GetEllipticalRadialGridGenerator()
    {return this->r;}
PlasmaVolume* EquationSystem::GetPlasmaVolume()
    {return this->PV;}
ConfinementTime* EquationSystem::GetConfinementTime()
    {return this->CT;}
NeutralInflux *EquationSystem::GetNeutralInflux()
    {return this->NI;}
    
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
    
    
    

