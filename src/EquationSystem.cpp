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
    ConfinementTime *CT, NeutralInflux *NI
) : DREAM::EquationSystem(emptygrid, rgrid, ht_type, hottailGrid, re_type, runawayGrid),
    CT(CT), NI(NI) {}
    
void EquationSystem::SetConfinementTime(ConfinementTime *CT)
    { this->CT = CT; }
    
void EquationSystem::SetNeutralInflux(NeutralInflux *NI)
    { this->NI = NI; }
    
EllipticalRadialGridGenerator* EquationSystem::GetEllipticalRadialGridGenerator()
    {return this->r;}
SputteredRecycledCoefficient* EquationSystem::GetSputteredRecycledCoefficient()
    {return this->SRC;}
PlasmaVolume* EquationSystem::GetPlasmaVolume()
    {return this->PV;}
ConfinementTime* EquationSystem::GetConfinementTime()
    {return this->CT;}

void EquationSystem::SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator *r){
    this->r=r;
}
void EquationSystem::SetSputteredRecycledCoefficient(SputteredRecycledCoefficient *SRC){
    this->SRC=SRC;
}
void EquationSystem::SetPlasmaVolume(PlasmaVolume *PV){
    this->PV=PV;
}
    
    
    

