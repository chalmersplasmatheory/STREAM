#include "STREAM/EquationSystem.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
 EquationSystem::EquationSystem(
    FVM::Grid *emptygrid, FVM::Grid *rgrid,
    enum OptionConstants::momentumgrid_type ht_type, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid,
    ConfinementTime *CT, NeutralInflux *NI
) : EquationSystem(emptygrid, rgrid, hottailGrid, runawayGrid, ht_type, re_type),
    CT(CT), NI(NI) {}
    
void EquationSystem::SetConfinementTime(ConfinementTime *CT)
    { this->CT = CT; }
    
void EquationSystem::SetNeutralInflux(NeutralInflux *NI)
    { this->NI = NI; }
    
EllipticalRadialGridGenerator *GetEllipticalRadialGridGenerator()
    {return this->r}
SputteredRecycledCoefficient *GetSputteredRecycledCoefficient()
    {return this->SRC}
PlasmaVolume *GetPlasmaVolume()
    {return this->PV}
ConfinementTime *GetConfinementTime()
    {return this->CT}

void *SetEllipticalRadialGridGenerator(EllipticalRadialGridGenerator* r){
    this->r=r;
}
void *SetSputteredRecycledCoefficient(SputteredRecycledCoefficient *SRC){
    this->SRC=SRC;
}
void *SetPlasmaVolume(PlasmaVolume *PV){
    this->PV=PV;
}
    
    
    

