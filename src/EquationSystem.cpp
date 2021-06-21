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
