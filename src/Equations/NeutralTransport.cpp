#include "STREAM/Equations/NetralTransport.hpp"

using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */
NeutralTransport::NeutralTransport(NeutralInflux *NI, PlasmaVolume *PV, real_t vessel_vol) : NI(NI), PV(PV), vessel_vol(vessel_vol) {
    SetName("NeutralTransport");
}

real_t void GetNeutralTransport(real_t t, const len_t iIon){
    Gamma0 = this->NI->EvaluateNeutralInflux(t, iIon);
    V_p = PV->GetPlasmaVolume();
    V_ni= PV->GetNeutralVolume(iIon);
    
    return Gamma0/(vessel_vol-V_p+V_ni);
}
