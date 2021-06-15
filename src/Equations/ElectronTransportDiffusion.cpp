/**
 * Implementation of a transport operator taking the form
 *
 *   
 *
 * This operator should be applied to 'T_cold'.
 */

#include "DREAM/Constants.hpp"
#include "STREAM/Equations/ElectronTransportDiffusion.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;
using namespace STREAM;

/**
 * Constructor.
 */

ElectronTransportDiffusion::ElectronTransportDiffusion(
    FVM::Grid *grid, enum OptionConstants::momentumgrid_type mgtype,
    EllipticalRadialGridGenerator *radials, ConfinementTime *tauinv, FVM::UnknownQuantityHandler *unknowns
) : FVM::DiffusionTerm(grid), mgtype(mgtype), coefftauinv(tauinv) {

    SetName("ElectronTransportDiffusion"); 

    this->unknowns = unknowns;
    this->id_Ip    = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    this->id_Iwall = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi    = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    
    this->radials = radials;
    
    AddUnknownForJacobian(unknowns, this->id_Ip);
    AddUnknownForJacobian(unknowns, this->id_Iwall);
    AddUnknownForJacobian(unknowns, this->id_Tcold);
    AddUnknownForJacobian(unknowns, this->id_Wi); 
    AddUnknownForJacobian(unknowns, this->id_ni);

    AllocateDiffCoeff(); 
}

/**
 * Destructor.
 */
ElectronTransportDiffusion::~ElectronTransportDiffusion() {
    delete this->coefftauinv;
    delete [] this->dI_p;
    delete [] this->dI_wall;
    delete [] this->dT_cold;
    delete [] this->dW_i;
    delete [] this->dn_i;
}


/**
 * Allocate memory for the differentiation coefficient.
 */
void ElectronTransportDiffusion::AllocateDiffCoeff() {
    const len_t nr = this->grid->GetNr();
    this->dI_p = new real_t[nr+1];
    this->dI_wall = new real_t[nr+1];
    this->dT_cold = new real_t[nr+1];
    this->dW_i = new real_t[nr+1];
    this->dn_i = new real_t[nr+1];
}

/**
 * Called whenever the grid is rebuilt.
 */
bool ElectronTransportDiffusion::GridRebuilt() {
    this->FVM::DiffusionTerm::GridRebuilt();

    delete [] this->dI_p;
    delete [] this->dI_wall;
    delete [] this->dT_cold;
    delete [] this->dW_i;
    delete [] this->dn_i;
    AllocateDiffCoeff();
    
    return true;
}

/**
 * Rebuild the coefficients for this equation term.
 */
 
void ElectronTransportDiffusion::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns 
) {
    real_t a = radials->GetMinorRadius();
    
    const len_t nr = this->grid->GetNr();
    
    for (len_t ir = 0; ir < nr+1; ir++) {
        real_t tauinv        = this->coefftauinv->EvaluateConfinementTime(ir, t); 
        real_t dtauinvdIp    = this->coefftauinv->EvaluateConfinementTime_dIp(ir, t); 
        real_t dtauinvdIwall = this->coefftauinv->EvaluateConfinementTime_dIwall(ir, t); 
        real_t dtauinvdTcold = this->coefftauinv->EvaluateConfinementTime_dTe(ir, t); 
        real_t dtauinvdWi    = this->coefftauinv->EvaluateConfinementTime_dWi(ir, t); 
        real_t dtauinvdni    = this->coefftauinv->EvaluateConfinementTime_dni(ir, t);
         


        this->dI_p[ir]    = a * a * dtauinvdIp;
        this->dI_wall[ir] = a * a * dtauinvdIwall;
        this->dT_cold[ir] = a * a * dtauinvdTcold;
        this->dW_i[ir]    = a * a * dtauinvdWi;
        this->dn_i[ir]    = a * a * dtauinvdni;
        
        Drr(ir, 0, 0) += a * a * tauinv;
    }
}

/**
 * Set jacobian of diffusion coefficients for this diffusion term.
 *
 * derivId:    ID of the quantity with respect to which the coefficient should
 *             be differentiated.
 * nMultiples: (not used).
 */
void ElectronTransportDiffusion::SetPartialDiffusionTerm(
    len_t derivId, len_t
) {
    if (derivId != this->id_Ip && derivId != this->id_Iwall && derivId != this->id_Tcold && derivId != this->id_Wi && derivId != this->id_ni)
        return;

    ResetDifferentiationCoefficients();
    
    const len_t nr = this->grid->GetNr();
    
    if (derivId == id_Ip) {
        for (len_t ir = 0; ir < nr+1; ir++)
            dDrr(ir, 0, 0) = this->dI_p[ir];
    }
    else if (derivId == id_Iwall) {
        for (len_t ir = 0; ir < nr+1; ir++)
            dDrr(ir, 0, 0) = this->dI_wall[ir];
    }
    else if (derivId == id_Tcold) {
        for (len_t ir = 0; ir < nr+1; ir++)
            dDrr(ir, 0, 0) = this->dT_cold[ir];
    }
    else if (derivId == id_Wi) {
        for (len_t ir = 0; ir < nr+1; ir++)
            dDrr(ir, 0, 0) = this->dW_i[ir];
    }
    else if (derivId == id_ni) {
        for (len_t ir = 0; ir < nr+1; ir++)
            dDrr(ir, 0, 0) = this->dn_i[ir];
    }
}

