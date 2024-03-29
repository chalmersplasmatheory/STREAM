/**
 * Implementation of the ion rate equation
 *
 *   d n_i^(j) / dt =
 *      [I_i^(j-1) n_i^(j-1) n_cold Vhat_i^(j-1) - I_i^(j) n_i^(j) n_cold Vhat_i^(j) +
 *      R_i^(j+1) n_i^(j+1) n_cold Vhat_i^(j+1) - R_i^(j) n_i^(j) n_cold Vhat_i^(j)]/V_i^(j)
 *
 * where
 *
 *   I_i^(j)  = ionization rate coefficient for charge state 'j'
 *              of ion species 'i'.
 *   R_i^(j)  = radiative recombination rate for charge state 'j'
 *              of ion species 'i'.
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 * If using collfreq_mode FULL and kinetic ionization is used, 
 * ionization rates will here be set to 0 to avoid double counting.
 */

#include "DREAM/ADAS.hpp"
#include "DREAM/Constants.hpp"
#include "STREAM/Equations/IonRateEquation.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"
#include "STREAM/Settings/OptionConstants.hpp"

using namespace STREAM;
using namespace DREAM;


/**
 * Constructor.
 */
IonRateEquation::IonRateEquation(
    FVM::Grid *g, IonHandler *ihdl, const len_t iIon,
    ADAS *adas, FVM::UnknownQuantityHandler *unknowns, PlasmaVolume *volumes,
    bool addFluidIonization, bool addFluidJacobian, bool isAbl = false  
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), adas(adas), volumes(volumes),
    addFluidIonization(addFluidIonization), addFluidJacobian(addFluidJacobian) {
    
    SetName("IonRateEquation");

    this->unknowns  = unknowns;
    if(isAbl){
		this->id_ions   = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES_ABL);
    }else{
		this->id_ions   = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_ION_SPECIES);
		this->id_n_tot  = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_TOT);
    }

    this->id_n_cold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    this->id_T_cold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    this->id_lambda_i = unknowns->GetUnknownID(OptionConstants::UQTY_LAMBDA_I);
    this->id_Wi      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    this->id_Ni      = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);

    AllocateRateCoefficients();
}

/**
 * Destructor.
 */
IonRateEquation::~IonRateEquation() {
    DeallocateRateCoefficients();
}


/**
 * Allocate memory for storing the ionization rate coefficients.
 */
void IonRateEquation::AllocateRateCoefficients() {
    const len_t Nr  = this->grid->GetNr();

    this->Rec = new real_t*[(Zion+1)];
    this->PartialNRec = new real_t*[(Zion+1)];
    this->PartialTRec = new real_t*[(Zion+1)];
    this->Ion = new real_t*[(Zion+1)];
    this->PartialNIon = new real_t*[(Zion+1)];
    this->PartialTIon = new real_t*[(Zion+1)];

    // Diagnostic utilities
    this->posIonizTerm = new real_t*[(Zion+1)];
    this->negIonizTerm = new real_t*[(Zion+1)];
    this->posRecTerm = new real_t*[(Zion+1)];
    this->negRecTerm = new real_t*[(Zion+1)];
    this->posCXTerm = new real_t*[(Zion+1)];
    this->negCXTerm = new real_t*[(Zion+1)];

    this->Rec[0]         = new real_t[Nr*(Zion+1)];
    this->PartialNRec[0] = new real_t[Nr*(Zion+1)];
    this->PartialTRec[0] = new real_t[Nr*(Zion+1)];
    this->Ion[0]         = new real_t[Nr*(Zion+1)];
    this->PartialNIon[0] = new real_t[Nr*(Zion+1)];
    this->PartialTIon[0] = new real_t[Nr*(Zion+1)];

    this->posIonizTerm[0] = new real_t[Nr*(Zion+1)];
    this->negIonizTerm[0] = new real_t[Nr*(Zion+1)];
    this->posRecTerm[0] = new real_t[Nr*(Zion+1)];
    this->negRecTerm[0] = new real_t[Nr*(Zion+1)];
    this->posCXTerm[0] = new real_t[Nr*(Zion+1)];
    this->negCXTerm[0] = new real_t[Nr*(Zion+1)];

    for (len_t i = 1; i <= Zion; i++) {
        this->Rec[i]         = this->Rec[i-1] + Nr;
        this->PartialNRec[i] = this->PartialNRec[i-1] + Nr;
        this->PartialTRec[i] = this->PartialTRec[i-1] + Nr;
        this->Ion[i]         = this->Ion[i-1] + Nr;
        this->PartialNIon[i] = this->PartialNIon[i-1] + Nr;
        this->PartialTIon[i] = this->PartialTIon[i-1] + Nr;

        this->posIonizTerm[i] = this->posIonizTerm[i-1] + Nr;
        this->negIonizTerm[i] = this->negIonizTerm[i-1] + Nr;
        this->posRecTerm[i] = this->posRecTerm[i-1] + Nr;
        this->negRecTerm[i] = this->negRecTerm[i-1] + Nr;
        this->posCXTerm[i] = this->posCXTerm[i-1] + Nr;
        this->negCXTerm[i] = this->negCXTerm[i-1] + Nr;
    }
}

/**
 * Deallocate memory for the ionization rate coefficients.
 */
void IonRateEquation::DeallocateRateCoefficients() {
    delete [] this->Ion[0];
    delete [] this->PartialNIon[0];
    delete [] this->PartialTIon[0];
    delete [] this->Rec[0];
    delete [] this->PartialNRec[0];
    delete [] this->PartialTRec[0];

    delete [] this->Ion;
    delete [] this->PartialNIon;
    delete [] this->PartialTIon;
    delete [] this->Rec;
    delete [] this->PartialNRec;
    delete [] this->PartialTRec;

    delete [] this->posIonizTerm[0];
    delete [] this->negIonizTerm[0];
    delete [] this->posRecTerm[0];
    delete [] this->negRecTerm[0];
    delete [] this->posCXTerm[0];
    delete [] this->negCXTerm[0];

    delete [] this->posIonizTerm;
    delete [] this->negIonizTerm;
    delete [] this->posRecTerm;
    delete [] this->negRecTerm;
    delete [] this->posCXTerm;
    delete [] this->negCXTerm;
}

/**
 * Returns an ADAS rate interpolator in the CCD (charge-exchange)
 * coefficient for the ion species with the given ID.
 */
ADASRateInterpolator *IonRateEquation::GetCCD(const len_t iIon) {
    if (this->ions->IsTritium(iIon))
        return this->adas->GetCCD(1, 3);
	else if (this->ions->IsHydrogen(iIon))
		return this->adas->GetCCD(1, 1);
    else
        return this->adas->GetCCD(this->ions->GetZ(iIon));
}

/**
 * Method called whenever the grid is rebuilt.
 */
bool IonRateEquation::GridRebuilt() {
    this->IonEquationTerm<FVM::EquationTerm>::GridRebuilt();

    DeallocateRateCoefficients();
    AllocateRateCoefficients();

    return true;
}

/**
 * Rebuild rate coefficients.
 */
void IonRateEquation::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const len_t Nr = this->grid->GetNr();

    real_t *T = unknowns->GetUnknownData(id_T_cold);
    real_t *n = unknowns->GetUnknownData(id_n_cold);

    ADASRateInterpolator *acd = adas->GetACD(Zion);
    ADASRateInterpolator *scd = adas->GetSCD(Zion);

    // Iterate over charge state (0 ... Z)
    for (len_t i = 0; i < Nr; i++){
        for (len_t Z0 = 0; Z0 <= Zion; Z0++){
            Rec[Z0][i]         = acd->Eval(Z0, n[i], T[i]);
            PartialNRec[Z0][i] = acd->Eval_deriv_n(Z0, n[i], T[i]);
            PartialTRec[Z0][i] = acd->Eval_deriv_T(Z0, n[i], T[i]);
            Ion[Z0][i]         = 0; 
            PartialNIon[Z0][i] = 0;
            PartialTIon[Z0][i] = 0;

            posIonizTerm[Z0][i] = 0;
            negIonizTerm[Z0][i] = 0;
            posRecTerm[Z0][i] = 0;
            negRecTerm[Z0][i] = 0;
            posCXTerm[Z0][i] = 0;
            negCXTerm[Z0][i] = 0;
        }
    }
    // if not covered by the kinetic ionization model, set fluid ionization rates
    if(addFluidIonization || addFluidJacobian)
        for (len_t i = 0; i < Nr; i++){
            for (len_t Z0 = 0; Z0 <= Zion; Z0++){
                Ion[Z0][i]         = scd->Eval(Z0, n[i], T[i]);
                PartialNIon[Z0][i] = scd->Eval_deriv_n(Z0, n[i], T[i]);
                PartialTIon[Z0][i] = scd->Eval_deriv_T(Z0, n[i], T[i]);
            }
        }
}


/**
 * Build block of Jacobian matrix for the given charge state.
 * If accounting for ionization via the kinetic ionization term,
 * the jacobian may still be set here if addFluidJacobian as a 
 * computationally less expensive approximation.
 *
 * derivId: ID of unknown quantity with respect to which differentiation
 *          should be carried out.
 * uqtyId:  ID of unknown quantity to differentiate.
 * jac:     Jacobian matrix to build.
 * x:       Current value of the unknown quantity.
 * iIon:    Index of ion to build jacobian for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
bool IonRateEquation::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac,
    const real_t* nions,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    bool contributes = (derivId==uqtyId);
    if (derivId == uqtyId) 
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset, JACOBIAN);

    #define NI(J,V) \
        jac->SetElement(\
            rOffset+ir, ir, \
            (V) * nions[rOffset+ir+(J)*Nr] \
        )
    #define NI_Z(IZ,J,V) \
        jac->SetElement( \
            rOffset+ir, (IZ)*Nr+ir,\
            (V) * nions[rOffset+ir+(J)*Nr] \
        )
    bool setIonization = addFluidIonization || addFluidJacobian;
    const len_t Nr = this->grid->GetNr();

    if(derivId == id_T_cold) {
        contributes = true;
        #include "IonRateEquation.setDT.cpp"
    } else if(derivId == id_n_cold){
        contributes = true;
        #include "IonRateEquation.setDN.cpp"        
    } else if(derivId == id_Wi){
        #undef NI
        #define NI(J,V) NI_Z(iIon,(J),(V))

        contributes = true;
        #include "IonRateEquation.setDWI.cpp"        
    } else if (derivId == id_lambda_i) {
        #undef NI
        #define NI(J,V) NI_Z(iIon,(J),(V))

        contributes = true;
        #include "IonRateEquation.setDL.cpp"
    } else if(derivId == id_Ni){
        #undef NI
        #define NI(J,V) NI_Z(iIon,(J),(V))

        contributes = true;
        #include "IonRateEquation.setDTINI.cpp"        
    } else if (derivId == uqtyId) {     // Cross-terms from n_i
        #undef NI
        #define NI(J,V) NI_Z(iIon,(J),(V))

        contributes = true;
        #include "IonRateEquation.setDNI.cpp"
    }
    
    #undef NI_Z
    #undef NI

    return contributes;
}

/**
 * Build linear operator matrix for this equation.
 *
 * mat:     Linear operator matrix.
 * rhs:     Vector representing the equation right-hand-side.
 * iIon:    Index of ion to build matrix for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonRateEquation::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset, SetMode sm
) {
    bool setIonization = addFluidIonization || (sm==JACOBIAN&&addFluidJacobian);
    const real_t *nions = this->unknowns->GetUnknownData(id_ions);
    #define NI(J,V,DIAG) \
        do { \
            mat->SetElement(\
                rOffset+ir, rOffset+ir+(J)*Nr, \
                (V) \
            ); \
            if (sm!=JACOBIAN) this->DIAG ## Term[Z0][ir] += (V)*nions[rOffset+ir+(J)*Nr]; \
        } while (false)
    #   include "IonRateEquation.set.cpp"
    #undef NI
}

/**
 * Build function vector.
 *
 * vec:     Function vector to set elements of.
 * nions:   Ion densities.
 * iIon:    Index of ion to build matrix for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonRateEquation::SetCSVectorElements(
    real_t *vec, const real_t *nions,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    bool setIonization = addFluidIonization;
    #define NI(J,V,DIAG) \
        do { \
            vec[rOffset+ir] += (V) * nions[rOffset+ir+(J)*Nr]; \
            this->DIAG ## Term[Z0][ir] += (V) * nions[rOffset+ir+(J)*Nr]; \
        } while (false)
    #   include "IonRateEquation.set.cpp"
    #undef NI
}

