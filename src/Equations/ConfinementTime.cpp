#include "STREAM/Equations/ConfinementTime.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "STREAM/EquationSystem.hpp"
#include <cmath>
#include <algorithm>
#include <functional>

using namespace STREAM;
using namespace DREAM;
using namespace std;

 
/**
 * Constructor
 */

ConfinementTime::ConfinementTime(
	FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
	IonHandler *ions, ConnectionLength* CL, STREAM::EquationSystem *eqsys,  len_t D_index) {
    unknowns = u;
    radials  = r;
    
    this->ions = ions;
    this->CL = CL;
    this->eqsys = eqsys;
    this->D_index = D_index;
    
    
    id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    id_Wi    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    id_Ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    id_Ip    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    id_Ncold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    id_Efield= unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
    
    type = (enum STREAM::OptionConstants::Conf_Time_type)eqsys->GetSettings()->GetInteger("eqsys/tau_perp/tau_perp"); //Setting the type of confinement time
	maxPerpLaw = (bool)eqsys->GetSettings()->GetInteger("eqsys/tau_perp/maxPerpLaw"); // If a combination of confinement law is activated
}
/**
 * Evaluates the inverted confinement time
 */
real_t ConfinementTime::EvaluateConfinementTime(len_t ir) { 
    return EvaluatePerpendicularConfinementTime(ir) + EvaluateParallelConfinementTime(ir);
}

/**
 * Evaluates the parallel confinement time.
 */
real_t ConfinementTime::EvaluateParallelConfinementTime(len_t ir) {
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);
    
    return Lfinv * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/mi);
}

/**
 * Evaluates the perpendicular confinement time.
 */

/**
 * Returns a power list containing the coefficients for each confinement law
 * stored in a map in a .hpp file.
 *
 * _type: Type of the confinement time.
 * 
 * Returns: PowerList containing the coefficients.
 */
const PowerList ConfinementTime::GetCoeff(STREAM::OptionConstants::Conf_Time_type _type) const {
	return STREAM::LawCoefficients.at(_type);
}

/**
 * Returns the power list containing the coefficients for the current confinement law.
 *
 * Returns: PowerList containing the coefficients.
 */
const PowerList ConfinementTime::GetCoeff() const {
	return STREAM::LawCoefficients.at(type);
}

/**
 * Computes the ohmic power for a given radial index.
 *
 * ir: Radial index.
 *
 * Returns: Computed ohmic power.
 */
real_t ConfinementTime::GetOhmicPower(len_t ir) const {
	real_t V = eqsys->GetPlasmaVolume()->GetPlasmaVolume();
	real_t E_field = unknowns->GetUnknownData(id_Efield)[ir];	
	real_t sigma = eqsys->GetREFluid()->GetElectricConductivity(ir);
	real_t R0 = radials->GetMajorRadius();
	
	return E_field * E_field * sigma * V / R0; 
}

/**
 * Computes the Goldstone power scaling law.
 *
 * coeff: PowerList containing coefficients.
 * ir: Radial index.
 * ConvUnit: Conversion unit.
 *
 * Returns: Computed Goldstone power scaling.
 */
real_t ConfinementTime::Goldstone_scaling(const PowerList& coeff, len_t ir, 
		real_t const& ConvUnit) const {
	real_t a = radials->GetMinorRadius();
	real_t B0 = radials->GetMagneticField();
	real_t R0 = radials->GetMajorRadius();
	real_t kappa = radials->GetElongation();
	real_t Ip = ConversionFactor * unknowns->GetUnknownData(id_Ip)[ir];
	real_t ne = ConvUnit * unknowns->GetUnknownData(id_Ncold)[ir];
	real_t Pm = GetOhmicPower(ir) * ConversionFactor;
	// Atomic mass of deuterium
	real_t A = DREAM::Constants::mD / DREAM::Constants::mu; 

	real_t prod = coeff[0];
	const std::array<const real_t, 8> val = {B0, Ip, ne, a, R0, kappa, A, Pm};
	for (size_t i = 0; i < val.size(); ++i) 	{
		prod *= pow(val[i], coeff[i+1]);
	}
	return prod;
}

/**
 * Computes Bohm confinement time.
 *
 * ir: Radial index.
 *
 * Returns: Bohm confinement time.
 */
real_t ConfinementTime::Bohm_ConfinementTime(len_t ir) const  {
	real_t a = radials->GetMinorRadius();
	real_t B = radials->GetMagneticField();
	real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
	return T_cold / (8*a*a*B);
}

/**
 * Computes INTOR confinement time.
 *
 * [R. J. Goldston, PPCF 26 (1984) https://doi.org/10.1088/0741-3335/26/1A/308]
 * ir: Radial index.
 *
 * Returns: INTOR confinement time.
 */
real_t ConfinementTime::INTOR_ConfinementTime(len_t ir) const  {
	real_t mu0Ip = DREAM::Constants::mu0 * this->unknowns->GetUnknownData(id_Ip)[0];
	real_t qR0 = eqsys->GetFluidGrid()->GetRadialGrid()->SafetyFactorNormalized(0, mu0Ip);
	real_t R0 = radials->GetMajorRadius();
	real_t tau_0 = sqrt(qR0 / R0);
	real_t C0 = 5e-21;
	real_t ne = unknowns->GetUnknownData(id_Ncold)[ir];
	real_t a = radials->GetMinorRadius();
	return 1.0 / (C0 * tau_0 * ne * a * a);
}

/**
 * Computes RLW confinement time.
 *
 * ir: Radial index.
 *
 * Returns: RLW confinement time.
 */
real_t ConfinementTime::RLW_ConfinementTime(len_t ir) const {
	real_t temp = Goldstone_scaling(GetCoeff(), ir);
	real_t Zeff = ions->GetZeff(ir);
	real_t a = radials->GetMinorRadius();
	real_t R0 = radials->GetMajorRadius();
	real_t kappa = radials->GetElongation();
	real_t Ip = ConversionFactor * unknowns->GetUnknownData(id_Ip)[ir];

	return 1.0 / (temp * sqrt(sqrt(Zeff)) + 0.012 * Ip * sqrt(R0) * a * sqrt(kappa) / sqrt(Zeff));
}

/**
 * Computes ITER89 offset linear confinement time.
 *
 * ir: Radial index.
 *
 * Returns: ITER89 OL confinement time.
 */
real_t ConfinementTime::ITER89_OL_ConfinementTime(len_t ir) const {		
	real_t P = GetOhmicPower(ir) * ConversionFactor;
	real_t WOH = Goldstone_scaling(GetCoeff(), ir, 1e-20);
	real_t tau_inc = Goldstone_scaling({0.04, nullset, 0.5, nullset, 0.8, 0.3, 0.6, 0.5, nullset}, ir, Conv_n20);
	
	return 1.0 / (WOH / P + tau_inc);
}

/**
 * Computes OS OL confinement time.
 *
 * ir: Radial index.
 *
 * Returns: OS OL confinement time.
 */
real_t ConfinementTime::OS_OL_ConfinementTime(len_t ir) const {
	real_t Zeff = ions->GetZeff(ir);
	real_t mu0Ip = DREAM::Constants::mu0 * this->unknowns->GetUnknownData(id_Ip)[0];
	real_t qR0 = eqsys->GetFluidGrid()->GetRadialGrid()->SafetyFactorNormalized(0, mu0Ip);
	real_t R0 = radials->GetMajorRadius();
	real_t q = qR0 / R0;
	real_t gq = pow(3 * q * (q + 5) / ((q + 2) * (q + 7)), 0.6);
	real_t fZ = pow(Zeff, 0.4) * pow((15.0 - Zeff) / 20.0, 0.6);
	real_t P = GetOhmicPower(ir) * ConversionFactor;
	real_t WOH = fZ * gq * Goldstone_scaling(GetCoeff(), ir, 1e-20);
	real_t tau_inc = Goldstone_scaling({0.085, nullset, nullset, nullset, 2.0, nullset, 1.0, 0.5, nullset}, ir, Conv_n20);

	return 1.0 / (WOH / P + tau_inc);
}

/**
 * Computes RL OL confinement time.
 *
 * ir: Radial index.
 *
 * Returns: RL OL confinement time.
 */
real_t ConfinementTime::RL_OL_ConfinementTime(len_t ir) const {
	real_t Zeff = ions->GetZeff(ir);
	real_t fZ = sqrt(sqrt(Zeff));
	real_t P = GetOhmicPower(ir) * ConversionFactor;
	real_t WOH = fZ * Goldstone_scaling(GetCoeff(), ir, 1e-20);
	real_t tau_inc = pow(Zeff, -0.5) * Goldstone_scaling({0.014, nullset, 1.0, nullset, 1.0, 0.5, 0.5, 0.5, nullset}, ir, Conv_n20);

	return 1.0 / (WOH / P + tau_inc);
}

/**
 * Evaluates the confinement time based on the confinement type.
 *
 * ir: Radial index.
 *
 * Returns: Evaluated perpendicular confinement time.
 */
const real_t ConfinementTime::EvaluatePerpendicularConfinementTimeType(len_t ir) const {
	switch(type) {
		case STREAM::OptionConstants::CONF_TIME_BOHM:
			return Bohm_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_INTOR:
			return INTOR_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_RLW:
			return RLW_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_ITER89:
		case STREAM::OptionConstants::CONF_TIME_IPB98:
		case STREAM::OptionConstants::CONF_TIME_GOLDSTONE:
		case STREAM::OptionConstants::CONF_TIME_KAYE_BIG:
		case STREAM::OptionConstants::CONF_TIME_CY:
			return 1.0 / Goldstone_scaling(GetCoeff(), ir, Conv_n20);
		case STREAM::OptionConstants::CONF_TIME_ITER97:
		case STREAM::OptionConstants::CONF_TIME_EIV1:
		case STREAM::OptionConstants::CONF_TIME_EIV2:
			return 1.0 / Goldstone_scaling(GetCoeff(), ir);
		case STREAM::OptionConstants::CONF_TIME_ITER89_OL:
			return ITER89_OL_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_OS_OL:
			return OS_OL_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_RL_OL:
			return RL_OL_ConfinementTime(ir);
		default:
			throw DREAM::SettingsException("Unrecognized equation type for '%s' : %d", "tau_perp", type);
	}
}

/**
 * Evaluates the perpendicular confinement time, optionally taking the maximum between the current 
 * scaling law and the Bohm confinement time.
 *
 * ir: Radial index.
 *
 * Returns: Evaluated perpendicular confinement time.
 */
real_t ConfinementTime::EvaluatePerpendicularConfinementTime(len_t ir) {
	real_t bohmTime = Bohm_ConfinementTime(ir);
	real_t perpTime = EvaluatePerpendicularConfinementTimeType(ir);

	if (maxPerpLaw) {
		if (bohmTime < perpTime) {
			currentTypeUsed = STREAM::OptionConstants::CONF_TIME_BOHM;
			maxPerpLawUsed = true;
		} 
		else {
			currentTypeUsed = type;
			maxPerpLawUsed = false;
		}
		return std::min(bohmTime, perpTime);
	}
	else {
		currentTypeUsed = type;
		maxPerpLawUsed = false;
		return perpTime;
	}
}

/**
 * Returns the plasma elongation.
 *
 * Returns: Plasma elongation (kappa).
 */
real_t ConfinementTime::KappaOut() {
	return radials->GetElongation();
}

/**
 * Returns the current confinement type being used.
 *
 * Returns: Current confinement type.
 */
real_t ConfinementTime::ConfTimeTypeOut() {
	return currentTypeUsed;
}

/**
 * Returns whether the maximum perpendicular confinement law is being used.
 *
 * Returns: Boolean indicating if the maximum perpendicular confinement law is used.
 */
real_t ConfinementTime::MaxPerpLawOut() {
	return maxPerpLawUsed;
}


/**
 * Evaluates the derivative of the inverted confinement time with respect to the plasma current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIp(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t ec = DREAM::Constants::ec;

    real_t dLfinv_dIp = this->CL->EvaluateInverseConnectionLength_dIp(ir);
    
    return dLfinv_dIp * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the wall current
 */
real_t ConfinementTime::EvaluateConfinementTime_dIwall(len_t ir){
    len_t nr = radials->GetNr();
    
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);

    real_t ec = DREAM::Constants::ec;

    real_t dLfinv_dIwall = this->CL->EvaluateInverseConnectionLength_dIwall(ir);
    
    return dLfinv_dIwall * sqrt((ec*T_cold+2.0/3.0*W_i/N_i)/mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the electron temperature
 */
real_t ConfinementTime::EvaluateConfinementTime_dTcold(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);
    
    real_t a = radials->GetMinorRadius();
    real_t B = radials->GetMagneticField(); 
    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);

    return 1.0/(8*a*a*B) + 1.0/2.0*ec * Lfinv * 1.0/sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the ion energy
 */
real_t ConfinementTime::EvaluateConfinementTime_dWi(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);
    
    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);

    return Lfinv * 1.0/3.0*1.0/N_i * 1.0 / sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*mi);
}

/**
 * Evaluates the derivative of the inverted confinement time with respect to the total ion density
 */
real_t ConfinementTime::EvaluateConfinementTime_dNi(len_t ir){
    len_t nr = radials->GetNr();
    real_t T_cold    = unknowns->GetUnknownData(id_Tcold)[ir];
    real_t W_i    = unknowns->GetUnknownData(id_Wi)[D_index*nr+ir];
    real_t N_i    = unknowns->GetUnknownData(id_Ni)[D_index*nr+ir];
    real_t mi     = this->ions->GetIonSpeciesMass(this->D_index);
    
    real_t ec = DREAM::Constants::ec;

    real_t Lfinv = this->CL->EvaluateInverseConnectionLength(ir);

    return -Lfinv * 1.0/3.0*W_i/(N_i*N_i) * 1.0 / sqrt((ec*T_cold+2.0/3.0*W_i/N_i)*mi);
}


