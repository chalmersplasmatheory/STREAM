#include "STREAM/Equations/ConfinementTime.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "STREAM/Grid/EllipticalRadialGridGenerator.hpp"
//#include "STREAM/Equations/PlasmaVolume.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "STREAM/EquationSystem.hpp"
//#include "STREAM/Settings/SimulationGenerator.hpp"
#include <cmath>
//#include "STREAM/Equations/Coefficients.hpp"

using namespace STREAM;
using namespace DREAM;
using namespace std;

 
/**
 * Constructor
 */
ConfinementTime::ConfinementTime(
	FVM::UnknownQuantityHandler *u, EllipticalRadialGridGenerator *r,
	IonHandler *ions, ConnectionLength* CL, STREAM::EquationSystem *eqsys,  len_t D_index) //PlasmaVolume *V,
{
    unknowns = u;
    radials  = r;
    this->ions = ions;
    this->CL = CL;
    this->eqsys = eqsys;
    this->D_index = D_index;
    
    //STREAM::SimulationGenerator::DefineOptions_ConfinementTime(eqsys->GetSettings());
    
    id_Tcold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_T_COLD);
    id_Wi    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_WI_ENER);
    id_Ni    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_NI_DENS);
    id_Ip    = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_I_P);
    id_Ncold = unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_N_COLD);
    id_Efield= unknowns->GetUnknownID(DREAM::OptionConstants::UQTY_E_FIELD);
    
    //type = (enum STREAM::OptionConstants::Conf_Time_type)eqsys->GetSettings()->GetInteger("eqsys/tau_perp/tau_perp");
}

/**
 * Evaluates the inverted confinement time
 */
real_t ConfinementTime::EvaluateConfinementTime(len_t ir){ 
    return EvaluatePerpendicularConfinementTime(ir)
        + EvaluateParallelConfinementTime(ir);
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

const PowerList ConfinementTime::GetCoeff()
{
	return COEFFICIENTS::LawCoefficients.at(type);
}

real_t ConfinementTime::GetOhmicPower(len_t ir)
{
	real_t V = eqsys->GetPlasmaVolume()->GetPlasmaVolume();
	real_t E_field = unknowns->GetUnknownData(id_Efield)[ir];	
	real_t sigma = eqsys->GetREFluid()->GetElectricConductivity(ir);
	real_t R0 = radials->GetMajorRadius();
	
	return E_field * E_field * sigma * V / R0;
}

real_t ConfinementTime::Goldstone_scaling(std::array<const double, 9> const& coeff, len_t ir, real_t const& ConvUnit)
{
	constexpr real_t ConversionFactor = 1e-6; // Convert in Mega for Ampere and Watt
	
	real_t a = radials->GetMinorRadius();
	real_t B0 = radials->GetMagneticField();
	real_t R0 = radials->GetMajorRadius();
	real_t kappa = radials->GetElongation();
	
	real_t Ip = ConversionFactor * unknowns->GetUnknownData(id_Ip)[ir];
		
	real_t ne = ConvUnit * unknowns->GetUnknownData(id_Ncold)[ir]; 
	
	real_t Pm = GetOhmicPower(ir) * ConversionFactor;
	
	real_t A = 2.0141017926; // ATOMIC MASS OF DEUTERIUM ----> Possible Improvement using IonHandler
	
	real_t prod = coeff[0];
	
	const std::array<const real_t, 8> val = {B0, Ip, ne, a, R0, kappa, A, Pm};
	for (size_t i = 0; i < val. size(); ++i)
	{
		prod *= pow(val[i], coeff[i+1]);
	}
	//return C * pow(B0,alpha_1) * pow(Ip,alpha_2) * pow(ne,alpha_3) * pow(a,alpha_4) * pow(R0,alpha_5) * pow(kappa,alpha_6) * pow(A,alpha_7) * pow(Pm, alpha_8);
	return prod;
}

real_t ConfinementTime::Bohm_ConfinementTime(len_t ir)
{
	real_t a      = radials->GetMinorRadius();
	real_t B      = radials->GetMagneticField();
	real_t T_cold = unknowns->GetUnknownData(id_Tcold)[ir];
	return T_cold / (8*a*a*B);
}

real_t ConfinementTime::INTOR_ConfinementTime(len_t ir)
{
	real_t q = eqsys->GetSettings()->GetReal("timestep/safetyfactor");
	real_t tau_0 = q*q;
	real_t C0 = 5e-19;
	real_t ne = unknowns->GetUnknownData(id_Ncold)[ir];
	real_t a = radials->GetMinorRadius();
	return 1.0 / (C0 * tau_0 * ne * a * a);
}

real_t ConfinementTime::ITER89_ConfinementTime(len_t ir, real_t const& ConvUnit)
{
	/*
	real_t C = 0.048;
	real_t alpha_B0 = 0.2;
	real_t alpha_Ip = 0.85;
	real_t alpha_ne = 0.1;
	real_t alpha_a = 0.3;
	real_t alpha_R0 = 1.2;
	real_t alpha_kappa = 0.5;
	real_t alpha_A = 0.5;
	real_t alpha_Pm = -0.5;
	
	return 1.0 / Goldstone_scaling(C, alpha_B0, alpha_Ip, alpha_ne, alpha_a, alpha_R0, alpha_kappa, alpha_A, alpha_Pm, ir, ConvUnit);
	*/
	return 1.0 / Goldstone_scaling(GetCoeff(), ir, ConvUnit);
}

real_t ConfinementTime::ITER97_ConfinementTime(len_t ir)
{
	/*
	real_t C = 0.023;
	real_t alpha_B0 = 0.03;
	real_t alpha_Ip = 0.96;
	real_t alpha_ne = 0.4;
	real_t alpha_a = -0.06;
	real_t alpha_R0 = 1.89;
	real_t alpha_kappa = 0.64;
	real_t alpha_A = 0.20;
	real_t alpha_Pm = -0.73;
	*/
	return 1.0 / Goldstone_scaling(GetCoeff(), ir);
}

real_t ConfinementTime::IPB98_ConfinementTime(len_t ir, real_t const& ConvUnit)
{
	/*
	real_t C = 0.145;
	real_t alpha_B0 = 0.15;
	real_t alpha_Ip = 0.93;
	real_t alpha_ne = 0.41;
	real_t alpha_a = 0.58;
	real_t alpha_R0 = 1.39;
	real_t alpha_kappa = 0.78;
	real_t alpha_A = 0.19;
	real_t alpha_Pm = -0.69;
	*/
	return 1.0 / Goldstone_scaling(GetCoeff(), ir, ConvUnit);
}

real_t ConfinementTime::Goldstone_ConfinementTime(len_t ir)
{
	/*
	real_t C = 0.03;
	real_t alpha_B0 = 0.0;
	real_t alpha_Ip = 1.0;
	real_t alpha_ne = 0.0;
	real_t alpha_a = -0.37;
	real_t alpha_R0 = 1.75;
	real_t alpha_kappa = 0.5;
	real_t alpha_A = 0.5;
	real_t alpha_Pm = -0.5;
	*/
	//return 1.0 / Goldstone_scaling(C, alpha_B0, alpha_Ip, alpha_ne, alpha_a, alpha_R0, alpha_kappa, alpha_A, alpha_Pm, ir, 1e-20);
	return 1.0 / Goldstone_scaling(GetCoeff(), ir, 1e-20);
}

real_t ConfinementTime::ITER89_OL_ConfinementTime(len_t ir)
{
	constexpr real_t ConversionFactor = 1e-6;
	real_t nullset = 0.0;
	/*
	real_t C = 0.064;
	real_t alpha_B0 = 0.35;
	real_t alpha_Ip = 0.8;
	real_t alpha_ne = 0.6;
	real_t alpha_a = 0.6;
	real_t alpha_R0 = 1.6;
	real_t alpha_kappa = 0.5;
	real_t alpha_A = 0.2;
	real_t alpha_Pm = nullset;*/
		
	real_t P = GetOhmicPower(ir) * ConversionFactor;
	real_t WOH = Goldstone_scaling(GetCoeff(), ir, 1e-20);
	real_t tau_inc = Goldstone_scaling({0.04, nullset, 0.5, nullset, 0.8, 0.3, 0.6, 0.5, nullset}, ir, 1e-20);
	
	return 1.0 / (WOH / P + tau_inc);
}

real_t ConfinementTime::EIV1_ConfinementTime(len_t ir)
{
	//return 1.0 / Goldstone_scaling(0.041, 0.16, 0.73, 0.43, -2.18, 1.97, 1.14, 0.39, -0.66, ir);
	return 1.0 / Goldstone_scaling(GetCoeff(), ir);
}

real_t ConfinementTime::EIV2_ConfinementTime(len_t ir)
{
	/*
	real_t C        = 2e-2; 
	real_t alpha_B0 = 0.07; real_t alpha_a     = -3.0;
	real_t alpha_Ip = 0.64; real_t alpha_R0    = 2.12;
	real_t alpha_ne = 0.29; real_t alpha_kappa = 1.35;
	real_t alpha_A  = 0.62; real_t alpha_Pm    = -0.4;
		
	return 1.0 / Goldstone_scaling(C, alpha_B0, alpha_Ip, alpha_ne, alpha_a, alpha_R0, alpha_kappa, alpha_A, alpha_Pm, ir);
	*/
	return 1.0 / Goldstone_scaling(GetCoeff(), ir);
}

real_t ConfinementTime::Kaye_Big_ConfinementTime(len_t ir)
{
	/*
	real_t C        = 0.105; 
	real_t alpha_B0 = 0.3;  real_t alpha_a     = 0.8;
	real_t alpha_Ip = 0.85; real_t alpha_R0    = 0.5;
	real_t alpha_ne = 0.1;  real_t alpha_kappa = 0.25;
	real_t alpha_A  = 0.5;  real_t alpha_Pm    = -0.5;*/
		
	return 1.0 / Goldstone_scaling(GetCoeff(), ir, 1e-20);
}

real_t ConfinementTime::CY_ConfinementTime(len_t ir)
{
	/*
	real_t C        = 0.125; 
	real_t alpha_B0 = 0.35; real_t alpha_a     = 1.1;
	real_t alpha_Ip = 0.65; real_t alpha_R0    = 0.4;
	real_t alpha_ne = 0.15; real_t alpha_kappa = 0.3;
	real_t alpha_A  = 0.5;  real_t alpha_Pm    = -0.5;*/
		
	return 1.0 / Goldstone_scaling(GetCoeff(), ir, 1e-20);
}

real_t ConfinementTime::OS_OL_ConfinementTime(len_t ir)
{
	constexpr real_t ConversionFactor = 1e-6;
	real_t nullset  = 0.0;
	/*real_t C        = 0.064; real_t Zeff        = ions->GetZeff(ir); 
	real_t alpha_B0 = 0.2;   real_t alpha_a     = 0.4;
	real_t alpha_Ip = 1.0;   real_t alpha_R0    = 1.6;
	real_t alpha_ne = 0.6;   real_t alpha_kappa = 0.2;
	real_t alpha_A  = 0.5;   real_t alpha_Pm    = nullset;*/
	
	real_t Zeff        = ions->GetZeff(ir);
	
	real_t q = eqsys->GetSettings()->GetReal("timestep/safetyfactor");
	real_t gq = pow(3 * q * (q + 5) / ((q + 2) * (q + 7)), 0.6);
	real_t fZ = pow(Zeff, 0.4) * pow((15.0-Zeff)/20.0, 0.6);
	real_t P = GetOhmicPower(ir) * ConversionFactor;
	real_t WOH = fZ*gq*Goldstone_scaling(GetCoeff(), ir, 1e-20);
	real_t tau_inc = Goldstone_scaling({0.085, nullset, nullset, nullset, 2.0, nullset, 1.0, 0.5, nullset}, ir, 1e-20);
	
	return 1.0 / (WOH / P + tau_inc);
}

real_t ConfinementTime::RL_OL_ConfinementTime(len_t ir)
{
	constexpr real_t ConversionFactor = 1e-6;
	real_t nullset  = 0.0;
	/*
	real_t C        = 0.17; real_t Zeff        = ions->GetZeff(ir); 
	real_t alpha_B0 = 0.5;  real_t alpha_a     = 2 * 2.75 / 3.0;
	real_t alpha_Ip = 0.5;  real_t alpha_R0    = 2.75 / 3.0;
	real_t alpha_ne = 0.75; real_t alpha_kappa = 2.75 / 3.0;
	real_t alpha_A  = 0.5;  real_t alpha_Pm    = nullset;*/
	
	
	real_t Zeff        = ions->GetZeff(ir); 
	real_t fZ = pow(Zeff, 0.25);
	real_t P = GetOhmicPower(ir) * ConversionFactor;
	real_t WOH = fZ*Goldstone_scaling(GetCoeff(), ir, 1e-20);
	real_t tau_inc = pow(Zeff, -0.5) * Goldstone_scaling({0.014, nullset, 1.0, nullset, 1.0, 0.5, 0.5, 0.5, nullset}, ir, 1e-20);
	
	return 1.0 / (WOH / P + tau_inc);
}

real_t ConfinementTime::EvaluatePerpendicularConfinementTime(len_t ir) 
{
	//type = (enum STREAM::OptionConstants::Conf_Time_type)eqsys->GetSettings()->GetInteger("eqsys/tau_perp/tau_perp");
	switch(type) 
	{
		case STREAM::OptionConstants::CONF_TIME_BOHM :
			return Bohm_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_INTOR :
			return INTOR_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_ITER89 :
			return ITER89_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_ITER97 :
			return ITER97_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_IPB98 :
			return IPB98_ConfinementTime(ir, 1e-20);
		case STREAM::OptionConstants::CONF_TIME_ITER89_OL :
			return ITER89_OL_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_EIV1 :
			return EIV1_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_EIV2 :
			return EIV2_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_GOLDSTONE :
			return Goldstone_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_KAYE_BIG :
			return Kaye_Big_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_CY :
			return CY_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_OS_OL :
			return OS_OL_ConfinementTime(ir);
		case STREAM::OptionConstants::CONF_TIME_RL_OL :
			return RL_OL_ConfinementTime(ir);
		default :
			throw DREAM::SettingsException("Unrecognized equation type for '%s' : %d", "tau_perp", type);
	}
}

real_t ConfinementTime::KappaOut()
{
	return radials->GetElongation();
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

