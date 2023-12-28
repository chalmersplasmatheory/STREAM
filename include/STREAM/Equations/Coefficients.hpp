#ifndef _STREAM_EQUATIONS_COEFFICIENTS_HPP
#define _STREAM_EQUATIONS_COEFFICIENTS_HPP


#include <map>
#include <array>
#include "STREAM/Settings/OptionConstants.hpp"


typedef std::array<const double,9> PowerList;

namespace COEFFICIENTS
{	
	constexpr double nullset = 0.0;
	
	const std::map<const enum STREAM::OptionConstants::Conf_Time_type, const PowerList> LawCoefficients = 
	{
			// {CONF_TIME_LAW, {C, alpha_B0, alpha_Ip, alpha_ne, alpha_a, alpha_R0, alpha_kappa, alpha_A, alpha_Pm}}
			{STREAM::OptionConstants::CONF_TIME_ITER89, {0.048, 0.2, 0.85, 0.1, 0.3, 1.2, 0.5, 0.5, -0.5}},
			{STREAM::OptionConstants::CONF_TIME_ITER97, {0.0264, 0.03, 0.96, 0.4, -0.06, 1.89, 0.64, nullset, -0.73}},
			{STREAM::OptionConstants::CONF_TIME_IPB98, {0.145, 0.15, 0.93, 0.41,0.58, 1.39, 0.78, 0.19, -0.69}},
			{STREAM::OptionConstants::CONF_TIME_ITER89_OL, {0.064, 0.35, 0.8, 0.6, 0.6, 1.6, 0.5, 0.2, nullset}},
			{STREAM::OptionConstants::CONF_TIME_EIV1, {0.041, 0.16, 0.73, 0.43, -2.18, 1.97, 1.14, 0.39, -0.66}},
			{STREAM::OptionConstants::CONF_TIME_EIV2, {0.02, 0.07, 0.64, 0.29, -3.0, 2.12, 1.35, 0.62, -0.4}},
			{STREAM::OptionConstants::CONF_TIME_GOLDSTONE, {0.03, nullset, 1.0, nullset, -0.37, 1.75, 0.5, 0.5, -0.5}},
			{STREAM::OptionConstants::CONF_TIME_KAYE_BIG, {0.105, 0.3, 0.85, 0.1, 0.8, 0.5, 0.25, 0.5, -0.5}},
			{STREAM::OptionConstants::CONF_TIME_CY, {0.125, 0.35, 0.65, 0.15, 1.1, 0.4, 0.3, 0.5, -0.5}},
			{STREAM::OptionConstants::CONF_TIME_OS_OL, {0.064, 0.2, 1.0, 0.6, 0.4, 1.6, 0.2, 0.5, nullset}},
			{STREAM::OptionConstants::CONF_TIME_RL_OL, {0.17, 0.5, 0.5, 0.75, 2 * 2.75 / 3.0, 2.75 / 3.0, 2.75 / 3.0, 0.5, nullset}}
	};
}

#endif