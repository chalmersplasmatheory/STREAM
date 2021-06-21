/**
 * Definition of OptionConstants constants.
 */

#include "STREAM/Settings/OptionConstants.hpp"

using namespace STREAM;

/**
 * NAMES OF UNKNOWN QUANTITIES
 */
const char *OptionConstants::UQTY_LAMBDA_I                = "lambda_i";
const char *OptionConstants::UQTY_ION_HEAT_TRANSPORT      = "Ion heat transport";
const char *OptionConstants::UQTY_ION_TRANSPORT           = "Ion  transport";

/**
 * DESCRIPTIONS OF UNKNOWN QUANTITIES
 */
const char *OptionConstants::UQTY_LAMBDA_I_DESC           = "Mean free path of neutrals of each species [m]";
const char *OptionConstants::UQTY_ION_HEAT_TRANSPORT_DESC = "Transport of ion heat / energy [J/s]";
const char *OptionConstants::UQTY_ION_TRANSPORT_DESC      = "Transport of ions [1/m^3 s]";
