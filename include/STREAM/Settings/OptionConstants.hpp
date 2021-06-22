#ifndef _STREAM_OPTION_CONSTANTS_HPP
#define _STREAM_OPTION_CONSTANTS_HPP

namespace STREAM {
    class OptionConstants {
    public:
        #include "OptionConstants.enum.hpp"
        
        // CONSTANTS
        // When adding a new constant here, remember
        // to also define it in 'src/Settings/Constants.cpp'.
        // Please, also maintain alphabetical order.
        static const char
            *UQTY_LAMBDA_I;
            
        // Descriptions of unknown quantities
        static const char
            *UQTY_LAMBDA_I_DESC;
    };
}
#endif/*_STREAM_OPTION_CONSTANTS_HPP*/
