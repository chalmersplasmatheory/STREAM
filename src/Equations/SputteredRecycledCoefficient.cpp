#include "STREAM/Equations/SputteredRecycledCoefficient"

using namespace STREAM;
using namespace DREAM;
using namespace std;

/**
 * Constructor
 */
SputteredRecycledCoefficient::SputteredRecycledCoefficient(){}

void SputteredRecycledCoefficient::add(len_t upper, len_t lower, real_t coefficient){
    if (coefficientTable.find(upper) != coefficientTable.end()){
        if( coefficientTable.find(upper)->second.find(lower) != coefficientTable.find(upper)->second.end()){
            coefficientTable.find(upper)->second.erase(lower);
        }
        coefficientTable.find(upper)->second.insert({{lower,coefficient}});
    }
    coefficientTable.insert({{upper,{{lower},{koefficient}}}});
}

real_t SputteredRecycledCoefficient::get(len_t upper, len_t lower){
    if (coefficientTable.find(upper) != coefficientTable.end()){
        if (coefficientTable.find(upper)->second.find(lower) != coefficientTable.find(upper)->second.end()){
            return coefficientTable.find(upper)->second.find(lower)->second;
        } else {
            return 0; 
        }
    } else {
        return 0; 
    }
}
