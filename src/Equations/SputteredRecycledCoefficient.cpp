#include "STREAM/Equations/SputteredRecycledCoefficient.hpp"

using namespace STREAM;
//using namespace DREAM;
using namespace std;

/**
 * Constructor
 */
SputteredRecycledCoefficient::SputteredRecycledCoefficient(std::unordered_map<len_t /*iIon1*/, std::unordered_map<len_t /*iIon2*/, real_t /*Y*/> > *coefficientTable) : coefficientTable(coefficientTable) {}

void SputteredRecycledCoefficient::AddSRCoefficient(len_t upper, len_t lower, real_t coefficient){
    if (coefficientTable->find(upper) != coefficientTable->end()){ 
        if( coefficientTable->find(upper)->second.find(lower) != coefficientTable->find(upper)->second.end()){ 
            coefficientTable->find(upper)->second.erase(lower);
        }
        coefficientTable->find(upper)->second.insert(pair<len_t,real_t>(lower,coefficient)); 
    } else {
        unordered_map<len_t,real_t> innerUM;
        innerUM.insert(pair<len_t,real_t>(lower,coefficient));
        coefficientTable->insert(pair<len_t, unordered_map<len_t,real_t>>(upper,innerUM));
    }
}

real_t SputteredRecycledCoefficient::GetSRCoefficient(len_t upper, len_t lower){
    if (coefficientTable->find(upper) != coefficientTable->end() && coefficientTable->find(upper)->second.find(lower) != coefficientTable->find(upper)->second.end()){ 
        return coefficientTable->find(upper)->second.find(lower)->second;
    } else {
        return 0; 
    }
}
