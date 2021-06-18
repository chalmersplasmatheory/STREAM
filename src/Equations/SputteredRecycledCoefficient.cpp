#include "STREAM/Equations/SputteredRecycledCoefficient.hpp"

using namespace STREAM;
//using namespace DREAM;
using namespace std;

/**
 * Constructor
 */
SputteredRecycledCoefficient::SputteredRecycledCoefficient(){} // Denna ska vara tom va? Det enda vi vill ha är tabellen

void SputteredRecycledCoefficient::AddSRCoefficient(len_t upper, len_t lower, real_t coefficient){
    if (coefficientTable.find(upper) != coefficientTable.end()){ // Ska det vara . innan vissa funktioner?
        if( coefficientTable.find(upper)->second.find(lower) != coefficientTable.find(upper)->second.end()){ //Är detta rätt sätt att nå den inre unordered_map?
            coefficientTable.find(upper)->second.erase(lower);
        }
        coefficientTable.find(upper)->second.insert(pair<len_t,real_t>(lower,coefficient)); //Är detta rätt sätt att lägga till i den inre unordered_map?
    } else {
        unordered_map<len_t,real_t> innerUM;
        innerUM.insert(pair<len_t,real_t>(lower,coefficient));
        coefficientTable.insert(pair<len_t, unordered_map<len_t,real_t>>(upper,innerUM));//Är detta rätt sätt att lägga till ett till element??
    }
}

real_t SputteredRecycledCoefficient::GetSRCoefficient(len_t upper, len_t lower){
    if (coefficientTable.find(upper) != coefficientTable.end() && coefficientTable.find(upper)->second.find(lower) != coefficientTable.find(upper)->second.end()){ // Korrekt?
        return coefficientTable.find(upper)->second.find(lower)->second;
    } else {
        return 0; 
    }
}
