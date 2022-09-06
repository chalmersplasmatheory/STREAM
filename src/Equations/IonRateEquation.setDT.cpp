/**
 * Sets temperature jacobian of IonRateEquation
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    
    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){ 
            // d/dT_cold[I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)]            
            if (Z0 == 1){
                NI(-1, PartialTIon[Z0-1][ir]*n_cold[ir]*V_n/V_p);
            } else if (Z0 > 1){
                NI(-1, PartialTIon[Z0-1][ir] * n_cold[ir]);
            }
            // d/dT_cold[-I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)]
            if (Z0 == 0){
                NI(0, -n_cold[ir] * PartialTIon[Z0][ir] * V_n/V_n_tot);
            }else if (Z0 < Z){
                NI(0, - PartialTIon[Z0][ir] * n_cold[ir] );
            }           
        }
        // d/dT_cold[R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)]
        if (Z0 == 0){
            NI(+1, n_cold[ir] * V_p * PartialTRec[Z0+1][ir] /V_n_tot); 
        }else if (Z0 < Z) {
            NI(+1, PartialTRec[Z0+1][ir] * n_cold[ir]);
        }

        // d/dT_cold[-R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)]
        if (Z0 > 0){                        
            NI(0, -PartialTRec[Z0][ir] * n_cold[ir]);
        }
        
    }
