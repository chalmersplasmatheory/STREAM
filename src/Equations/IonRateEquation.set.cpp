/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); //Is the argument iIon the right one here?
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    
    //const len_t ionidx = this->ions->GetIndex(iIon, Z0);

    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)
            if (Z0 == 1){
                NI(-1, Ion[Z0-1][ir] * n_cold[ir] * V_n/V_p);
            } else if (Z0 > 1){
                NI(-1, Ion[Z0-1][ir] * n_cold[ir]);
            }

            // -I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
            if (Z0 == 0){
                NI(0, -Ion[Z0][ir] * n_cold[ir] * V_n/V_n_tot);
            }else{
                NI(0, -Ion[Z0][ir] * n_cold[ir]);
            }
        }
        
        // R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)
        if (Z0 == 0){
            NI(+1, Rec[Z0+1][ir] * n_cold[ir] * V_p/V_n_tot);
        } else if (Z0 < Z){
            NI(+1, Rec[Z0+1][ir] * n_cold[ir]);
        }

        // -R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
        // Does not contribute when Z0=0 since there is no recombination for neutrals
        if (Z0 > 0){                        
            NI(0, -Rec[Z0][ir] * n_cold[ir]);
        }
    }
