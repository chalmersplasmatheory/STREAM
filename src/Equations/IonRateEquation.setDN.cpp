/**
 * Sets electron-density jacobian of IonRateEquation
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    const real_t dV_ndN = this->volumes->GetNeutralVolume_dn(iIon);
    const real_t dV_n_totdN = this->volumes->GetTotalNeutralVolume_dn(iIon);
    
    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // d/dn_cold[I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)]
            if (Z0 == 1){
                NI(-1, Ion[Z0-1][ir] * V_n/V_p + PartialNIon[Z0-1][ir] * n_cold[ir] * V_n /V_p + Ion[Z0-1][ir] * n_cold[ir]/V_p * dV_ndN ); 
            } else if (Z0 > 1){
                NI(-1, Ion[Z0-1][ir] + PartialNIon[Z0-1][ir] * n_cold[ir]);
            }   

            // d/dn_cold[-I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)]
            if (Z0 == 0){
                NI(0, -Ion[Z0][ir] * V_n/V_n_tot - PartialNIon[Z0][ir] * n_cold[ir] * V_n /V_n_tot - Ion[Z0][ir] * n_cold[ir]/V_n_tot * dV_ndN  + Ion[Z0][ir] * n_cold[ir] * V_n /(V_n_tot*V_n_tot) * dV_n_totdN );
            }else{
                NI(0, -Ion[Z0][ir] - PartialNIon[Z0][ir] * n_cold[ir] );
            }
        }
        
        // d/dn_cold[R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)]
        if (Z0 == 0){
            NI(+1, Rec[Z0+1][ir] * V_p/V_n_tot + PartialNRec[Z0+1][ir] * n_cold[ir] * V_p /V_n_tot - Rec[Z0+1][ir] * n_cold[ir] * V_p/(V_n_tot*V_n_tot) * dV_n_totdN ); 
        } else if (Z0 < Z){
            NI(+1, Rec[Z0+1][ir] + PartialNRec[Z0+1][ir] * n_cold[ir] );
        }
        
        // d/dn_cold[-R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)]
        if (Z0 > 0){                        
            NI(0, -Rec[Z0][ir] - PartialNRec[Z0][ir] * n_cold[ir] );
        }
    }
