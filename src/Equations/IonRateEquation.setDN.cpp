/**
 * Sets electron-density jacobian of IonRateEquation
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const real_t *T_cold = this->unknowns->GetUnknownData(id_T_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    const len_t NZ = this->ions->GetNZ();
    
    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // d/dn_cold[I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)]
            if (Z0 == 1){
                NI(-1, Ion[Z0-1][ir] * V_n/V_p + PartialNIon[Z0-1][ir] * n_cold[ir] * V_n /V_p); 
            } else if (Z0 > 1){
                NI(-1, Ion[Z0-1][ir] + PartialNIon[Z0-1][ir] * n_cold[ir]);
            }   

            // d/dn_cold[-I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)]
            if (Z0 == 0){
                NI(0, -Ion[Z0][ir] * V_n/V_n_tot - PartialNIon[Z0][ir] * n_cold[ir] * V_n /V_n_tot);
            }else{
                NI(0, -Ion[Z0][ir] - PartialNIon[Z0][ir] * n_cold[ir] );
            }
        }
        
        // d/dn_cold[R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)]
        if (Z0 == 0){
            NI(+1, Rec[Z0+1][ir] * V_p/V_n_tot + PartialNRec[Z0+1][ir] * n_cold[ir] * V_p /V_n_tot); 
        } else if (Z0 < Z){
            NI(+1, Rec[Z0+1][ir] + PartialNRec[Z0+1][ir] * n_cold[ir] );
        }
        
        // d/dn_cold[-R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)]
        if (Z0 > 0){                        
            NI(0, -Rec[Z0][ir] - PartialNRec[Z0][ir] * n_cold[ir] );
        }
        
        // d/dn_cold(Positive charge-exchange term)
        if (Z == 1){
            if (Z0 == 1){
                for (len_t iz=0; iz<NZ; iz++){ 
                    if(iz==iIon) 
                        continue;
                    len_t Zi = ions->GetZ(iz); 
                    const len_t IonOffset = ions->GetIndex(iz,0); 
                    ADASRateInterpolator *ccd = adas->GetCCD(Zi); 
                    for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ 
                        real_t PartialNRcx = ccd->Eval_deriv_n(Z0i-1, n_cold[ir], T_cold[ir]); 
                        NI(-1, PartialNRcx * V_n/V_p * nions[(IonOffset+Z0i)*Nr+ir]); 
                    }
                }
            }
        }else if (Z0 < Z){  
            for (len_t iz=0; iz<NZ; iz++){ 
                if(ions->GetZ(iz)!=1) 
                    continue;
                const len_t Doffset = ions->GetIndex(iz,0); 
                ADASRateInterpolator *ccd = adas->GetCCD(Z); 
                real_t PartialNRcx = ccd->Eval_deriv_n(Z0+1-1, n_cold[ir], T_cold[ir]); 
                const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                if (Z0 == 0){
                    NI(+1, PartialNRcx * V_n_D/V_n_tot * nions[Doffset*Nr + ir]); 
                }else{
                    NI(+1, PartialNRcx * V_n_D/V_p * nions[Doffset*Nr + ir]);
                }
            }
        }
        
        // d/dn_cold(Negative charge-exchange term)
        if (Z == 1){
            if(Z0 == 0){
                for (len_t iz=0; iz<NZ; iz++){
                    if(iz==iIon)
                        continue;
                    len_t Zi = ions->GetZ(iz); 
                    const len_t IonOffset = ions->GetIndex(iz,0);
                    ADASRateInterpolator *ccd = adas->GetCCD(Zi);
                    for(len_t Z0i=1; Z0i<Zi+1; Z0i++){
                        real_t PartialNRcx = ccd->Eval_deriv_n(Z0i-1, n_cold[ir], T_cold[ir]);
                        NI(0, -PartialNRcx * V_n/V_n_tot * nions[(IonOffset+Z0i)*Nr+ir]); 
                    }
                }
            }
        } else if (Z0 >= 1){  
            for (len_t iz=0; iz<NZ; iz++){ 
                if(ions->GetZ(iz)!=1) 
                    continue;
                const len_t Doffset = ions->GetIndex(iz,0); 
                ADASRateInterpolator *ccd = adas->GetCCD(Z); 
                real_t PartialNRcx = ccd->Eval_deriv_n(Z0-1, n_cold[ir], T_cold[ir]); 
                const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                NI(0, -PartialNRcx * V_n_D/V_p * nions[Doffset*Nr + ir]); 
                
            }
        }
    }
