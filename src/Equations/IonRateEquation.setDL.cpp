/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const real_t *T_cold = this->unknowns->GetUnknownData(id_T_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    const real_t dV_n = this->volumes->GetNeutralVolume_dLambdai(iIon); 
    const real_t dV_n_tot = this->volumes->GetTotalNeutralVolume_dLambdai(iIon);
    const len_t NZ = this->ions->GetNZ();

    if (V_n > 0) {
        for (len_t ir = 0; ir < Nr; ir++) {
            if(setIonization){
                // I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)
                if (Z0 == 1){
                    NI(-1, Ion[Z0-1][ir] * n_cold[ir] * dV_n/V_p);
                } else if (Z0 > 1){
                    NI(-1, Ion[Z0-1][ir] * n_cold[ir]); //Isn't this 0 if only V_n and V_n_tot depend on lambdai?
                }

                // -I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
                if (Z0 == 0){
                    NI(0, -Ion[Z0][ir] * n_cold[ir] * (1/V_n - V_n/(V_n_tot*V_n_tot)*dV_n_tot)); //Shouldn't it be 1/V_n_tot in the first term?
                }else{
                    NI(0, -Ion[Z0][ir] * n_cold[ir]); //Isn't this 0 if only V_n and V_n_tot depend on lambdai?
                }
            }
            
            // R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)
            if (Z0 == 0){
                NI(+1, -Rec[Z0+1][ir] * n_cold[ir] * V_p/(V_n_tot*V_n_tot)*dV_n_tot);
            } else if (Z0 < Z){
                NI(+1, Rec[Z0+1][ir] * n_cold[ir]); //Isn't this 0 if only V_n and V_n_tot depend on lambdai?
            }

            // -R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
            // Does not contribute when Z0=0 since there is no recombination for neutrals
            if (Z0 > 0){                        
                NI(0, -Rec[Z0][ir] * n_cold[ir]); //Isn't this 0 if only V_n and V_n_tot depend on lambdai?
            }
            
            // d/dlambda_i(Positive charge-exchange term)
            if (Z == 1){
                if (Z0 == 1){
                    for (len_t iz=0; iz<NZ; iz++){ 
                        if(iz==iIon) 
                            continue;
                        len_t Zi = ions->GetZ(iz); 
                        const len_t IonOffset = ions->GetIndex(iz,0); 
                        ADASRateInterpolator *ccd = adas->GetCCD(Zi); 
                        for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ 
                            real_t Rcx = ccd->Eval(Z0i, n_cold[ir], T_cold[ir]); 
                            NI(-1, Rcx * dV_n/V_p * nions[IonOffset+Z0i*Nr+ir]); 
                        }
                    }
                }
            }else if (Z0 < Z){  
                for (len_t iz=0; iz<NZ; iz++){ 
                    if(ions->GetZ(iz)!=1) 
                        continue;
                    const len_t Doffset = ions->GetIndex(iz,0); 
                    ADASRateInterpolator *ccd = adas->GetCCD(Z); 
                    real_t Rcx = ccd->Eval(Z0+1, n_cold[ir], T_cold[ir]); 
                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz);
                    const real_t dV_n_D = this->volumes->GetNeutralVolume_dLambdai(iz);
                    if (Z0 == 0){
                        NI(+1, Rcx * (dV_n_D/V_n_tot - V_n_D * dV_n_tot /(V_n_tot * V_n_tot)) * nions[Doffset + ir]); 
                    }else{
                        NI(+1, Rcx * dV_n_D/V_p * nions[Doffset + ir]);
                    }
                }
            }
            
            // d/dlambda_i(Negative charge-exchange term)
            if (Z == 1){
                if(Z0 == 0){
                    for (len_t iz=0; iz<NZ; iz++){
                        if(iz==iIon)
                            continue;
                        len_t Zi = ions->GetZ(iz); 
                        const len_t IonOffset = ions->GetIndex(iz,0);
                        ADASRateInterpolator *ccd = adas->GetCCD(Zi);
                        for(len_t Z0i=1; Z0i<Zi+1; Z0i++){
                            real_t Rcx = ccd->Eval(Z0i, n_cold[ir], T_cold[ir]);
                            NI(0, -Rcx * (dV_n/V_n_tot - V_n * dV_n_tot/(V_n_tot*V_n_tot)) * nions[IonOffset+Z0i*Nr+ir]); 
                        }
                    }
                }
            } else if (Z0 > 1){  
                for (len_t iz=0; iz<NZ; iz++){ 
                    if(ions->GetZ(iz)!=1) 
                        continue;
                    const len_t Doffset = ions->GetIndex(iz,0); 
                    ADASRateInterpolator *ccd = adas->GetCCD(Z); 
                    real_t Rcx = ccd->Eval(Z0, n_cold[ir], T_cold[ir]); 
                    const real_t dV_n_D = this->volumes->GetNeutralVolume_dLambdai(iz); 
                    NI(0, -Rcx * dV_n_D/V_p * nions[Doffset + ir]); 
                    
                }
            }
        }
    }
