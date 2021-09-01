/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    //const real_t *T_cold = this->unknowns->GetUnknownData(id_T_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    const real_t dV_n = this->volumes->GetNeutralVolume_dLambdai(iIon); 
    const real_t dV_n_tot = this->volumes->GetTotalNeutralVolume_dLambdai(iIon);
    const len_t NZ = this->ions->GetNZ();
    const real_t *W_i = this->unknowns->GetUnknownData(id_Wi);
    const real_t *N_i = this->unknowns->GetUnknownData(id_Ni); 
    const real_t ec = DREAM::Constants::ec;  
    real_t Ti;

    if (V_n > 0) {
        for (len_t ir = 0; ir < Nr; ir++) {
            if(setIonization){
                // I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)
                if (Z0 == 1){
                    NI(-1, Ion[Z0-1][ir] * n_cold[ir] * dV_n/V_p);
                }

                // -I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
                if (Z0 == 0){
                    NI(0, -Ion[Z0][ir] * n_cold[ir] * (dV_n/V_n_tot - V_n/(V_n_tot*V_n_tot)*dV_n_tot)); 
                }
            }
            
            // R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)
            if (Z0 == 0){
                NI(+1, -Rec[Z0+1][ir] * n_cold[ir] * V_p/(V_n_tot*V_n_tot)*dV_n_tot);
            }

            // d/dlambda_i(Positive charge-exchange term)
            if (this->includeChargeExchange) {
                if (Z == 1){
                    if (Z0 == 1){
                        ADASRateInterpolator *ccdIon = GetCCD(iIon);

                        const len_t DOffset = ions->GetIndex(iIon, 0);
                        real_t nD1 = ions->GetIonDensity(ir, iIon, 1);
                        real_t WD  = this->unknowns->GetUnknownData(id_Wi)[iIon*Nr+ir];
                        real_t ND  = this->unknowns->GetUnknownData(id_Ni)[iIon*Nr+ir];
                        real_t TD;
                        if (ND <= 0) TD = 0;
                        else TD = 2.0/3.0 * WD / (DREAM::Constants::ec*ND);

                        real_t Rcx_ion = ccdIon->Eval(0, nD1, TD);

                        for (len_t iz=0; iz<NZ; iz++){ 
                            if(iz==iIon) 
                                continue;
                            len_t Zi = ions->GetZ(iz); 
                            const len_t IonOffset = ions->GetIndex(iz,0); 
                            ADASRateInterpolator *ccd = GetCCD(iz); 
                            for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ 
                                real_t ni = ions->GetIonDensity(ir, iz, Z0i);
                                real_t N_i_temp = N_i[iz*Nr+ir];
                                if (N_i_temp == 0)
                                    Ti = 0;
                                else
                                    Ti = 2.0/3.0*W_i[iz*Nr+ir]/(ec*N_i_temp);
                                real_t Rcx = ccd->Eval(Z0i-1, ni, Ti); 

                                NI_Z(iIon, -1, Rcx * dV_n/V_p * nions[(IonOffset+Z0i)*Nr+ir]); 

                                // D-T term (absent in DYON)
                                if (Zi == 1)
                                    NI_Z(DOffset+Z0, 0, Rcx_ion * dV_n/V_p * nions[(IonOffset+0)*Nr+ir]);
                            }
                        }
                    }
                }else if (Z0 < Z){  
                    ADASRateInterpolator *ccd = GetCCD(iIon); 
                    real_t ni = ions->GetIonDensity(ir, iIon, Z0);
                    real_t N_i_temp = N_i[iIon*Nr+ir];
                    if (N_i_temp == 0)
                        Ti = 0;
                    else
                        Ti = 2.0/3.0*W_i[iIon*Nr+ir]/(ec*N_i_temp);
                    real_t Rcx = ccd->Eval(Z0, ni, Ti); 
                    for (len_t iz=0; iz<NZ; iz++){ 
                        if(ions->GetZ(iz)!=1) 
                            continue;
                        const len_t Doffset = ions->GetIndex(iz,0); 
                        const real_t V_n_D = this->volumes->GetNeutralVolume(iz);
                        const real_t dV_n_D = this->volumes->GetNeutralVolume_dLambdai(iz);
                        if (Z0 == 0){
                            NI_Z(iz, +1, Rcx * (dV_n_D/V_n_tot - V_n_D * dV_n_tot /(V_n_tot * V_n_tot)) * nions[Doffset*Nr + ir]); 
                        }else{
                            NI_Z(iz, +1, Rcx * dV_n_D/V_p * nions[Doffset*Nr + ir]);
                        }
                    }
                }
                
                // d/dlambda_i(Negative charge-exchange term)
                if (Z == 1){
                    if(Z0 == 0){
                        ADASRateInterpolator *ccdIon = GetCCD(iIon);

                        const len_t DOffset = ions->GetIndex(iIon, 0);
                        real_t nD1 = ions->GetIonDensity(ir, iIon, 1);
                        real_t WD  = this->unknowns->GetUnknownData(id_Wi)[iIon*Nr+ir];
                        real_t ND  = this->unknowns->GetUnknownData(id_Ni)[iIon*Nr+ir];
                        real_t TD;
                        if (ND <= 0) TD = 0;
                        else TD = 2.0/3.0 * WD / (DREAM::Constants::ec*ND);

                        real_t Rcx_ion = ccdIon->Eval(0, nD1, TD);

                        for (len_t iz=0; iz<NZ; iz++){
                            if(iz==iIon)
                                continue;
                            len_t Zi = ions->GetZ(iz); 
                            const len_t IonOffset = ions->GetIndex(iz,0);
                            ADASRateInterpolator *ccd = GetCCD(iz);
                            for(len_t Z0i=1; Z0i<Zi+1; Z0i++){
                                real_t ni = ions->GetIonDensity(ir, iz, Z0i);
                                real_t N_i_temp = N_i[iz*Nr+ir];
                                if (N_i_temp == 0)
                                    Ti = 0;
                                else
                                    Ti = 2.0/3.0*W_i[iz*Nr+ir]/(ec*N_i_temp);

                                const real_t V_n_iz = this->volumes->GetNeutralVolume(iz);
                                const real_t dV_n_iz = this->volumes->GetNeutralVolume_dLambdai(iz);

                                real_t Rcx = ccd->Eval(Z0i-1, ni, Ti);
                                NI_Z(iIon, 0, -Rcx * (dV_n/V_n_tot - V_n * dV_n_tot/(V_n_tot*V_n_tot)) * nions[(IonOffset+Z0i)*Nr+ir]); 

                                // D-T CX term (absent in DYON)
                                if (Zi == 1) {
                                    NI_Z(DOffset, +1, -Rcx_ion * V_n_iz * dV_n_tot/(V_n_tot*V_n_tot) * nions[(IonOffset+Z0i)*Nr+ir]);
                                    NI_Z(iz+Z0i, +1, Rcx_ion * dV_n_iz/V_n_tot * nions[(IonOffset+Z0i)*Nr+ir]);
                                }
                            }
                        }
                    }
                } else if (Z0 >= 1){  
                    ADASRateInterpolator *ccd = GetCCD(iIon); 
                    real_t ni = ions->GetIonDensity(ir, iIon, Z0);
                    real_t N_i_temp = N_i[iIon*Nr+ir];
                    if (N_i_temp == 0)
                        Ti = 0;
                    else
                        Ti = 2.0/3.0*W_i[iIon*Nr+ir]/(ec*N_i_temp);
                    real_t Rcx = ccd->Eval(Z0-1, ni, Ti); 
                    for (len_t iz=0; iz<NZ; iz++){ 
                        if(ions->GetZ(iz)!=1) 
                            continue;
                        const len_t Doffset = ions->GetIndex(iz,0); 
                        const real_t dV_n_D = this->volumes->GetNeutralVolume_dLambdai(iz); 
                        NI_Z(iz, 0, -Rcx * dV_n_D/V_p * nions[Doffset*Nr + ir]); 
                        
                    }
                }
            }
        }
    }
