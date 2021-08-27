/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    //const real_t *T_cold = this->unknowns->GetUnknownData(id_T_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    const len_t NZ = this->ions->GetNZ();
    const real_t *W_i = this->unknowns->GetUnknownData(id_Wi);
    const real_t *N_i = this->unknowns->GetUnknownData(id_Ni);
    const real_t ec = DREAM::Constants::ec;  
    real_t Ti;
    
    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)
            if (Z0 == 1){
                NI(-1, Ion[Z0-1][ir] * n_cold[ir] * V_n/V_p, posIoniz);
            } else if (Z0 > 1){
                NI(-1, Ion[Z0-1][ir] * n_cold[ir], posIoniz);
            }

            // -I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
            if (Z0 == 0){
                NI(0, -Ion[Z0][ir] * n_cold[ir] * V_n/V_n_tot, negIoniz);
            }else{
                NI(0, -Ion[Z0][ir] * n_cold[ir], negIoniz);
            }
        }
        
        // R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)
        if (Z0 == 0){
            NI(+1, Rec[Z0+1][ir] * n_cold[ir] * V_p/V_n_tot, posRec);
        } else if (Z0 < Z){
            NI(+1, Rec[Z0+1][ir] * n_cold[ir], posRec);
        }

        // -R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
        // Does not contribute when Z0=0 since there is no recombination for neutrals
        if (Z0 > 0){                        
            NI(0, -Rec[Z0][ir] * n_cold[ir], negRec);
        }
        
        // Positive charge-exchange term
        if (this->includeChargeExchange) {
            if (Z == 1){ //Deuterium or Tritium
                //if (Z0 == 1){
                    for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    len_t Zi = ions->GetZ(iz); //Get Z for other ion
                        if(Z==Zi) //Skip if Deuterium/Tritium with itself or Tritium/Deuterium
                            continue;
                        const len_t IonOffset = ions->GetIndex(iz,0); //Get index of neutral state of other ion
                        ADASRateInterpolator *ccd = adas->GetCCD(Zi); //Get cx-coeff. for the other ion
                        for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ //Loop over all charge states of other ion
                            real_t ni = ions->GetIonDensity(ir, iz, Z0i);
                            real_t N_i_temp = N_i[iz*Nr+ir];
                            if (N_i_temp == 0)
                                Ti = 0;
                            else
                                Ti = 2.0/3.0*W_i[iz*Nr+ir]/(ec*N_i_temp);
                            real_t Rcx = ccd->Eval(Z0i-1, ni, Ti); //Evaluate cx-coeff. for the charge state

                            if (Z0 == 0)
                                // Apply to neutral deuterium (Z0=0)
                                NI(0, -Rcx * V_n/V_n_tot * nions[(IonOffset+Z0i)*Nr+ir], negCX); //First argument is 0 since we want the neutral density for D/T (and we have Z0=0 here)
                            else if (Z0 == 1)
                                // Apply to neutral deuterium (Z0-1 = 0)
                                NI(-1, Rcx * V_n/V_p * nions[(IonOffset+Z0i)*Nr+ir], posCX); //First argument in NI 0 because we want the neutral density for D/T (and we have Z0=1 here)
                        }
                    }
                //}
            }else if (Z0 < Z){  //Not Deuterium/Tritium
                for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    if(ions->GetZ(iz)!=1) //Don't add anything if the other ion is not D/T
                        continue;
                    const len_t Doffset = ions->GetIndex(iz,0); //Get index of neutral state of D
                    ADASRateInterpolator *ccd = adas->GetCCD(Z); //Get cx-coeff. for the ion that is not D/T (or should this be 1, 2+IsTritium(iz)?)
                    real_t ni = ions->GetIonDensity(ir, iIon, Z0);
                    real_t N_i_temp = N_i[iIon*Nr+ir];
                    if (N_i_temp == 0)
                        Ti = 0;
                    else
                        Ti = 2.0/3.0*W_i[iIon*Nr+ir]/(ec*N_i_temp);
                        
                    real_t Rcx = ccd->Eval(Z0+1-1, ni, Ti); //Evaluate cx-coeff. for charge state 
                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                    if (Z0 == 0){
                        NI(+1, Rcx * V_n_D/V_n_tot * nions[Doffset*Nr + ir], posCX); 
                    }else{
                        NI(+1, Rcx * V_n_D/V_p * nions[Doffset*Nr + ir], posCX);
                    }
                }
            }
            
            // Negative charge-exchange term
            if (Z != 1 && Z0 >= 1){  //Not Deuterium/Tritium. Z0>1 since this term not present if Z0=0
                for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    if(ions->GetZ(iz)!=1) //Don't add anything if the other ion is not D/T
                        continue;
                    const len_t Doffset = ions->GetIndex(iz,0); //Get index of neutral state of D
                    ADASRateInterpolator *ccd = adas->GetCCD(Z); //Get cx-coeff. for the ion that is not D/T (or should this be 1, 2+IsTritium(iz)?)
                    real_t ni = ions->GetIonDensity(ir, iIon, Z0);
                    real_t N_i_temp = N_i[iIon*Nr+ir];
                    if (N_i_temp == 0)
                        Ti = 0;
                    else
                        Ti = 2.0/3.0*W_i[iIon*Nr+ir]/(ec*N_i_temp);
                    real_t Rcx = ccd->Eval(Z0-1, ni, Ti); //Evaluate cx-coeff. for charge state 
                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                    NI(0, -Rcx * V_n_D/V_p * nions[Doffset*Nr + ir], negCX); 
                    
                }
            }
        }
    }
