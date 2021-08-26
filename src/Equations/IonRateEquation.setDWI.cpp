/**
 * Implementation of derivate w.r.t. ion heat in the charge-exchange term.
 */    
    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    //const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
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
    // Positive charge-exchange term
        if (this->includeChargeExchange) {
            if (Z == 1){ //Deuterium or Tritium
                //if (Z0 == 1){
                    for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                        if(iz==iIon) //Skip if Deuterium/Tritium with itself
                            continue;
                        len_t Zi = ions->GetZ(iz); //Get Z for other ion
                        const len_t IonOffset = ions->GetIndex(iz,0); //Get index of neutral state of other ion
                        ADASRateInterpolator *ccd = adas->GetCCD(Zi); //Get cx-coeff. for the other ion
                        for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ //Loop over all charge states of other ion
                            real_t ni = ions->GetIonDensity(ir, iz, Z0i);
                            real_t N_i_temp = N_i[iz*Nr+ir];
                            if (N_i_temp == 0)
                                continue;

                            Ti = 2.0/3.0*W_i[iz*Nr+ir]/(ec*N_i_temp);
                            real_t PartialTRcx = ccd->Eval_deriv_T(Z0i-1, ni, Ti); //Evaluate cx-coeff. for the charge state

                            if (Z0 == 0)
                                // Apply to neutral deuterium (Z0=0)
                                NI(0, -PartialTRcx * 2.0/(3.0*ec*N_i[iz*Nr+ir]) * V_n/V_n_tot * nions[(IonOffset+Z0i)*Nr+ir]); //First argument is 0 since we want the neutral density for D/T (and we have Z0=0 here)
                            else if (Z0 == 1)
                                // Apply to neutral deuterium (Z0-1 = 0)
                                NI(-1, PartialTRcx * 2.0/(3.0*ec*N_i[iz*Nr+ir]) * V_n/V_p * nions[(IonOffset+Z0i)*Nr+ir]); //First argument in NI 0 because we want the neutral density for D/T (and we have Z0=1 here)
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
                        continue;

                    Ti = 2.0/3.0*W_i[iIon*Nr+ir]/(ec*N_i_temp);
                    real_t PartialTRcx = ccd->Eval_deriv_T(Z0+1-1, ni, Ti); //Evaluate cx-coeff. for charge state 
                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                    if (Z0 == 0){
                        NI(+1, PartialTRcx * 2.0/(3.0*ec*N_i[iz*Nr+ir]) * V_n_D/V_n_tot * nions[Doffset*Nr + ir]); 
                    }else{
                        NI(+1, PartialTRcx * 2.0/(3.0*ec*N_i[iz*Nr+ir]) * V_n_D/V_p * nions[Doffset*Nr + ir]);
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
                        continue;

                    Ti = 2.0/3.0*W_i[iIon*Nr+ir]/(ec*N_i_temp);
                    real_t PartialTRcx = ccd->Eval_deriv_T(Z0-1, ni, Ti); //Evaluate cx-coeff. for charge state 
                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                    NI(0, -PartialTRcx * 2.0/(3.0*ec*N_i[iz*Nr+ir]) * V_n_D/V_p * nions[Doffset*Nr + ir]); 
                    
                }
            }
        }
    }
