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
    const len_t NZ = this->ions->GetNZ();
    
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
        
        // Positive charge-exchange term
        if (Z == 1){ //Deuterium or Tritium
            if (Z0 == 1){
                for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    if(iz==iIon) //Skip if Deuterium/Tritium with itself
                        continue;
                    len_t Zi = ions->GetZ(iz); //Get Z for other ion
                    const len_t IonOffset = ions->GetIndex(iz,0); //Get index of neutral state of other ion
                    ADASRateInterpolator *ccd = adas->GetCCD(Zi); //Get cx-coeff. for the other ion
                    for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ //Loop over all charge states of other ion
                        real_t Rcx = ccd->Eval(Z0i, n_cold[ir], T_cold[ir]); //Evaluate cx-coeff. for the charge state
                        NI(-1, Rcx * V_n/V_p * nions[IonOffset+Z0i*Nr+ir]); //First argument in NI 0 because we want the neutral density for D/T (and we have Z0=1 here)
                    }
                }
            }
        }else if (Z0 < Z){  //Not Deuterium/Tritium
            for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                if(ions->GetZ(iz)!=1) //Don't add anything if the other ion is not D/T
                    continue;
                const len_t Doffset = ions->GetIndex(iz,0); //Get index of neutral state of D
                ADASRateInterpolator *ccd = adas->GetCCD(Z); //Get cx-coeff. for the ion that is not D/T (or should this be 1, 2+IsTritium(iz)?)
                real_t Rcx = ccd->Eval(Z0+1, n_cold[ir], T_cold[ir]); //Evaluate cx-coeff. for charge state 
                const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                if (Z0 == 0){
                    NI(+1, Rcx * V_n_D/V_n_tot * nions[Doffset + ir]); 
                }else{
                    NI(+1, Rcx * V_n_D/V_p * nions[Doffset + ir]);
                }
            }
        }
        
        // Negative charge-exchange term
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
                        NI(0, -Rcx * V_n/V_n_tot * nions[IonOffset+Z0i*Nr+ir]); //First argument is 0 since we want the neutral density for D/T (and we have Z0=0 here)
                    }
                }
            }
        } else if (Z0 > 1){  //Not Deuterium/Tritium. Z0>1 since this term not present if Z0=0
            for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                if(ions->GetZ(iz)!=1) //Don't add anything if the other ion is not D/T
                    continue;
                const len_t Doffset = ions->GetIndex(iz,0); //Get index of neutral state of D
                ADASRateInterpolator *ccd = adas->GetCCD(Z); //Get cx-coeff. for the ion that is not D/T (or should this be 1, 2+IsTritium(iz)?)
                real_t Rcx = ccd->Eval(Z0, n_cold[ir], T_cold[ir]); //Evaluate cx-coeff. for charge state 
                const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                NI(0, -Rcx * V_n_D/V_p * nions[Doffset + ir]); 
                
            }
        }
    }
