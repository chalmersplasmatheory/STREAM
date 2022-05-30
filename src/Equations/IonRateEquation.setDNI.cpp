/**
 * Implementation of derivate w.r.t. ion densities in the charge-exchange
 * term. These derivates are NOT w.r.t. the density of the ion to which
 * the term is applied, but w.r.t. to the other ion density which is being
 * multiplied with.
 */
    
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); 
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);
    const len_t NZ = this->ions->GetNZ();
    const real_t *W_i = this->unknowns->GetUnknownData(id_Wi);
    const real_t *N_i = this->unknowns->GetUnknownData(id_Ni);
    const real_t ec = DREAM::Constants::ec;
    real_t Ti;

    if (this->includeChargeExchange) {
        //ADASRateInterpolator *ccdIon = GetCCD(iIon);
        ADASRateInterpolator *ccdIon;// = GetCCD(iIon);
        if(ions->IsTritium(iIon)){
    	    ccdIon = adas->GetCCD(1,3);
	} else if (ions->IsHydrogen(iIon)) {
	    ccdIon = adas->GetCCD(1,1);
        } else { 
           ccdIon = GetCCD(iIon);
        }
        for (len_t ir = 0; ir < Nr; ir++) {
            // Positive charge-exchange term
            real_t WA = this->unknowns->GetUnknownData(id_Wi)[iIon*Nr+ir];
            real_t NA = this->unknowns->GetUnknownData(id_Ni)[iIon*Nr+ir];
            real_t TA;
            if (NA <= 0) TA = 0;
            else TA = 2.0/3.0 * WA / (DREAM::Constants::ec*NA);
            if (Z == 1){ //Deuterium or Tritium
                const len_t DOffset = ions->GetIndex(iIon, 0);
                real_t nD1 = ions->GetIonDensity(ir, iIon, 1);
                
                real_t Rcx_ion = ccdIon->Eval(0, nD1, TA);
                real_t PartialnRcx_ion = ccdIon->Eval_deriv_n(0, nD1, TA);

                for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    if(iz==iIon) //Skip if Deuterium/Tritium with itself
                        continue;
                    len_t Zi = ions->GetZ(iz); //Get Z for other ion
                    const len_t IonOffset = ions->GetIndex(iz,0); //Get index of neutral state of other ion
                    //ADASRateInterpolator *ccd = GetCCD(iz); //Get cx-coeff. for the other ion
                    ADASRateInterpolator *ccd;
                    if(ions->IsTritium(iz)){
		        ccd = adas->GetCCD(1,3);
	            } else if (ions->IsHydrogen(iz)) {
		        ccd = adas->GetCCD(1,1);
	            } else { 
		        ccd = GetCCD(iz);
	            }
                    for(len_t Z0i=1; Z0i<Zi+1; Z0i++){ //Loop over all charge states of other ion
                        real_t ni = ions->GetIonDensity(ir, iz, Z0i);
                        real_t N_i_temp = N_i[iz*Nr+ir];
                        if (N_i_temp == 0)
                            Ti = 0;
                        else
                            Ti = 2.0/3.0*W_i[iz*Nr+ir]/(ec*N_i_temp);
                        real_t Rcx = ccd->Eval(Z0i-1, ni, Ti); //Evaluate cx-coeff. for the charge state
                        real_t PartialnRcx = ccd->Eval_deriv_n(Z0i-1, ni, Ti);

                        const real_t V_n_iz = this->volumes->GetNeutralVolume(iz);

                        if (Z0 == 0) {
                            // Apply to neutral deuterium (Z0=0)
                            NI_Z(IonOffset+Z0i, 0, -Rcx * V_n/V_n_tot - PartialnRcx * V_n/V_n_tot * nions[(IonOffset+Z0i)*Nr+ir]); //First argument is 0 since we want the neutral density for D/T (and we have Z0=0 here)

                            // D-T term (absent in DYON)
                            if (Zi == 1) {
                                NI_Z(IonOffset+Z0i, 1, Rcx_ion*V_n_iz/V_n_tot);
                                NI_Z(DOffset+1, 1, PartialnRcx_ion*V_n_iz/V_n_tot * nions[(IonOffset+0)*Nr+ir]);
                            }
                        } else if (Z0 == 1) {
                            // Apply to neutral deuterium (Z0-1 = 0)
                            NI_Z(IonOffset+Z0i, -1, Rcx * V_n/V_p + PartialnRcx * V_n/V_p * nions[(IonOffset+Z0i)*Nr+ir]); //First argument in NI 0 because we want the neutral density for D/T (and we have Z0=1 here)

                            // D-T term (absent in DYON)
                            if (Zi == 1) {
                                NI_Z(IonOffset+0, 0, -Rcx_ion * V_n_iz/V_p);
                                NI_Z(DOffset+1, 0, -PartialnRcx_ion * V_n_iz/V_p * nions[(IonOffset+0)*Nr+ir]);
                            }
                        }
                    }
                }
            }else if (Z0 < Z){  //Not Deuterium/Tritium
                real_t nZ0 = ions->GetIonDensity(ir, iIon, Z0);
                
                real_t Rcx_ion = ccdIon->Eval(Z0, nZ0, TA);
                for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    if(ions->GetZ(iz)!=1) //Don't add anything if the other ion is not D/T
                        continue;
                    const len_t Doffset = ions->GetIndex(iz,0); //Get index of neutral state of D

                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                    if (Z0 == 0){
                        NI_Z(Doffset, +1, Rcx_ion * V_n_D/V_n_tot); 
                    }else{
                        NI_Z(Doffset, +1, Rcx_ion * V_n_D/V_p);
                    }
                }
            }
            
            // Negative charge-exchange term
            if (Z != 1 && Z0 >= 1){  //Not Deuterium/Tritium. Z0>1 since this term not present if Z0=0
                real_t nZ0 = ions->GetIonDensity(ir, iIon, Z0);
                
                real_t Rcx_ion = ccdIon->Eval(Z0, nZ0, TA);
                for (len_t iz=0; iz<NZ; iz++){ //Loop over all other ion species
                    if(ions->GetZ(iz)!=1) //Don't add anything if the other ion is not D/T
                        continue;
                    const len_t Doffset = ions->GetIndex(iz,0); //Get index of neutral state of D
                    
                    const real_t V_n_D = this->volumes->GetNeutralVolume(iz); 
                    NI_Z(Doffset, 0, -Rcx_ion * V_n_D/V_p); 
                    
                }
            }
        }
    }
