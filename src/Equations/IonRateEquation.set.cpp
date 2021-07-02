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
        
        // R_ik,cx^(j+1) * n_i^(j) * n_k^(0) * Vhat_i^(j+1)/V_i^(j) 
        if (Z == 1 && ions->IsTritium(ir)){ //ok? Z and not Z0, right?
            ADASRateInterpolator *ccd = adas->GetCCD(Z, 3); //Is this the best place to get it? 
            real_t Rcx = ccd->Eval(Z+1, n_cold[ir], T_cold[ir]);
            NI(+1, Rcx * V_p/V_n_tot); //Is the first argument in NI +1 here? How to get the n_k^(0) that we need to multiply with. unknowns->GetUnknownData(id_n_i)[ir] maybe, but then we need to make sure that it does not have Z=1 before we sum over all other ions?
        }else if (Z == 1){
            ADASRateInterpolator *ccd = adas->GetCCD(Z);
            real_t Rcx = ccd->Eval(Z+1, n_cold[ir], T_cold[ir]); //This does not exist here right, can we have Z+1, when Z is maximal ionization?
            NI(+1, Rcx * V_p/V_n_tot); //Need to multiply with n_k^(0) with k not D/T
        } else {
            ADASRateInterpolator *ccd = adas->GetCCD(Z0);
            real_t Rcx = ccd->Eval(Z0+1, n_cold[ir], T_cold[ir]);
            NI(+1, Rcx * V_p/V_n_tot); //Need to multiply with n_k^(0) with k equal to D
        }
        
        // -R_ik,cx^(j) * n_i^(j) * n_k^(0) * Vhat_i^(j)/V_i^(j)
        if (Z == 1 && ions->IsTritium(ir)){
            ADASRateInterpolator *ccd = adas->GetCCD(Z, 3);
            real_t Rcx = ccd->Eval(Z, n_cold[ir], T_cold[ir]);
            NI(0, -Rcx ); //Is the first argument in NI 0 here? Multiply with n_k^(0) for k not D/T
        } else if (Z == 1){
            ADASRateInterpolator *ccd = adas->GetCCD(Z);
            real_t Rcx = ccd->Eval(Z, n_cold[ir], T_cold[ir]);
            NI(0, -Rcx ); //Multiply with n_k^(0) for k not D/T
        } else{
            ADASRateInterpolator *ccd = adas->GetCCD(Z0);
            real_t Rcx = ccd->Eval(Z0, n_cold[ir], T_cold[ir]);
            NI(0, -Rcx );//Multiply with n_k^(0) for k equal to D/T
        }
    }
