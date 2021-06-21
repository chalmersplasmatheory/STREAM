/**
 * Sets temperature jacobian of IonRateEquation
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const real_t V_p = this->volumes->GetPlasmaVolume();
    const real_t V_n = this->volumes->GetNeutralVolume(iIon); //Is the argument iIon the right one here?
    const real_t V_n_tot = this->volumes->GetTotalNeutralVolume(iIon);

    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // I_i^(j-1) n_cold * n_i^(j-1) * Vhat_i^(j-1)/V_i^(j)             
            if (Z0 == 1){
                NI(-1,); //TODO: PartialTIon[Z0-1][ir] * n_cold[ir] if lambda_i>a and for lambda_i<=a PartialTIon[Z0-1][ir] * n_cold[ir] * 1/V_p * dV_n/dT_cold. where dV_n/dT_cold is equal to 4pi^2*R*kappa*(a-lambda_i)*dlambda_i/dT_cold, and dlambda_i/dT_cold=-v_i/(n_e (I_i^(0))²) * PartialTIon[0][ir] (i forgot the triangularity here so it's even longer.) Maybe the derivative dV_n/dT_cold is better calculated in the PlasmaVolume class?
            } else if (Z0 > 1){
                NI(-1, PartialTIon[Z0-1][ir] * n_cold[ir]);
            }
            // -I_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
            if (Z0 == 0){
                NI(0, ); //TODO: d/dT_cold(-Ion[Z0][ir] * n_cold[ir] * V_n/V_n_tot)
            }else{
                NI(0, - PartialTIon[Z0][ir] * n_cold[ir] );
            }           
        }
        // R_i^(j+1) n_cold * n_i^(j+1) * Vhat_i^(j+1)/V_i^(j)
        if (Z0 == 0){
            NI(+1, ); //TODO: d/dT_e(Rec[Z0+1][ir] * n_cold[ir] * V_p/V_n_tot)
        }else if (Z0 < Z) {
            NI(+1, PartialTRec[Z0+1][ir] * n_cold[ir]);
        }

        // -R_i^(j) n_cold * n_i^(j) * Vhat_i^(j)/V_i^(j)
        if (Z0 > 0){                        
            NI(0, -PartialTRec[Z0][ir] * n_cold[ir]);
        }
    }
