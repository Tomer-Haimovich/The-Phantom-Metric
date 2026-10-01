// PHANTOM METRIC: BACKGROUND DENSITY (REPLACING CDM)
double rho_phantom = pba->Omega0_cdm * pow(pba->H0, 2) / (a * a * a);
pvecback[pba->index_bg_rho_cdm] = rho_phantom;