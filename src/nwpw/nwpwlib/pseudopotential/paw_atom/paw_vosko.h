#ifndef _PAW_VOSKO_H_
#define _PAW_VOSKO_H_
/*
   $Id: paw_vosko.h,v 1.2 2004-10-14 22:05:03 bylaska Exp $
*/


extern void paw_init_vosko();
extern void paw_generate_corr_pot(double **rho, double **Vc);
extern double paw_get_correlation_energy(double **rho);
extern void paw_generate_corr_pot_LDA(double *rho);
extern double* paw_get_corr_pot_LDA();
extern double paw_get_correlation_energy_LDA(double *rho);

#endif


