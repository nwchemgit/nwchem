#ifndef _PAW_DIRAC_EXCHANGE_H_
#define _PAW_DIRAC_EXCHANGE_H_
/*
   $Id$
*/


extern void paw_init_dirac_exchange();
extern double paw_get_exchange_energy(double **rho);
extern void paw_generate_exchange_pot_LDA(double *rho);
extern double paw_get_exchange_energy_LDA(double *rho);
extern double* paw_get_exchange_potential();

#endif

