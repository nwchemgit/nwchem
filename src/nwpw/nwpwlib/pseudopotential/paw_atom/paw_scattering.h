#ifndef	_PAW_SCATTERING_H_
#define _PAW_SCATTERING_H_
/*
   $Id: paw_scattering.h,v 1.3 2004-10-14 22:05:03 bylaska Exp $
*/


extern void paw_init_paw_scattering_set();
extern void paw_init_paw_scattering();
extern void paw_end_paw_scattering();

extern void paw_solve_paw_scattering(int l, double r, double e, double* psi,double* psi_prime);

#endif
