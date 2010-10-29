#ifndef _PAW_HARTREE_H_
#define _PAW_HARTREE_H_
/*
   $Id$
*/


extern void paw_init_hartree();
extern void paw_generate_hartree_pot(double *n);
extern double paw_get_hartree_energy(double *n);
extern double* paw_get_hartree_pot();
#endif


