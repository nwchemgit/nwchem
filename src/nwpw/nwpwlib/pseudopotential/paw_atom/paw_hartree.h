#ifndef _PAW_HARTREE_H_
#define _PAW_HARTREE_H_
/*
   $Id: paw_hartree.h,v 1.2 2004-10-14 22:05:03 bylaska Exp $
*/


extern void paw_init_hartree();
extern void paw_generate_hartree_pot(double *n);
extern double paw_get_hartree_energy(double *n);
extern double* paw_get_hartree_pot();
#endif


