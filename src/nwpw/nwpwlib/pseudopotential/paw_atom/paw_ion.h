#ifndef _PAW_ION_H_
#define _PAW_ION_H_
/*
   $Id: paw_ion.h,v 1.2 2004-10-14 22:05:03 bylaska Exp $
*/


extern void   paw_init_ion(double Z);
extern double paw_get_ion_energy(double *dn);
extern double* paw_get_ion_pot();
extern double paw_get_ion_charge();

#endif


