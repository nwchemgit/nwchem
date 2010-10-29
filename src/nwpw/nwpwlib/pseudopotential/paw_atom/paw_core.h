#ifndef _PAW_CORE_H_
#define _PAW_CORE_H_
/*
   $Id$
*/


extern void  paw_init_core();
extern double* paw_get_pointer_ps_core_density();
extern double paw_generate_core_kin_energy();
extern double* paw_get_pointer_core_density();
extern double paw_get_core_charge();
extern double paw_get_ps_core_charge();
extern void  paw_set_core();
extern double paw_get_core_kinetic_energy();
extern void  paw_print_core_information(FILE *fp);

#endif

