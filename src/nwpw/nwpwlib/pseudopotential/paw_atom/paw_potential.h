#ifndef _PAW_POTENTIAL_H_
#define _PAW_POTENTIAL_H_
/*
   $Id$
*/


extern void   paw_Thomas_Fermi(double Z, double Vtmp[]);
extern void   paw_init_potential();
extern void   paw_find_kohn_sham_potential(double *rho,double *V_ks);
extern void    paw_set_kohn_sham_potential(double *rho);
extern double* paw_get_kohn_sham_potential();


extern void paw_init_paw_potential(int,double,double,double*,double*);
extern double* paw_get_paw_potential(int);
extern void paw_update_paw_potential(int*, int , double , double ,double*);
extern double* paw_get_ref_pot();
extern void paw_generate_pseudopot();
extern double* paw_get_pointer_pseudopotential();
extern double paw_get_potential_matching_radius();
extern void  paw_print_paw_potential_information(FILE *fp);
extern void paw_print_paw_potential_to_file(char* atom_name);


#endif

