#ifndef	_PAW_ORBITALS_H_
#define _PAW_ORBITALS_H_
/*
   $Id$
*/


/* Solver type: Solve_Type */
#define Schrodinger     -8001
#define Pauli           -8002
#define Dirac           -8003
#define ZORA            -8004


extern void   paw_init_orbitals_from_file(FILE *fp);

extern void   paw_solve_orbital(int i,double *V);

extern int    paw_bound_state_test(int l, double *v);

extern void   paw_print_orbital_information(FILE *fp);

extern void   paw_solve_unoccupied_orbitals();

extern void   paw_solve_occupied_orbitals();

extern void   paw_solve_scattering_orbitals();

extern void   paw_print_orbitals_to_file(char* output);

extern int paw_get_Nscat();


extern double* paw_get_psi_prime(int);
extern double** paw_get_pointer_psi_prime_array();

extern double* paw_get_psi(int);
extern double** paw_get_pointer_psi_array();

extern double paw_get_fill(int);
extern double* paw_get_pointer_fill_array();

extern int paw_get_l(int);
extern int* paw_get_pointer_l_array();

extern double paw_get_e(int);

extern int paw_get_Ntotal();

extern int paw_get_orbital_index(int prin_n, int orb_l);

extern char* paw_orbital_type_name(int i);

extern int paw_get_orb_type(int i);

extern void paw_generate_density(double* rho);

extern double* paw_get_pointer_density();

#endif

