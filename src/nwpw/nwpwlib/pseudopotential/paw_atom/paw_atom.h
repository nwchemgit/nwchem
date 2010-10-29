#ifndef	_PAW_ATOM_H_
#define _PAW_ATOM_H_

/*
   $Id$
*/



extern void   paw_init_atom(char *aname, char *infile);
extern void   paw_solve_atom();
extern void   paw_print_atom();

extern double paw_get_ion_energy();

extern void paw_save_atom();
extern void paw_plot_orbitals();
extern void paw_generate_paw_basis();

extern double paw_get_atom_kinetic_energy();
extern char* paw_get_atom_name();
extern char* paw_get_comment();
extern double paw_get_atom_total_energy();


extern void   paw_init_paw_atom(char *aname);
extern void   paw_solve_paw_atom();
extern void   paw_print_paw_atom();




#endif

