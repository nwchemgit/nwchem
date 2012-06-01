#ifndef	_ATOM_H_
#define _ATOM_H_
/* atom.h -
   author - Eric Bylaska and Patrick Nichols

*/
#include <stdio.h>

/* Solver type: Solve_Type */
#define	Schrodinger	-8001
#define Pauli		-8002
#define Dirac		-8003
#define ZORA            -8004
/* atom.c */
extern void init_Atom(char *filename);
extern void end_Atom();
extern void Thomas_Fermi(double Z, double Vtmp[]);
extern void solve_Atom(void);
extern void solve_Scattering_State_Atom(int nt, int lt, 
	double et, double rmax);
extern char *solver_Name_Atom(void);
extern int solver_Type_Atom(void);
extern void print_Atom(FILE *fp);
extern void set_Solver_Atom(int solver);
extern double E_Atom(void);
extern double eigenvalue_Atom(int i);
extern double *rho_Atom(void);
extern double *rho_core_Atom(void);
extern double *rho_valence_Atom(void);
extern double *Vall_Atom(void);
extern double *r_psi_Atom(int i);
extern double *r_psi_prime_Atom(int i);
extern int n_Atom(int i);
extern int l_Atom(int i);
extern int s_Atom(int i);
extern int lmax_Atom(void);
extern double fill_Atom(int i);
extern int Ncore_Atom(void);
extern int Nvalence_Atom(void);
extern double peak_Atom(int i);
extern int turning_point_Atom(int i);
extern double Zion_Atom(void);
extern double Amass_Atom(void);
extern int state_Atom(int nt, int lt);
extern int state_RelAtom(int nt, int lt, int st);
extern char *name_Atom(void);
extern char *spin_Name(int i);
extern int isRelativistic_Atom(void);
#endif
/* $Id$ */
