#ifndef _RPSP_H_
#define _RPSP_H_
#include <stdio.h>

/* rpsp.c */
extern void init_RelPsp(char *filename);
extern void solve_RelPsp(void);
extern char *solver_Name_RelPsp(void);
extern void print_RelPsp(FILE *fp);
extern void set_Solver_RelPsp(int solver);
extern int Vanderbilt_RelPsp(void);
extern int NormConserving_RelPsp(void);
extern double E_RelPsp(void);
extern double eigenvalue_RelPsp(int i);
extern double *rho_RelPsp(void);
extern double *rho_semicore_RelPsp(void);
extern double *drho_semicore_RelPsp(void);
extern double r_semicore_RelPsp(void);
extern double *Beta_RelPsp(int i, int l);
extern double *r_psi_il_RelPsp(int i, int l);
extern double *r_hard_psi_il_RelPsp(int i, int l);
extern int ns_RelPsp(int l);
extern double D0_RelPsp(int i, int j, int l);
extern double q_RelPsp(int i, int j, int l);
extern double *Vlocal_RelPsp(void);
extern double *V_RelPsp(int i);
extern double *r_psi_RelPsp(int i);
extern int n_RelPsp(int i);
extern int l_RelPsp(int i);
extern int lmax_RelPsp(void);
extern double fill_RelPsp(int i);
extern int Nvalence_RelPsp(void);
extern double peak_RelPsp(int i);
extern double rcut_RelPsp(int i);
extern double rcut_il_RelPsp(int i, int l);
extern double Zion_RelPsp(void);
extern int state_RelPsp(int nt, int lt, int st);
extern char *comment_RelPsp(void);
#endif
/* $Id$ */
