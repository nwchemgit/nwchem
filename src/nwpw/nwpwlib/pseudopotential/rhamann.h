#ifndef _RHAMANN_H_
#define _RHAMANN_H_
/* rhamann.c */
extern void Suggested_Param_RelHamann (int *num_states_psp,
				       int n_psp[], int l_psp[], int s_psp[],
				       double e_psp[], double fill_psp[],
				       double rcut_psp[]);

extern void solve_RelHamann (int num_psp, int *n_psp, int *l_psp,
			     int *s_psp, double *e_psp, double *fill_psp,
			     double *rcut_psp, double **r_psi_psp,
			     double **r_psi_prime_psp, double *rho_psp,
			     double *rho_semicore, double **V_psp,
			     double *eall_psp, double *eh_psp, double *ph_psp,
			     double *ex_psp, double *px_psp, double *ec_psp,
			     double *pc_psp);

#endif
/* $Id$ */
