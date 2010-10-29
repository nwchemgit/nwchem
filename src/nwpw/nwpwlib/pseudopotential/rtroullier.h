/*
 $Id$
*/
#ifndef _TROULLIER_H_
#define _TROULLIER_H_
/* troullier.h -
   Author - Patrick Nichols

*/

extern void Suggested_Param_RelTroullier(int *num_states_psp, int *n_psp, 
int *l_psp, int *s_psp, double *e_psp, double *fill_psp, double *rcut_psp);

extern void solve_RelTroullier(int num_psp, int *n_psp, int *l_psp, 
int *s_psp, double *e_psp, double *fill_psp, double *rcut_psp, double **r_psi_psp, 
double **r_psi_prime_psp, double *rho_psp, double *rho_semicore, double **V_psp, 
double *eall_psp, double *eh_psp, double *ph_psp, double *ex_psp, double *px_psp, 
double *ec_psp, double *pc_psp);

#endif
