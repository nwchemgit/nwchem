#ifndef _DIRAC_H_
#define _DIRAC_H_
/* dirac.c */
extern void R_Dirac (int n, int l, int s2, double Z, const double *v,
                         int *mch, double *Eig, double *u, double *uprime);
extern void R_Dirac_Fixed_E (int n, int l, int s2, double Z, const double *v,
                                 int match, double E, double *u, double *uprime);
extern void R_Dirac_FixedLogDeriv(int n, int l, int s2, double Z, const double *v, 
	int match, double u_logderiv, double *Eig, double *u, double *uprime);

#endif
/* $Id$ */
