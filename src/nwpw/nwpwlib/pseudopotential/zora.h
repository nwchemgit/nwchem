#ifndef _ZORA_H_
#define _ZORA_H_
/* dirac.c */
extern void R_ZORA (int n, int l, int s2, double Z, const double *v,
                         int *mch, double *Eig, double *u, double *uprime);
extern void R_ZORA_Fixed_E (int n, int l, int s2, double Z, const double *v,
                                 int match, double E, double *u, double *uprime);
#endif
/* $Id$ */
