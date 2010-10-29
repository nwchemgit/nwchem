#ifndef _PAW_SCHRODINGER_H_
#define _PAW_SCHRODINGER_H_

/*
   $Id$
*/

/* Schrodinger.h - 6/9/95
   author     - Eric Bylaska

   This file contains routines for integrating the radial
   Schodinger equation.

*/

extern int   paw_R_Schrodinger(int n,
                                   int l,
                                   double* v,
                                   double *Eig,
                                   double* u,
                                   double* uprime);

extern  int paw_R_Schrodinger_Fixed_E(
        int l,
        double *v,
        int match,
        double E,
        double *u,
        double *uprime
    );
extern int paw_R_Schrodinger_Fixed_Logderiv(
        int n,
        int l,
        double *v,
        int match,
        double u_logderiv,
        double *Eig,
        double *u,
        double *uprime
    );


extern  void paw_R_Schrodinger_Fixed_E1(
        int l,
        double *v,
        double *f,
        int match,
        double E,
        double *u,
        double *uprime
    );

extern int paw_R_Schrodinger_Fixed_Logderiv1(
        int n,
        int l,
        double *v,
        int match,
        double u_logderiv,
        double *Eig,
        double *u,
        double *uprime
    );


#endif



