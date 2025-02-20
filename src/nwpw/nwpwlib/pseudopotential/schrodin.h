/*
 $Id$
*/
#ifndef _SCHRODINGER_H_
#define _SCHRODINGER_H_
/* Schrodinger.h - 6/9/95
   author     - Eric Bylaska

   This file contains routines for integrating the radial
   Schodinger equation.

*/

extern void   R_Schrodinger(int, int, double *, int *, double *, double *, double *);
extern void   R_Schrodinger_Fixed_E(int, int, double *, int , double , double *, double *);
extern void   R_Schrodinger_Fixed_Logderiv(int, int, double *, int , double, double *, double *, double *);

#endif
