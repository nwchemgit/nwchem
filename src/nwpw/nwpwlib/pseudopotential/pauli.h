/*
 $Id$
*/
#ifndef _PAULI_H_
#define _PAULI_H_
/* pauli.h - 6/9/95
   author     - Eric Bylaska

   This file contains routines for integrating the radial
   Pauli equation.

*/

extern void   R_Pauli(int, int, double ,double *, int *, double *, double *, double *);
extern void   R_Pauli_Fixed_E(int, int, double ,double *, int , double , double *, double *);

#endif
