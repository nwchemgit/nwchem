/*
 $Id$
*/
#ifndef _PRED_CORR_H_
#define _PRED_CORR_H_

/* Pred_Corr.h - 6/9/95
   author      - Eric Bylaska

   This file contains the 5th order predictor-corrector
formulas for integrating inward and outward.
This file also contains 5th order derivatives.

*/
extern	double  Predictor_In(int, double *, double *);
extern	double	Predictor_Out(int, double *, double *);

extern	double	Corrector_In(int, double *, double *);
extern  double Corrector_In_F(int, double *);
extern	double	Corrector_Out(int, double *, double *);

extern	double	Derivative5_1(int, double *);
extern	double	Derivative5_2(int, double *);
extern	double	Derivative5_3(int, double *);
extern	double	Derivative5_4(int, double *);
extern	double	Derivative5_5(int, double *);

extern	double	Derivative7_1(int, double *);
extern	double	Derivative7_2(int, double *);
extern	double	Derivative7_3(int, double *);
extern	double	Derivative7_4(int, double *);
extern	double	Derivative7_5(int, double *);
extern	double	Derivative7_6(int, double *);
extern	double	Derivative7_7(int, double *);
extern	double	Laplacian7_4(int, double *);

#endif
