#ifndef _PAW_PRED_CORR_H_
#define _PAW_PRED_CORR_H_
/*
   $Id$
*/


/* paw_Pred_Corr.h - 6/9/95
   author      - Eric Bylaska

   This file contains the 5th order predictor-corrector
formulas for integrating inward and outward.
This file also contains 5th order derivatives.

*/
extern	double  paw_Predictor_In(int, double *, double *);
extern	double	paw_Predictor_Out(int, double *, double *);
extern double   paw_Predictor_Out_F(int, double *);

extern	double	paw_Corrector_In(int, double *, double *);
extern  double  paw_Corrector_In_F(int, double *);
extern	double	paw_Corrector_Out(int, double *, double *);

extern	double	paw_Derivative5_1(int, double *);
extern	double	paw_Derivative5_2(int, double *);
extern	double	paw_Derivative5_3(int, double *);
extern	double	paw_Derivative5_4(int, double *);
extern	double	paw_Derivative5_5(int, double *);

#endif


