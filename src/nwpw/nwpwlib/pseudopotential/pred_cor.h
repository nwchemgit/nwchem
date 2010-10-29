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
extern	double  Predictor_In();
extern	double	Predictor_Out();

extern	double	Corrector_In();
extern  double Corrector_In_F();
extern	double	Corrector_Out();

extern	double	Derivative5_1();
extern	double	Derivative5_2();
extern	double	Derivative5_3();
extern	double	Derivative5_4();
extern	double	Derivative5_5();

extern	double	Derivative7_1();
extern	double	Derivative7_2();
extern	double	Derivative7_3();
extern	double	Derivative7_4();
extern	double	Derivative7_5();
extern	double	Derivative7_6();
extern	double	Derivative7_7();
extern	double	Laplacian7_4();

#endif
