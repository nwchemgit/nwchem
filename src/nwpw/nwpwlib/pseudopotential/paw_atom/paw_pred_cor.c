/*
   $Id$
*/


/* paw_Pred_Cor.c - 6/9/95
   author      - Eric Bylaska

   This file contains the 5th order predictor-corrector
formulas for integrating inward and outward.
This file also contains 5th order derivatives.
These routines also ask for y, so that y[i] can
be added to F(...), this is done because some
orders have correctors of the form

Corrector_In
    y[i-1] <- y[i+5] + F(...)
*/

#include        "paw_pred_cor.h"
#define one_over_24     4.16666666667e-2

/********************************
 *                              *
 *       paw_Predictor_In           *
 *                              *
 ********************************/

/*
   Where F is 4 point predictor.

   y[i-1] <-- y[i] + F(f[i],f[i+1],...,f[i+3]);

   Entry - i,
           y,f;
   Exit  - returns the value of y[?] + F(....)
*/

double  paw_Predictor_In(int i,double  y[], double f[])
{

    double tmp;

    tmp = y[i]  - one_over_24*(  55.0*f[i]
                                 - 59.0*f[i+1]
                                 + 37.0*f[i+2]
                                 -  9.0*f[i+3]);

    return tmp;
} /* paw_Predictor_In */




/********************************
 *                              *
 *       paw_Predictor_Out          *
 *                              *
 ********************************/

/*
   Where F is 4 point predictor.

   y[i+1] <-- y[i] + F(f[i],f[i-1],...,f[i-3]);

   Entry - i,
           y,f;
   Exit  - returns the value of F(....)


*/

double  paw_Predictor_Out(int i,double y[],double f[])
{

    double tmp;

    tmp =  y[i] + one_over_24*(  55.0*f[i]
                                 - 59.0*f[i-1]
                                 + 37.0*f[i-2]
                                 -  9.0*f[i-3]);

    return tmp;

} /* paw_Predictor_Out */


double  paw_Predictor_Out_F(int i,double f[])
{

    double tmp;

    tmp =  one_over_24*(  55.0*f[i]
                          - 59.0*f[i-1]
                          + 37.0*f[i-2]
                          -  9.0*f[i-3]);

    return tmp;

} /* paw_Predictor_Out_F */

/********************************
 *                              *
 *       paw_Corrector_In           *
 *                              *
 ********************************/

/*
   Where F is 4 point Corrector.

   y[i-1] <-- y[i] + F(f[i-1],f[i],f[i+1],f[i+2]);

   Entry - i,
           y,f;
   Exit  - returns the value of y[i] + F(....)


*/

double  paw_Corrector_In(int i,double y[],double f[])
{

    double tmp;

    tmp = y[i] - one_over_24*(   9.0*f[i-1]
                                 + 19.0*f[i]
                                 -  5.0*f[i+1]
                                 +  1.0*f[i+2]);

    return tmp;

} /* paw_Corrector_In */

/********************************
 *                              *
 *       paw_Corrector_In_F         *
 *                              *
 ********************************/

/*
   Where F is 4 point Corrector.

   y[i-1] <-- y[i] + F(f[i-1],f[i],f[i+1],f[i+2]);

   Entry - i,
           f;
   Exit  - returns the value of  F(....)


*/

double  paw_Corrector_In_F(int i,double f[])
{

    double tmp;

    tmp = -one_over_24*(   9.0*f[i-1]
                           + 19.0*f[i]
                           -  5.0*f[i+1]
                           +  1.0*f[i+2]);

    return tmp;

} /* paw_Corrector_In_F */



/********************************
 *                              *
 *       paw_Corrector_Out          *
 *                              *
 ********************************/

/*
   Where F is 4 point Corrector.

   y[i+1] <-- y[i] + F(f[i+1],f[i],f[i-1],f[i-2]);

   Entry - i,
           y,f;
   Exit  - returns the value of F(....)


*/

double  paw_Corrector_Out(int i,double y[],double f[])
{

    double tmp;

    tmp = y[i] + one_over_24*(   9.0*f[i+1]
                                 + 19.0*f[i]
                                 -  5.0*f[i-1]
                                 +  1.0*f[i-2]);

    return tmp;

} /* paw_Corrector_Out */




/********************************
 *                              *
 *       paw_Derivative5_1          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i],f[i+1],...,f[i+4])


*/

double  paw_Derivative5_1(int i, double f[])
{

    double tmp;

    tmp =  one_over_24*( -50.0*f[i]
                         + 96.0*f[i+1]
                         - 72.0*f[i+2]
                         + 32.0*f[i+3]
                         -  6.0*f[i+4]);

    return tmp;

} /* paw_Derivative5_1 */



/********************************
 *                              *
 *       paw_Derivative5_2          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-1],f[i],...,f[i+3])


*/

double  paw_Derivative5_2(int i,double f[])
{
    double tmp;

    tmp =  one_over_24*(  -6.0*f[i-1]
                          - 20.0*f[i]
                          + 36.0*f[i+1]
                          - 12.0*f[i+2]
                          +  2.0*f[i+3]);

    return tmp;

} /* paw_Derivative5_2 */








/********************************
 *                              *
 *       paw_Derivative5_3          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-2],f[i-1],...,f[i+2])


*/

double  paw_Derivative5_3( int i,double f[])
{

    double tmp;

    tmp =  one_over_24*(   2.0*f[i-2]
                           - 16.0*f[i-1]

                           + 16.0*f[i+1]
                           -  2.0*f[i+2]);

    return tmp;

} /* paw_Derivative5_3 */




/********************************
 *                              *
 *       paw_Derivative5_4          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-3],f[i-2],...,f[i+1])


*/

double  paw_Derivative5_4(int i,double f[])
{

    double tmp;

    tmp =  one_over_24*(  -2.0*f[i-3]
                          + 12.0*f[i-2]
                          - 36.0*f[i-1]
                          + 20.0*f[i]
                          +  6.0*f[i+1]);

    return tmp;

} /* paw_Derivative5_4 */




/********************************
 *                              *
 *       paw_Derivative5_5          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-4],f[i-2],...,f[i])


*/

double  paw_Derivative5_5(int i,double  f[])
{

    double tmp;

    tmp =  one_over_24*(   6.0*f[i-4]
                           - 32.0*f[i-3]
                           + 72.0*f[i-2]
                           - 96.0*f[i-1]
                           + 50.0*f[i]);

    return tmp;

} /* paw_Derivative5_5 */



