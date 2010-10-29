/*
 $Id$
   Pred_Cor.c - 6/9/95
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

#include        "pred_cor.h"

#define one_over_24     4.16666666667e-2
#define one_over_60	1.66666666667e-2
#define	one_over_180	5.55555555556e-3
#define	one_over_720	1.38888888889e-3

/********************************
 *                              *
 *       Predictor_In           *
 *                              *
 ********************************/

/*
   Where F is 5 point predictor.

   y[i-1] <-- y[i] + F(f[i],f[i+1],...,f[i+3],f[i+4]);

   Entry - i,
           y,f;
   Exit  - returns the value of y[?] + F(....)
*/

double  Predictor_In(i,y,f)
int     i;
double  y[],
f[];
{

    double tmp;

    tmp = y[i]  - one_over_24*(  55.0*f[i]
                                 - 59.0*f[i+1]
                                 + 37.0*f[i+2]
                                 -  9.0*f[i+3]);
    /*
       tmp = y[i]  - one_over_720*( 1901.0*f[i]
                                  - 1387.0*f[i+1]
                                  +  327.0*f[i+2]
                                  -  637.0*f[i+3]
                                  +  251.0*f[i+4]);
    */

    return tmp;
} /* Predictor_In */




/********************************
 *                              *
 *       Predictor_Out          *
 *                              *
 ********************************/

/*
   Where F is 4 point predictor.

   y[i+1] <-- y[i] + F(f[i],f[i-1],...,f[i-3],f[i-4]);

   Entry - i,
           y,f;
   Exit  - returns the value of F(....)


*/

double  Predictor_Out(i,y,f)
int     i;
double  y[],
f[];
{

    double tmp;

    tmp =  y[i] + one_over_24*(  55.0*f[i]
                                 - 59.0*f[i-1]
                                 + 37.0*f[i-2]
                                 -  9.0*f[i-3]);
    /*
       tmp =  y[i] + one_over_720*( 1901.0*f[i]
                                  - 1387.0*f[i-1]
                                  +  327.0*f[i-2]
                                  -  637.0*f[i-3]
                                  +  251.0*f[i-4]);
    */


    return tmp;

} /* Predictor_Out */



/********************************
 *                              *
 *       Corrector_In           *
 *                              *
 ********************************/

/*
   Where F is 5 point Corrector.

   y[i-1] <-- y[i] + F(f[i-1],f[i],f[i+1],f[i+2],f[i+3]);

   Entry - i,
           y,f;
   Exit  - returns the value of y[i] + F(....)


*/

double  Corrector_In(i,y,f)
int     i;
double  y[],
f[];
{

    double tmp;

    tmp = y[i] - one_over_24*(   9.0*f[i-1]
                                 + 19.0*f[i]
                                 -  5.0*f[i+1]
                                 +  1.0*f[i+2]);
    /*
       tmp = y[i] - one_over_720*(  251.0*f[i-1]
                                 +  646.0*f[i]
                                 -  264.0*f[i+1]
                                 +  106.0*f[i+2]
                                 -   19.0*f[i+3]);
    */

    return tmp;

} /* Corrector_In */

/********************************
 *                              *
 *       Corrector_In_F         *
 *                              *
 ********************************/

/*
   Where F is 5 point Corrector.

   y[i-1] <-- y[i] + F(f[i-1],f[i],f[i+1],f[i+2],f[i+3]);

   Entry - i,
           f;
   Exit  - returns the value of  F(....)


*/

double  Corrector_In_F(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp = -one_over_24*(   9.0*f[i-1]
                           + 19.0*f[i]
                           -  5.0*f[i+1]
                           +  1.0*f[i+2]);
    /*
       tmp = -one_over_720*( 251.0*f[i-1]
                           + 646.0*f[i]
                           - 264.0*f[i+1]
                           + 106.0*f[i+2]
                           -  19.0*f[i+3]);
    */

    return tmp;

} /* Corrector_In_F */



/********************************
 *                              *
 *       Corrector_Out          *
 *                              *
 ********************************/

/*
   Where F is 5 point Corrector.

   y[i+1] <-- y[i] + F(f[i+1],f[i],f[i-1],f[i-2],f[i-3]);

   Entry - i,
           y,f;
   Exit  - returns the value of F(....)


*/

double  Corrector_Out(i,y,f)
int     i;
double  y[],
f[];
{

    double tmp;

    tmp = y[i] + one_over_24*(   9.0*f[i+1]
                                 + 19.0*f[i]
                                 -  5.0*f[i-1]
                                 +  1.0*f[i-2]);
    /*
       tmp = y[i] + one_over_720*( 251.0*f[i+1]
                                 + 646.0*f[i]
                                 - 264.0*f[i-1]
                                 + 106.0*f[i-2]
                                 -  19.0*f[i-3]);
    */

    return tmp;

} /* Corrector_Out */


/********************************
 *                              *
 *       Derivative7_1          *
 *                              *
 ********************************/

/*  This routine takes the 7 point derivative

    h*f'[i] = dx(f[i],f[i+1],...,f[i+6])


*/

double  Derivative7_1(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*( -147.0*f[i]
                         + 360.0*f[i+1]
                         - 450.0*f[i+2]
                         + 400.0*f[i+3]
                         - 225.0*f[i+4]
                         + 72.0*f[i+5]
                         - 10.0*f[i+6]);

    return tmp;
}

/********************************
 *                              *
 *       Derivative7_2          *
 *                              *
 ********************************/

/*  This routine takes the 7 point derivative

    h*f'[i] = dx(f[i-1],f[i],...,f[i+5])


*/

double  Derivative7_2(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*( -10.0*f[i-1]
                         -  77.0*f[i]
                         + 150.0*f[i+1]
                         - 100.0*f[i+2]
                         +  50.0*f[i+3]
                         -  15.0*f[i+4]
                         +   2.0*f[i+5]);

    return tmp;
}


/********************************
 *                              *
 *       Derivative7_3          *
 *                              *
 ********************************/

/*  This routine takes the 7 point derivative

    h*f'[i] = dx(f[i-2],f[i-1],...,f[i+4])


*/

double  Derivative7_3(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*(  +2.0*f[i-2]
                          -  24.0*f[i-1]
                          -  35.0*f[i]
                          +  80.0*f[i+1]
                          -  30.0*f[i+2]
                          +   8.0*f[i+3]
                          -   1.0*f[i+4]);

    return tmp;
}

/********************************
 *                              *
 *       Derivative7_4          *
 *                              *
 ********************************/

/*  This routine takes the 7 point derivative

    h*f'[i] = dx(f[i-3],f[i-2],...,f[i+2],f[i+3])


*/

double  Derivative7_4(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*(  -1.0*f[i-3]
                          +   9.0*f[i-2]
                          -  45.0*f[i-1]

                          +  45.0*f[i+1]
                          -   9.0*f[i+2]
                          +   1.0*f[i+3]);

    return tmp;
}


/********************************
 *                              *
 *       Derivative7_5          *
 *                              *
 ********************************/

/*  This routine takes the 7 point derivative

    h*f'[i] = dx(f[i+2],f[i+1],...,f[i-4])


*/

double  Derivative7_5(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*(  -2.0*f[i+2]
                          +  24.0*f[i+1]
                          +  35.0*f[i]
                          -  80.0*f[i-1]
                          +  30.0*f[i-2]
                          -   8.0*f[i-3]
                          +   1.0*f[i-4]);

    return tmp;
}

/********************************
 *                              *
 *       Derivative7_6          *
 *                              *
 ********************************/

/*  This routine takes the 7 point derivative

    h*f'[i] = dx(f[i+1],f[i],...,f[i-5])


*/

double  Derivative7_6(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*( +10.0*f[i+1]
                         +  77.0*f[i]
                         - 150.0*f[i-1]
                         + 100.0*f[i-2]
                         -  50.0*f[i-3]
                         +  15.0*f[i-4]
                         -   2.0*f[i-5]);
    return tmp;
}




/********************************
 *                              *
 *       Derivative7_7          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i],f[i-1],...,f[i-6])


*/

double  Derivative7_7(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_60*( +147.0*f[i]
                         - 360.0*f[i-1]
                         + 450.0*f[i-2]
                         - 400.0*f[i-3]
                         + 225.0*f[i-4]
                         -  72.0*f[i-5]
                         +  10.0*f[i-6]);
    return tmp;
}


/********************************
 *                              *
 *       Laplacian7_4           *
 *                              *
 ********************************/

/*  This routine takes the 7 point laplacian

    h*f'[i] = dx(f[i-3],f[i-2],...,f[i+2],f[i+3])


*/

double  Laplacian7_4(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_180*(  +2.0*f[i-3]
                           -  27.0*f[i-2]
                           + 270.0*f[i-1]
                           - 490.0*f[i]
                           + 270.0*f[i+1]
                           -  27.0*f[i+2]
                           +   2.0*f[i+3]);

    return tmp;
}



/********************************
 *                              *
 *       Derivative5_1          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i],f[i+1],...,f[i+4])


*/

double  Derivative5_1(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_24*( -50.0*f[i]
                         + 96.0*f[i+1]
                         - 72.0*f[i+2]
                         + 32.0*f[i+3]
                         -  6.0*f[i+4]);

    return tmp;

} /* Derivative5_1 */



/********************************
 *                              *
 *       Derivative5_2          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-1],f[i],...,f[i+3])


*/

double  Derivative5_2(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_24*(  -6.0*f[i-1]
                          - 20.0*f[i]
                          + 36.0*f[i+1]
                          - 12.0*f[i+2]
                          +  2.0*f[i+3]);

    return tmp;

} /* Derivative5_2 */








/********************************
 *                              *
 *       Derivative5_3          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-2],f[i-1],...,f[i+2])


*/

double  Derivative5_3(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_24*(   2.0*f[i-2]
                           - 16.0*f[i-1]

                           + 16.0*f[i+1]
                           -  2.0*f[i+2]);

    return tmp;

} /* Derivative5_3 */




/********************************
 *                              *
 *       Derivative5_4          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-3],f[i-2],...,f[i+1])


*/

double  Derivative5_4(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_24*(  -2.0*f[i-3]
                          + 12.0*f[i-2]
                          - 36.0*f[i-1]
                          + 20.0*f[i]
                          +  6.0*f[i+1]);

    return tmp;

} /* Derivative5_4 */




/********************************
 *                              *
 *       Derivative5_5          *
 *                              *
 ********************************/

/*  This routine takes the 5 point derivative

    h*f'[i] = dx(f[i-4],f[i-2],...,f[i])


*/

double  Derivative5_5(i,f)
int     i;
double  f[];
{

    double tmp;

    tmp =  one_over_24*(   6.0*f[i-4]
                           - 32.0*f[i-3]
                           + 72.0*f[i-2]
                           - 96.0*f[i-1]
                           + 50.0*f[i]);

    return tmp;

} /* Derivative5_5 */

