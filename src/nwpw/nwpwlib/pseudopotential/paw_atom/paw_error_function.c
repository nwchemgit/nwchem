/*
   $Id$
*/

/************************************
  REVISION LOG ENTRY
  Revision By: Marat Valiev
  Revised on 1/8/99 1:03:48 PM
  Comments: ...
 ************************************/

/*************************************************

	Name    : paw_error_function.c

  Purpose : Erf related subroutines



  Created  : Marat Valiev, 7/22/98

*************************************************/

#include        <stdlib.h>
#include        <stdio.h>
#include        <string.h>
#include        <math.h>



#define  B1  0.0705230784
#define  B2  0.0422820123
#define  B3  0.0092705272
#define  B4  0.0001520143
#define  B5  0.0002765672
#define  B6  0.0000430638




/****************************************
 Function name	  : paw_my_erf
 Description	    : to calculate error function for
                    double precision scalar argument.
 Return type		  : double
 Argument         : double x
 Author     		  : Marat Valiev
 Algorithm        : original source unknown,
                    replicated from
                    Ryoichi Kawai's Cray Code
 Date & Time		  : 1/8/99 1:04:47 PM
****************************************/
double  paw_my_erf(double x)

{
    double f;

    f = (1.0+x*(B1+x*(B2+x*(B3+x*(B4+x*(B5+x*B6))))));

    f = 1.0/pow(f,16.0);

    f = 1.0 - f;

    return f;
}


