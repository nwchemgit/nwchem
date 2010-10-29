/*
   $Id$
*/

#include        <stdlib.h>
#include        <stdio.h>
#include        <string.h>
#include        <math.h>
#include        "paw_my_constants.h"

#define MAXIT 400


/****************************************
 Function name	  : paw_bisection
 Description	    :
 Return type		  : double
 Argument         : *my_func)(double)
 Argument         : double x1
 Argument         : double x2
 Argument         : double eps
 Author     		  : Marat Valiev
 Date & Time		  : 9/30/98
****************************************/
double paw_bisection(double (*my_func)(double), double x1, double x2, double eps)
{

    int j;
    int root_found;
    double dx,f,fmid,xmid,root;

    root_found = False;

    fmid = my_func(x2);
    f    = my_func(x1);

    if (f*fmid >=0.0)
    {

        printf(" paw_bisection.c : root must be bracketed");
        exit(99);

    }

    if (f < 0.0)
    {
        root = x1;
        dx   = x2-x1;
    }
    else
    {
        root = x2;
        dx   = x1 - x2;
    }

    for (j=0;j<MAXIT;++j)
    {
        dx   = 0.5*dx;
        xmid = root + dx;
        fmid = my_func(xmid);

        if (fmid <= 0.0)
            root = xmid;

        if (fabs(dx) < eps || fmid == 0.0 )
        {
            root_found = True;
            break;
        }
    }

    if ( !(root_found))
    {
        printf("unable to find root within alotted number of iterations");
        exit(1);
    }

    return root;

}

