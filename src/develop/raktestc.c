/*
 $Id$
 */

#include "typesf2c.h"
#include <stdio.h>
#if defined(SOLARIS)
#define WORKS_FOR 1
#endif

#if defined(CRAY) || defined(WIN32)
#define is_this_val_okay_ IS_THIS_VAL_OKAY
#endif

/* */

#if defined(WORKS_FOR)
#include <ieeefp.h>
#endif

int FATR is_this_val_okay_(double *value)
{
#if defined(WORKS_FOR)
/*
 -1 unknown
  0 zero
  1 nan
  2 infinity
  3 nonzero
 */
    double val;
    fpclass_t fp_res;
    int ret_val;

    val = *value;
    fp_res = fpclass (val);

    if (fp_res == FP_NZERO || fp_res == FP_PZERO)
	{
	    ret_val = 0;
	}
    else if (fp_res == FP_SNAN || fp_res == FP_QNAN) 
	{
	    ret_val = 1;
	}
    else if (fp_res == FP_NINF || fp_res == FP_PINF)
	{
	    ret_val = 2;
	}
    else if (fp_res == FP_NDENORM || fp_res == FP_PDENORM || fp_res == FP_NNORM || fp_res == FP_PNORM)
	{
	    ret_val = 3;
	}
    else
	{
	    ret_val = -1;
	}
    return ret_val;
#else
    return (int) 3;
#endif
}
