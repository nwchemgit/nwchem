/*
 $Id: spline.h,v 1.1 2001-08-30 16:58:37 bylaska Exp $
*/
#ifndef _SPLINE_H_
#define _SPLINE_H_
/* spline.h -
    Taken from Numerical recipies, with slight modifications as
suggested by hamman's code.
*/

extern void	init_Linear(char*);
extern	int	nrl_Linear();
extern	double	drl_Linear();
extern void	Log_to_Linear();
extern void	normalize_Linear(double*);

#endif

