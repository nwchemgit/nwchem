/*
 $Id: spline.h,v 1.2 2004-05-24 13:43:19 bylaska Exp $
*/
#ifndef _SPLINE_H_
#define _SPLINE_H_
/* spline.h -
    Taken from Numerical recipies, with slight modifications as
suggested by hamman's code.
*/

extern void	init_Linear(char*);
extern void	end_Linear();
extern	int	nrl_Linear();
extern	double	drl_Linear();
extern void	Log_to_Linear();
extern void	normalize_Linear(double*);

#endif

