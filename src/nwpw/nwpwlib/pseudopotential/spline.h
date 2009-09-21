/*
 $Id: spline.h,v 1.3 2005-03-07 20:50:48 bylaska Exp $
*/
#ifndef _SPLINE_H_
#define _SPLINE_H_
/* spline.h -
*/

extern void	init_Linear(char*);
extern void	end_Linear();
extern	int	nrl_Linear();
extern	double	drl_Linear();
extern void	Log_to_Linear();
extern void	Log_to_Linear_zero();
extern void	normalize_Linear(double*);

#endif

