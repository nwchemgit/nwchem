/*
 $Id$
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
extern void	normalize_Linear2(double*,double*);

#endif

