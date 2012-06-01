/*
 $Id$
*/
#ifndef _LOG_GRID_H_
#define _LOG_GRID_H_
/* LogGrid.h - 6/9/95
   author - Eric Bylaska

   This file contains the data structure for handeling numerics
on a logarithmic grid.  The grid is defined from 0 to 45.0,
with the grid points defined by:

        r(i) = (a**i)*r0

        with r0 = 0.00625/Z
             a  = 1.0247
             i = 0,1,...,N; N = log(7200.0*Z)/AL

*/


extern void     init_LogGrid();
extern void     end_LogGrid();
extern double  *alloc_LogGrid();
extern void	dealloc_LogGrid();
extern double  *r_LogGrid();
extern int      N_LogGrid();
extern double   log_amesh_LogGrid();
extern double   amesh_LogGrid();
extern double	Integrate_LogGrid();
extern double  Integrate_LogGrid_na_nb(int,int,double*);
extern double	Integrate2_LogGrid();
extern void	Zero_LogGrid();
extern void	Copy_LogGrid();
extern double	Norm_LogGrid();
extern void	Derivative_LogGrid();
extern void	Plot_LogGrid(char*,double*);
extern int  	index_r_LogGrid();

#ifdef WIN32
/* Microsoft C does not implement rint */
#define rint(x) floor(x)

#endif
#include <math.h>
#if defined(CRAY) &&!defined(__crayx1)
#include <fp.h>
#endif

#endif
