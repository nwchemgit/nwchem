#ifndef _PAW_LOG_GRID_H_
#define _PAW_LOG_GRID_H_
/*
   $Id$
*/



extern void     paw_init_LogGrid_from_file();
extern double   paw_r0_LogGrid();
extern double   *paw_alloc_LogGrid();
extern void	    paw_dealloc_LogGrid(double *grid);
extern double   *paw_r_LogGrid();
extern double   *paw_r2_LogGrid();
extern double   *paw_r3_LogGrid();
extern double   *paw_scratch_LogGrid();
extern int      paw_N_LogGrid();
extern double   paw_log_amesh_LogGrid();
extern double   paw_amesh_LogGrid();
extern double	  paw_Integrate_LogGrid();
extern double	  paw_Def_Integr();
extern void	    paw_Zero_LogGrid();
extern void	    paw_Copy_LogGrid();
extern void     paw_Copy_spin_LogGrid(double **gridnew, double **gridold);
extern double	  paw_Norm_LogGrid();
extern void	    paw_Derivative_LogGrid();
extern double	  paw_Def_Integr();
extern double*	paw_Indef_Integr();
extern double   paw_dot_product(double*, double*);
extern double   paw_dot_product1(int ,double*, double*);
extern int      paw_get_grid_index(double r);
extern int      paw_get_grid_index(double r);
extern void     paw_print_loggrid_information();
extern double paw_dot_product1(int n, double *f, double *g);

/* this one was missing before */
extern void  paw_end_LogGrid();

#endif


