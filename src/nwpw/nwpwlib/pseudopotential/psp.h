#ifndef	_PSP_H_
#define _PSP_H_
/* psp.h -
   author - Eric Bylaska

*/


/* Solver type: Solve_Type */
#define	Hamann		-8401
#define	Troullier	-8402
#define Vanderbilt      -8403

extern void	init_Psp(char*);
extern void	end_Psp();
extern void	solve_Psp();
extern void	print_Psp();
extern double	E_Psp();
extern double	eigenvalue_Psp(int);
extern double	*V_Psp();
extern double	*rho_Psp();
extern double	*rho_semicore_Psp();
extern double	*drho_semicore_Psp();
extern double	r_semicore_Psp();
extern double	*r_psi_Psp(int );
extern int	n_Psp(int);
extern int	ns_Psp(int);
extern int	l_Psp(int);
extern int	lmax_Psp();
extern double	fill_Psp(int);
extern int	Nvalence_Psp();
extern double	peak_Psp(int);
extern double	rcut_Psp(int);
extern double	Zion_Psp();
extern int	state_Psp(int, int);
extern int      Vanderbilt_Psp();
extern double   *Vlocal_Psp();
extern double   *Beta_Psp();
extern double   *r_hard_psi_il_Psp();
extern double   *r_psi_il_Psp();
extern double   D0_Psp();
extern double   q_Psp();
extern double   rcut_il_Psp(int,int);
extern int      Vanderbilt_Psp();
extern int      NormConserving_Psp();
extern char     *comment_Psp();

extern int      kb_extra_Psp();
extern int      kb_expansion_Psp(int);
extern double   *r_psi_extra_Psp(int );


/* used for setting solver parameters */
extern void	set_Solver_Psp();

#endif
/* $Id$ */
