#ifndef	_ATOM_H_
#define _ATOM_H_
/* atom.h -
   author - Eric Bylaska

*/


/* Solver type: Solve_Type */
#define	Schrodinger	-8001
#define Pauli		-8002
#define Dirac		-8003

extern void	init_Atom(char*);
extern void	solve_Atom();
extern void	solve_Scattering_State_Atom();
extern void	print_Atom();
extern double	E_Atom();
extern double	eigenvalue_Atom(int);
extern double	*Vall_Atom();
extern double	*rho_Atom();
extern double   *rho_valence_Atom();
extern double	*rho_core_Atom();
extern double	*r_psi_Atom(int );
extern double	*r_psi_prime_Atom(int );
extern int	n_Atom(int);
extern int	l_Atom(int);
extern int	lmax_Atom();
extern double	fill_Atom(int);
extern int	Ncore_Atom();
extern int	Nvalence_Atom();
extern double	peak_Atom(int);
extern int	turning_point_Atom(int);
extern double	Zion_Atom();
extern double	Amass_Atom();
extern int	state_Atom(int, int);
extern char	*name_Atom();

/* used for setting solver parameters */
extern void	set_Solver_Atom();

#endif
