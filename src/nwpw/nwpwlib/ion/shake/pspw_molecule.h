/*
   pspw_molecule.h
   author - Eric Bylaska
*/
#ifndef _PSPW_MOLECULE_H_
#define _PSPW_MOLECULE_H_


/********************************/
/* the molecule list data structure */
/********************************/
#include 	"pspw_atom.h"

/* use a linked list for Molecular List and Atom List */
typedef	struct molecule_struct {
	struct molecule_struct	*next;
	
 	Atom_List_Type		atom;
        Bond_List_Type          bond;
        int			cyclic;
} *Molecule_List_Type;


extern int	pspw_molecule_add();
extern int 	pspw_molecule_size();
extern void	pspw_molecule_add_atom(int m, int a);
extern void	pspw_molecule_cyclic(int m, int cyclic);
extern void	pspw_molecule_init();

/*
extern void	pspw_molecule_read(char *filename);
extern void	pspw_molecule_data(int *m,
                           int *asize,
                           int *alist,
                           int *cyclic);
extern void	pspw_molecule_end();
*/

#endif
/* $Id$ */
