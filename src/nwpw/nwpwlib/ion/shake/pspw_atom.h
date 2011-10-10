/*
   pspw_atom.h
   author - Eric Bylaska
*/
#ifndef _PSPW_ATOM_H_
#define _PSPW_ATOM_H_

/********************************/
/* the atom list data structure */
/********************************/

/* use a linked list for Molecular List and Atom List */
typedef	struct atom_struct {
	struct atom_struct	*next;
	int			a;
} *Atom_List_Type;



extern void	pspw_atom_init(Atom_List_Type *atom);
extern void	pspw_atom_add(Atom_List_Type *atom, int a);
extern int 	pspw_atom_size(Atom_List_Type atom);
extern int 	pspw_atom(Atom_List_Type atom, int indx);
extern void	pspw_atom_end(Atom_List_Type atom);

#endif
/* $Id$ */
