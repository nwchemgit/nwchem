/*
   pspw_bond.h
   author - Eric Bylaska
*/
#ifndef _PSPW_BOND_H_
#define _PSPW_BOND_H_

/********************************/
/* the bond list data structure */
/********************************/

/* use a linked list for Molecular List and Atom List */
typedef	struct bond_struct {
	struct bond_struct	*next;
	double			a;
} *Bond_List_Type;


extern void		pspw_bond_init(Bond_List_Type *bond);
extern void		pspw_bond_add(Bond_List_Type *bond, double a);
extern int 		pspw_bond_size(Bond_List_Type bond);
extern double	pspw_bond(Bond_List_Type bond, int indx);
extern void		pspw_bond_end(Bond_List_Type bond);

#endif
/* $Id$ */
