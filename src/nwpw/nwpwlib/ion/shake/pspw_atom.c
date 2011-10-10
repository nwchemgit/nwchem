/*
   pspw_atom.c
   author - Eric Bylaska
*/

#include	<stdlib.h>
#include	<stdio.h>
#include	"pspw_atom.h"


#define	NIL	((void *) 0)



/********************************
 *				*
 *         pspw_atom_init	*
 *				*
 ********************************/
void	pspw_atom_init(Atom_List_Type *atom)
{
   *atom = NIL;

} /* pspw_atom_init */

/********************************
 *				*
 *         pspw_atom_add	*
 *				*
 ********************************/
void	pspw_atom_add(Atom_List_Type *atom, int a)
{
   Atom_List_Type prev,cur,node;

   node       = (Atom_List_Type) malloc(sizeof(struct atom_struct));
   node->a    = a;
   node->next = NIL;

   /* first atom on list */
   if (*atom == NIL)
   {
      *atom = node;
   }
   /* add to end of list */
   else
   {
      cur  = *atom;
      while (cur != NIL)
      {
         prev = cur;
         cur  = cur->next;
      }
      prev->next = node;
   }

}

/********************************
 *				*
 *         pspw_atom_size	*
 *				*
 ********************************/
int 	pspw_atom_size(Atom_List_Type atom)
{
   Atom_List_Type cur;
   int            count;

   cur   = atom;
   count = 0;
   while (cur != NIL)
   {
      cur = cur->next;
      ++count;
   }
   return count;
}

/********************************
 *				*
 *         pspw_atom		*
 *				*
 ********************************/
int 	pspw_atom(Atom_List_Type atom, int indx)
{
    Atom_List_Type	cur;
    int 		i;

    if (pspw_atom_size(atom) < indx) 
    {
       printf("pspw_atom ERROR\n");
       exit(99);
    }

    cur = atom;
    for (i=1; i<(indx); ++i)
       cur = cur->next;

    return cur->a;
}

/********************************
 *				*
 *         pspw_atom_end	*
 *				*
 ********************************/
void	pspw_atom_end(Atom_List_Type atom)
{
    Atom_List_Type cur;

    while (atom != NIL)
    {
       cur = atom;
       atom = atom->next;
       free(cur);
    }
}
    
/* $Id$ */
