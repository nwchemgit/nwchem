/*
   pspw_bond.c
   author - Eric Bylaska
*/

#include	<stdlib.h>
#include	<stdio.h>
#include	"pspw_bond.h"


#define	NIL	((void *) 0)



/********************************
 *				*
 *         pspw_bond_init	*
 *				*
 ********************************/
void	pspw_bond_init(Bond_List_Type *bond)
{
   *bond = NIL;

} /* pspw_bond_init */

/********************************
 *				*
 *         pspw_bond_add	*
 *				*
 ********************************/
void	pspw_bond_add(Bond_List_Type *bond, double a)
{
   Bond_List_Type prev,cur,node;

   node       = (Bond_List_Type) malloc(sizeof(struct bond_struct));
   node->a    = a;
   node->next = NIL;

   /* first bond on list */
   if (*bond == NIL)
   {
      *bond = node;
   }
   /* add to end of list */
   else
   {
      cur  = *bond;
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
 *         pspw_bond_size	*
 *				*
 ********************************/
int 	pspw_bond_size(Bond_List_Type bond)
{
   Bond_List_Type cur;
   int            count;

   cur   = bond;
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
 *         pspw_bond		*
 *				*
 ********************************/
double 	pspw_bond(Bond_List_Type bond, int indx)
{
    Bond_List_Type	cur;
    int 		i;

    if (pspw_bond_size(bond) < indx) 
    {
       printf("pspw_bond ERROR\n");
       exit(99);
    }

    cur = bond;
    for (i=1; i<(indx); ++i)
       cur = cur->next;

    return cur->a;
}

/********************************
 *				*
 *         pspw_bond_end	*
 *				*
 ********************************/
void	pspw_bond_end(Bond_List_Type bond)
{
    Bond_List_Type cur;

    while (bond != NIL)
    {
       cur = bond;
       bond = bond->next;
       free(cur);
    }
}
    
/* $Id$ */
