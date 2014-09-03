/*
   pspw_molecule.c
   author - Eric Bylaska
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typesf2c.h"

#include	"pspw_atom.h"
#include        "pspw_bond.h"
#include	"pspw_molecule.h"

#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define pspw_molecule_read_ PSPW_MOLECULE_READ
#define pspw_molecule_data_ PSPW_MOLECULE_DATA
#define pspw_molecule_end_ 	PSPW_MOLECULE_END
#define pspw_molecule_msize_ PSPW_MOLECULE_MSIZE
#endif




#define	NIL	((void *) 0)

static Molecule_List_Type	molecule;


/********************************
 *				*
 *         pspw_molecule_init	*
 *				*
 ********************************/
void	pspw_molecule_init()
{
   molecule = NIL;

} /* pspw_molecule_init */

/********************************
 *				*
 *         pspw_molecule_add	*
 *				*
 ********************************/
int	pspw_molecule_add()
{
   Molecule_List_Type prev,cur,node;
   int	m;

   node       = (Molecule_List_Type) malloc(sizeof(struct molecule_struct));
   node->next = NIL;
   pspw_atom_init(&(node->atom));
   pspw_bond_init(&(node->bond));
   node->cyclic = 0;

   /* first atom on list */
   if (molecule == NIL)
   {
      molecule = node;
   }
   /* add to end of list */
   else
   {
      cur  = molecule;
      while (cur != NIL)
      {
         prev = cur;
         cur  = cur->next;
      }
      prev->next = node;
   }

   m = pspw_molecule_size();

   return m;
}

/********************************
 *				*
 *         pspw_molecule_size	*
 *				*
 ********************************/
int 	pspw_molecule_size()
{
   Molecule_List_Type cur;
   int            count;

   cur   = molecule;
   count = 0;
   while (cur != NIL)
   {
      cur = cur->next;
      ++count;
   }
   return count;
}

/************************************************
 *						*
 *         pspw_molecule_add_atom		*
 *						*
 ************************************************/
void	pspw_molecule_add_atom(int m, int a)
{
    Molecule_List_Type	cur;
    int 		i;

    if (pspw_molecule_size() < m) 
    {
       printf("pspw_molecule_add_atom ERROR %d %d\n",pspw_molecule_size(),m);
       exit(99);
    }

    cur = molecule;
    for (i=1; i<(m); ++i)
       cur = cur->next;

    pspw_atom_add(&(cur->atom), a);


}

/************************************************
 *						*
 *         pspw_molecule_add_bond		*
 *						*
 ************************************************/
void	pspw_molecule_add_bond(int m, double a)
{
    Molecule_List_Type	cur;
    int 		i;

    if (pspw_molecule_size() < m) 
    {
       printf("pspw_molecule_add_bond ERROR %d %d\n",pspw_molecule_size(),m);
       exit(99);
    }
    cur = molecule;
    for (i=1; i<(m); ++i)
       cur = cur->next;

    pspw_bond_add(&(cur->bond), a);


}


/************************************************
 *						*
 *         pspw_molecule_cyclic			*
 *						*
 ************************************************/
void	pspw_molecule_cyclic(int m, int cyclic)
{
    Molecule_List_Type	cur;
    int 		i;

    if (pspw_molecule_size() < m) 
    {
       printf("ERROR\n");
       exit(99);
    }

    cur = molecule;
    for (i=1; i<(m); ++i)
       cur = cur->next;

    cur->cyclic = cyclic;
}

/************************************************
 *												*
 *         pspw_molecule_data  					*
 *												*
 ************************************************/
void FATR pspw_molecule_data_
		(Integer *m,
		 Integer  *aasize,
		 Integer  *alist,
	 DoublePrecision  *blist,
		 Integer  *cyclic)
{
    Molecule_List_Type	cur;
    int 		i,asize,bsize;

    if (pspw_molecule_size() < *m) 
    {
       printf("ERROR\n");
       exit(99);
    }

    cur = molecule;
    for (i=1; i<(*m); ++i)
       cur = cur->next;

    *cyclic = cur->cyclic;
    asize   = pspw_atom_size(cur->atom);
    *aasize = asize;
    for (i=0; i<asize; ++i)
       alist[i] = pspw_atom(cur->atom, i+1);

    bsize = asize-1;
    if (*cyclic) bsize = bsize + 1; 
    if (*cyclic==2) bsize = 1;
    for (i=0; i<bsize; ++i)
       blist[i] = pspw_bond(cur->bond, i+1);
}
   

/********************************
 *				*
 *      pspw_molecule_msize	*
 *				*
 ********************************/
void FATR pspw_molecule_msize_(Integer *msize)
{
   *msize = pspw_molecule_size();
    
}


/********************************
 *				*
 *         pspw_molecule_end	*
 *				*
 ********************************/
void FATR pspw_molecule_end_()
{
    Molecule_List_Type cur;

    while (molecule != NIL)
    {
       cur = molecule;
       molecule = molecule->next;
       pspw_atom_end(cur->atom);
       pspw_bond_end(cur->bond);
       free(cur);
    }
}
    
/********************************
 *				*
 *      pspw_molecule_read	*
 *				*
 ********************************/

void FATR pspw_molecule_read_
#if defined(USE_FCD)
( _fcd fcd_filename, Integer *n1)
{
 const char *filename = _fcdtocp(fcd_filename);
			  
#else
(char *filename, Integer *n1)
{
#endif


   FILE *fp;
   int value;
   int j,m,msize;
   int i,asize,bsize;
   double x;
   char cyclic;

   char *file = (char *) malloc(*n1+1);
   (void) strncpy(file, filename, *n1);
   file[*n1] = '\0';

   pspw_molecule_init();
   fp = fopen(file,"r+");
   value = 1;
   msize = 0;
   while (value)
   {
      value = fscanf(fp,"%d",&j);
      if (value != EOF) 
      {
         ++msize;
         asize = 1;
         m     = pspw_molecule_add();
         while (value==1)
         {
            pspw_molecule_add_atom(m,j);
            value = fscanf(fp,"%d",&j);
            ++asize;
         }
         --asize;
         fscanf(fp,"%c",&cyclic);
         if      (cyclic == 'd') pspw_molecule_cyclic(m,2);
         else if (cyclic == 'c') pspw_molecule_cyclic(m,1);
         else               pspw_molecule_cyclic(m,0);
         bsize = asize - 1;
         if (cyclic == 'c') bsize = bsize+1;
         if (cyclic == 'd') bsize = 1;
         for (i=0; i<bsize; ++i)
	 {
            if (!fscanf(fp,"%lf",&x))
	    {
               printf("Error reading bond lengths/contraints\n");
               exit(88);
            }
            pspw_molecule_add_bond(m,x);
	 }
         value = 1;
      }
      else
      {
         value = 0;
      }
   }
   fclose(fp);
   free(file);
}
/* $Id$ */
