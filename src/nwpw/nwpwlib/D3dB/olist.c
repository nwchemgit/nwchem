/* olist.c -
$Id$
   Author - Eric Bylaska

   This file contains an ordered list data structure.

*/

#include	<stdio.h>
#include	<stdlib.h>
#include	"olist.h"

void	create_olist(olist,size)
OList_Type	*olist;
int		size;
{
   int i;
   olist->list = (int *) malloc(size*sizeof(int));

   olist->max_index = 0;
   for (i=0; i<size; ++i)
      olist->list[i] = 0;
}

void	insert_olist(olist,item)
OList_Type	*olist;
int		item;
{
   int	index,j,
        max_index;

   olist->max_index = olist->max_index + 1;
   max_index = olist->max_index;

   index = 0;
   while (((olist->list[index]) < item) && (index < (max_index-1)))
      ++index;

   for (j=(max_index-1); j>index; --j)
      olist->list[j] = olist->list[j-1];
   olist->list[index] = item;

}

int	index_olist(olist,item)
OList_Type	*olist;
int		item;
{
   int index;

   index = 0;
   while ( (olist->list[index]) < item)
      ++index;

   return index;
}

void	destroy_olist(olist)
OList_Type	*olist;
{
   free(olist->list);
}

void	print_olist(olist)
OList_Type	*olist;
{
   int i;

   printf("%d: ",olist->max_index);
   for (i=0; i<(olist->max_index); ++i)
      printf("%d ", olist->list[i]);
   printf("\n");
 }
