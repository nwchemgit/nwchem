/*
   pspw_molecule.c
   author - Eric Bylaska
*/

#include	<stdlib.h>
#include	<stdio.h>
#include	"pspw_molecule.h"




main()
{
   int m,m1,msize;
   int i,asize,indx[50],cyclic;

   pspw_molecule_read("MOLECULE");

   msize = pspw_molecule_size();
   printf("msize: %d\n\n",msize);
   for (m=1; m<=msize; ++m)
   {
       pspw_molecule_data(&m,&asize,indx,&cyclic);
       printf("asize: %d\n",asize);
       for (i=0; i<asize; ++i)
          printf("%d ",indx[i]);
       if (cyclic) printf("c\n");
       else        printf("l\n");
   }
   pspw_molecule_end();

}
/* $Id$ */
