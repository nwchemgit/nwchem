/* hilbert.c -
$Id$
   Author - Eric Bylaska

   This file contains 2d hilbert mapping routines

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "typesf2c.h"
#include "olist.h"

#if defined(WIN32) && !defined(__MINGW32__)
#define hilbert2d_map_ HILBERT2D_MAP
#endif


#define bottom_left     0
#define bottom_right    1
#define top_left        2
#define top_right       3

#define right   0
#define left    1
#define up      2
#define down    3

int     hilbert2d(int i, int j, int level);
int     hilbert_dir(int i, int j, int level, int high, int* start);

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define pspsolve_ PSPSOLVE
#endif

void FATR hilbert2d_map_(sizex_ptr,sizey_ptr,map)
Integer     *sizex_ptr,*sizey_ptr;
Integer     map[];
{
    int i,j,size,sizex,sizey;
    int ii,jj,iii,jjj;
    int level,count;
    double dx,dy,x2,y2,dx2,dy2;
    OList_Type olist2;
    int *map2,*mapii,*mapjj;

   sizex = *sizex_ptr;
   sizey = *sizey_ptr;


   size = sizex;
   if (sizey>size) size = sizey;

   /* get the level of map */
   count = 1;
   level = 0;
   while (count < size)
   {
      ++level;
      count = count*2;
   }

   map2  = (int *) malloc(count*count*sizeof(int));
   mapii = (int *) malloc(count*count*sizeof(int));
   mapjj = (int *) malloc(count*count*sizeof(int));


   create_olist(&olist2,count*count);
   for (jj=0; jj<count; ++jj)
   for (ii=0; ii<count; ++ii)
   {

      map2[ii+jj*count] = hilbert2d(ii,jj,level);
      insert_olist(&olist2,map2[ii+jj*count]);
   }
   for (jj=0; jj<count; ++jj)
   for (ii=0; ii<count; ++ii)
      map2[ii+jj*count] = index_olist(&olist2,map2[ii+jj*count]);

   destroy_olist(&olist2);

   for (jj=0; jj<count; ++jj)
   for (ii=0; ii<count; ++ii)
   {
      mapii[map2[ii+jj*count]] = ii;
      mapjj[map2[ii+jj*count]] = jj;
   }


   dx2 = 1.0/((double) count);
   dy2 = 1.0/((double) count);
   dx =  1.0/((double) sizex);
   dy =  1.0/((double) sizey);

   for (j=0; j<sizex*sizey; ++j) map[j] = -9;

   iii = 0;
   for (jjj=0; jjj<count*count; ++jjj)
   {
      ii = mapii[jjj];
      jj = mapjj[jjj];

      x2 = dx2*(ii+0.5);
      y2 = dy2*(jj+0.5);
      i = rint((x2/dx) - 0.5);
      j = rint((y2/dy) - 0.5);

      if (map[i+j*sizex] == -9)
      {
        map[i+j*sizex] = iii;
        ++iii;
      }
   }
   free(mapii);
   free(mapjj);
   free(map2);

}


int     hilbert2d(i,j,level)
int     i,j;
int     level;
{
   int  start,direction;

   direction = hilbert_dir(i,j,level,level,&start);

   return start;
}


int     parent(i)
int     i;
{
   return(i/2);
}

int     corner(i,j)
int     i,j;
{
   return(2*(j%2) + (i%2));
}

int     hilbert_dir(i,j,level,high,start)
int     i,j;
int     level,high;
int     *start;
{
   int  direction,parent_direction,
        crnr,length,
        count;

   length = 1;
   for (count=0; count<(high-level); ++count)
      length = length*4;

   if (level == 0)
   {
      direction = right;
      *start    = 0;
   }
   else
   {
      parent_direction = hilbert_dir(parent(i),parent(j),
                                     level-1,high,
                                     start);
      crnr = corner(i,j);


      if (parent_direction == right)
      {
         if (crnr == bottom_left)
         {
            direction = up;
            *start    = *start + 0*length;
         }
         if (crnr == bottom_right)
         {
            direction = down;
            *start    = *start + 3*length;
         }
         if (crnr == top_left)
         {
            direction = right;
            *start    = *start + 1*length;
         }
         if (crnr == top_right)
         {
            direction = right;
            *start    = *start + 2*length;
         }
      }

      if (parent_direction == left)
      {
         if (crnr == bottom_left)
         {
            direction = left;
            *start    = *start + 2*length;
         }
         if (crnr == bottom_right)
         {
            direction = left;
            *start    = *start + 1*length;
         }
         if (crnr == top_left)
         {
            direction = up;
            *start    = *start + 3*length;
         }
         if (crnr == top_right)
         {
            direction = down;
            *start    = *start + 0*length;
         }
      }

      if (parent_direction == up)
      {
         if (crnr == bottom_left)
         {
            direction = right;
            *start    = *start + 0*length;
         }
         if (crnr == bottom_right)
         {
            direction = up;
            *start    = *start + 1*length;
         }
         if (crnr == top_left)
         {
            direction = left;
            *start    = *start + 3*length;
         }
         if (crnr == top_right)
         {
            direction = up;
            *start    = *start + 2*length;
         }
      }

      if (parent_direction == down)
      {
         if (crnr == bottom_left)
         {
            direction = down;
            *start    = *start + 2*length;
         }
         if (crnr == bottom_right)
         {
            direction = right;
            *start    = *start + 3*length;
         }
         if (crnr == top_left)
         {
            direction = down;
            *start    = *start + 1*length;
         }
         if (crnr == top_right)
         {
            direction = left;
            *start    = *start + 0*length;
         }
      }
   }

   return direction;
}
