/*
 * hcurve.c
 *
 *  Created on: Dec 13, 2010
 *      Author: bylaska
 *
 *      Implements the H-Curve 2d mapping by Niedermeier, Reinhardt and Sanders.  This program
 *      is based on the Java applet of Reinhard.
 *      (http://www2-fs.informatik.uni-tuebingen.de/~reinhard/hcurve.html)
 *
 *      R. Niedermeier , K. Reinhardt and P. Sanders . Towards optimal locality in
 *      mesh-indexings. Proceedings of the 11th International Symposium on Fundamentals
 *      of Computation Theory, number 1279 in LNCS, pages 364--375, Krakow, Poland,
 *      September 1997. Springer.
 *      (http://www2-fs.informatik.uni-tuebingen.de/~reinhard/fct97.pdf)
 */

#include <stdlib.h>
#include "typesf2c.h"

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define hcurve_map_ HCURVE_MAP
#endif


typedef struct rect_struct {
	struct rect_struct *lo,*ro,*lu,*ru,*sp;
	int x,y,g,h;
} Rect_Type;

extern int xcor(Rect_Type *, int);
extern int ycor(Rect_Type *, int);


static int nx,ny,nx1,ny1;
static int *bfel;
static int *bfeld;

static Rect_Type *fel;
static Rect_Type *feld;

extern Rect_Type *fels(int);
extern Rect_Type *felds(int, int);



int xcor(Rect_Type *r, int i) {

    if (r->x == 1) {return 0;} else {
     if (r->y == 1) {return i;} else {
      if (2*r->lo->y == r->y) {
        if (i < r->lo->g - r->lo->h) {return r->lo->x - 1- xcor(r->lo,i + r->lo->h);} else {
         if (i < r->lu->g + r->lo->g - r->lo->h) {
          return r->lu->x - 1 - xcor(r->lu,r->lu->g + r->lo->g -r->lo->h -i-1);
         } else {
          if (i < r->ru->g + r->lu->g + r->lo->g -r->lo->h) {
           return r->lu->x + xcor(r->ru,i - r->lu->g- r->lo->g + r->lo->h);
          } else {
           if (i < r->g - r->lo->h) {
            return r->lo->x + xcor(r->ro,r->g -r->lo->h  -i-1);
           } else {
            return r->lo->x - 1 - xcor(r->lo,i + r->lo->h -r->g);
      } }}}} else {
       if (2*r->lo->x == r->x) { return ycor(r->sp,r->g-i);
       } else {
        if (i < r->lo->h) {return xcor(r->lo,i);} else {
         if (i < r->lu->g + r->lo->h) {
          return r->lu->x - 1 - xcor(r->lu,r->lu->g + r->lo->h -i-1);
         } else {
          if (i < r->lu->g + r->ru->g + r->lo->h) {
           return r->lu->x + xcor(r->ru,i- r->lu->g -r->lo->h);
          } else {
           if (i < r->g - r->lo->g + r->lo->h) {
            return r->lo->x + xcor(r->ro,r->g -r->lo->g + r->lo->h  -i-1);
           } else {
            return xcor(r->lo,i + r->lo->g - r->g);
    }}}}}}}}
}

int ycor(Rect_Type *r, int i) {
    if (r->x == 1) {return i;} else {
     if (r->y == 1) {return 0;} else {
      if (2*r->lo->y == r->y) {
        if (i < r->lo->g - r->lo->h) {return r->lo->y - 1 - ycor(r->lo, i+r->lo->h);} else {
         if (i < r->lu->g + r->lo->g -r->lo->h) {
          return r->lo->y + ycor(r->lu, r->lu->g + r->lo->g -r->lo->h -i-1);
         } else {
          if (i < r->ru->g + r->lu->g + r->lo->g -r->lo->h) {
           return r->ro->y + ycor(r->ru, i - r->lu->g - r->lo->g + r->lo->h);
          } else {
           if (i < r->g - r->lo->h) {
            return r->ro->y - 1 - ycor(r->ro, r->g -r->lo->h  -i-1);
           } else {
            return r->lo->y - 1 - ycor(r->lo, i + r->lo->h - r->g);
      } }}}} else {
       if (2*r->lo->x == r->x) { return xcor(r->sp, r->g - i);
       } else {
        if (i < r->lo->h) {return ycor(r->lo, i);} else {
         if (i < r->lu->g + r->lo->h) {
          return r->lo->y + ycor(r->lu, r->lu->g + r->lo->h -i-1);
         } else {
          if (i < r->lu->g + r->ru->g + r->lo->h) {
           return r->ro->y + ycor(r->ru, i - r->lu->g - r->lo->h);
          } else {
           if (i < r->g - r->lo->g + r->lo->h) {
            return r->ro->y - 1 - ycor(r->ro, r->g -r->lo->g + r->lo->h  -i-1);
           } else {
            return ycor(r->lo, i + r->lo->g - r->g);
    }}}}}}}}


}



Rect_Type *fels(int xm) {
	Rect_Type *r;

    if (bfel[xm]) {return &fel[xm];} else {
      r = &fel[xm]; bfel[xm]=1;
      r->g=xm; r->x=xm; r->y=1; r->h=1;
      r->lo = felds(0,0); r->ro = felds(0,0);
      r->lu = felds(0,0); r->ru = felds(0,0);
      return r;
    }
}

Rect_Type *felds(int xm, int ym){

	Rect_Type *r;
	int xn, yn;

	if (bfeld[xm+nx1*ym]) {return &feld[xm+nx1*ym];} else {
		r = &feld[xm+nx1*ym]; bfeld[xm+nx1*ym]=1;
        r->x=xm; r->y=ym; xn=xm/2; yn=ym/2; r->g=xm*ym;
        if (xn==(xn/2)*2) {xn= xm -xn;};
        if (yn==(yn/2)*2) {yn= ym -yn;};
        if ((xm-xn)*(ym-yn)==0) {r->h=xm*ym-1;
        if (xm*ym==0) {r->h+=1;};
        r->lo = felds(0,0); r->ro = felds(0,0);
        r->lu = felds(0,0); r->ru = felds(0,0);
        } else {
         r->lo = felds(xn,yn);    r->ro = felds(xm-xn,yn);
         r->lu = felds(xn,ym-yn); r->ru = felds(xm-xn,ym-yn);
         if (xm == 3) {
          if (2*yn != ym) {
           r->lo = fels(3);       r->ro = felds(0,1);
           r->lu = felds(0,ym-1); r->ru = felds(3,ym-1);
          } else {
           if (2*(yn/2) != yn) {
            r->lo = felds(3,yn); r->ro = felds(0,yn);
            r->lu = felds(0,yn); r->ru = felds(3,yn);
           } else {
            r->lo = felds(3,yn); r->ro = felds(0,yn);
            r->lu = felds(3,yn); r->ru = felds(0,yn);
         }}} else {
          if ((xm == 2) && (2*yn == ym)) {
            r->lo = felds(2,yn);    r->ro = felds(0,yn);
            r->lu = felds(0,ym-yn); r->ru = felds(2,ym-yn);
          } else {
           if ((ym == 3) && (2*xn != xm)) {
            r->lo = fels(xm);      r->ro = felds(0,1);
            r->lu = felds(0,ym-1); r->ru = felds(xm,ym-1);
         }}}
         if ((2*r->lo->y) == ym) {r->h = (r->lo->g)-(r->lo->h);} else {r->h = r->lo->h;};
         r->h=r->h+r->lu->g+r->ru->h;
         if ((xm == 3) && ((2*(yn/2)) == yn) && ((2*yn) == ym)) {
          r->h = r->h - r->lu->h +1;
        }}
        if ((2*r->lo->x == xm) && (2*r->lo->y != ym)) {
         r->sp = felds(ym,xm);
         r->h = r->sp->g - r->sp->h;
        }
        return r;
     }
}



void FATR hcurve_map_(sizex_ptr,sizey_ptr,map)
Integer     *sizex_ptr,*sizey_ptr;
Integer     map[];
{

	int i,ii,jj;
	Rect_Type *r;

	nx = *sizex_ptr;
	ny = *sizey_ptr;
	nx1 = nx+1;
	ny1 = ny+1;

	bfel  = (int *) malloc(2*nx1*sizeof(int));
	bfeld = (int *) malloc(2*nx1*ny1*sizeof(int));
	for (i=0; i<2*nx1; ++i)       bfel[i]  = 0;
        for (i=0; i<(2*nx1*ny1); ++i) bfeld[i] = 0;

	fel  = (Rect_Type *) malloc(2*nx1*sizeof(Rect_Type));
	feld = (Rect_Type *) malloc(2*nx1*ny1*sizeof(Rect_Type));

	r = felds(nx,ny);

	for (i=0; i< (nx*ny); ++i) {
		ii = xcor(r,i);
		jj = ycor(r,i);
		map[ii+nx*jj] = i;
	}
	free(feld);
	free(fel);
	free(bfeld);
	free(bfel);
}

/* $Id$ */
