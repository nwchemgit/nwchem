/*
 $Id$
   spline.c -
    Based on algorithms in hamman's code.
*/

#include 	<stdlib.h>
#include	<stdio.h>
#include	<string.h>

#include	"get_word.h"
#include	"grids.h"
#include	"loggrid.h"

static	int	nrl=1501;
static	double	drl=0.02;
static  int     zeroflag=0;

static  int     *nl;

void	init_Linear(char *filename)
{
    FILE	*fp;
    char *w;

    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w != ((char *) EOF)) && (strcmp("<linear>",w) != 0))
        w = get_word(fp);
    if (w!=((char *) EOF))
    {
        fscanf(fp,"%d %lf",&nrl,&drl);
    }
    fclose(fp);
    nl   = (int *) malloc(nrl*sizeof(int));


    fp = fopen(filename,"r+");
    w = get_word(fp);
    while ((w != ((char *) EOF)) && (strcmp("<fixzero>",w) != 0))
        w = get_word(fp);
    if (w!=((char *) EOF))
    {
        fscanf(fp,"%d",&zeroflag);
    }
    fclose(fp);


}

void	end_Linear()
{
    free(nl);
}

int	nrl_Linear()
{
    return nrl;
}

double	drl_Linear()
{
    return drl;
}


/********************************
 *				*
 *	     Spline		*
 *				*
 ********************************/

void Spline(x,y,n,yp1,ypn,y2)
double 	x[],
y[];
int	n;
double	yp1;
double	ypn;
double	y2[];
{
    int	i,k;
    double sig,qn,un,p;
    double *u;

    u = alloc_Grid();
    if (yp1 > 0.99e30)
    {
        y2[0] = 0.0;
        u[0]  = 0.0;
    }
    else
    {
        y2[0] = -0.5;
        u[0] = 3.0/(x[1]-x[0]) * ((y[1]-y[0])/(x[1]-x[0]) - yp1);
    }

    for (i=1; i<(n-1); ++i)
    {
        sig = (x[i]-x[i-1])/(x[i+1] - x[i-1]);
        p   = sig*y2[i-1] + 2.0;
        y2[i] = (sig-1.0)/p;
        u[i] = ( 6.0 *
                 ((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]))
                 /(x[i+1]-x[i-1])
                 - sig*u[i-1]
               ) / p;
    }

    if (ypn > 0.99e30)
    {
        qn = 0.0;
        un = 0.0;
    }
    else
    {
        qn = 0.5;
        un = 3.0/(x[n-1]-x[n-2]) * (ypn - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }

    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2] + 1.0);
    for (k=n-2; k>=0; --k)
        y2[k] = y2[k]*y2[k+1] + u[k];

    dealloc_Grid(u);

} /* Spline */

/********************************
 *				*
 *	     Splint		*
 *				*
 ********************************/


double	Splint(xa,ya,y2a,n,nx,x)
double	xa[];
double	ya[];
double	y2a[];
int	n;
int	nx;
double	x;
{
    int khi,klo;
    double h,a,b;
    double y;

    khi = nx+1;
    klo = nx;

    while ( (xa[klo] > x) || ( xa[khi] < x))
    {
        /*
              printf("Error in Splint ");
              printf("%d ->  %le %le %le",klo,x,xa[klo],xa[khi]);
        */
        if (xa[klo] > x)
        {
            --klo;
            --khi;
            /*
                     printf("   <\n");
            */
        }
        if (xa[khi] < x)
        {
            ++klo;
            ++khi;
            /*
                     printf("   >\n");
            */
        }
    }
    h = xa[khi] - xa[klo];
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    y = a*ya[klo] + b*ya[khi]
        + ( (a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi] ) * (h*h)/6.0;

    return y;

} /* Splint */

/********************************
 *				*
 *	Log_to_Linear		*
 *				*
 ********************************/

void	Log_to_Linear(ulog,rl,ulin)
double	ulog[];
double	rl[];
double	ulin[];
{
    int i,Ngrid;
    /*
      int	nl[5000];
    */
    double r0,al;
    double *r;
    double *tmp;

    r = r_LogGrid();
    r0 = r[0];
    al = log_amesh_LogGrid();

    rl[0] = r[0];
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((double) i);
        nl[i] = rint(log(rl[i]/r0)/al -0.5);
    }

    Ngrid = N_LogGrid();
    tmp = alloc_Grid();

    Spline(r,ulog,Ngrid-4,0.0,0.0,tmp);
    ulin[0] = ulog[0];
    for (i=1; i<nrl; ++i)
        ulin[i] = Splint(r,ulog,tmp,Ngrid-4,nl[i],rl[i]);

    dealloc_Grid(tmp);

} /* Log_to_Linear */


/********************************
 *                              *
 *      Log_to_Linear_zero      *
 *                              *
 ********************************/

void    Log_to_Linear_zero(ulog,rl,ulin)
double  ulog[];
double  rl[];
double  ulin[];
{
    int i,Ngrid;
    /*
      int        nl[5000];
    */
    double r0,al;
    double *r;
    double *tmp;

    r = r_LogGrid();
    r0 = r[0];
    al = log_amesh_LogGrid();

    rl[0] = r[0];
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((double) i);
        nl[i] = rint(log(rl[i]/r0)/al -0.5);
    }

    Ngrid = N_LogGrid();
    tmp = alloc_Grid();

    Spline(r,ulog,Ngrid-4,0.0,0.0,tmp);
    ulin[0] = ulog[0];
    for (i=1; i<nrl; ++i)
        ulin[i] = Splint(r,ulog,tmp,Ngrid-4,nl[i],rl[i]);

    dealloc_Grid(tmp);

    if (zeroflag)
    {
        ulin[0] = ulin[1]
                  + (r0-rl[1])*(ulin[2]-ulin[1])/(rl[2]-rl[1]);
    }


} /* Log_to_Linear_zero */





double	nm2(int n, double *y, double h)
{
    int k;
    double sum,sum1,sum2;

    sum1 = 0.0;
    sum2 = 0.0;
    for (k=0; k<n; k=k+2)
        sum1 += y[k]*y[k];
    for (k=1; k<n; k=k+2)
        sum2 += y[k]*y[k];

    sum = 2.0*sum1 + 4.0*sum2 - y[0]*y[0] - y[n-1]*y[n-1];
    sum = sum*h/3.0;

    return sum;
} /* norm2 */

void	normalize_Linear(double	*wl)
{
    int	  k;
    double norm;

    norm = nm2(nrl,wl,drl);
    norm = 1.0/sqrt(norm);
    for (k=0; k<nrl; ++k)
        wl[k] = norm*wl[k];
}


void	normalize_Linear2(double *wl, double *ul)
{
    int	  k;
    double norm;

    norm = nm2(nrl,wl,drl);
    norm = 1.0/sqrt(norm);
    for (k=0; k<nrl; ++k)
    {
        wl[k] = norm*wl[k];
        ul[k] = norm*ul[k];
    }
}
