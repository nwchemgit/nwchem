/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "typesf2c.h"
#include "get_word.h"

extern double cpi_Splint();
extern void   cpi_Spline();


#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define cpi_parse_ CPI_PARSE
#endif

void FATR cpi_parse_
#if defined(USE_FCD)
(Integer *debug_ptr,
 Integer *lmax_ptr,
 Integer *locp_ptr,
 double  *rlocal_ptr,
 const _fcd fcd_sdir_name,
 Integer *n9,
 const _fcd fcd_dir_name,
 Integer *n0,
 const _fcd fcd_in_filename,
 Integer *n1,
 const _fcd fcd_out_filename,
 Integer *n2,
 const _fcd fcd_atom,
 Integer *n3)
{
    char *sdir_name    = _fcdtocp(fcd_sdir_name);
    char *dir_name     = _fcdtocp(fcd_dir_name);
    char *in_filename  = _fcdtocp(fcd_in_filename);
    char *out_filename = _fcdtocp(fcd_out_filename);
    char *atom         = _fcdtocp(fcd_atom);

#else
(debug_ptr,lmax_ptr,locp_ptr,rlocal_ptr,
 sdir_name,n9,dir_name,n0,in_filename,n1,out_filename,n2,atom,n3)
Integer	*debug_ptr;
Integer	*lmax_ptr;
Integer	*locp_ptr;
double 	*rlocal_ptr;
char	sdir_name[];
Integer	*n9;
char	dir_name[];
Integer	*n0;
char	in_filename[];
Integer	*n1;
char	out_filename[];
Integer	*n2;
char	atom[];
Integer	*n3;
{

#endif

    int      debug;
    int      lmax_out,locp_out;
    double   rlocal_out;

    int      lmax;


    double   Zion;      /* local psp parameters          */
    double over_fourpi;

    int      *nl;
    int      i,j,k,p,p1;
    int      Ngrid,nrl;
    double   *rgrid,*psi,*psp;
    double       *rl, *tmp, *tmp2, *sc_rho, *sc_rhol, *sc_drho, *sc_drhol,
    **psil,
    **pspl;
    double   r,ul,vl,amesh,al,drl,r_semicore,rmax;

    int      lmaxp;
    double   dum1,dum2,dum3,dum4,r0;
    int idum;


    char   *w,*tc;
    FILE   *fp;

    char     comment[255];
    int      argc,value;



    int		m9 = ((int) (*n9));
    int		m0 = ((int) (*n0));
    int		m1 = ((int) (*n1));
    int		m2 = ((int) (*n2));
    int		m3 = ((int) (*n3));
    char *infile  = (char *) malloc(m9+m1+5);
    char *outfile = (char *) malloc(m0+m2+5);
    char *atom_out = (char *) malloc(m3+5);

    char *full_filename = (char *) malloc(m9+25+5);


    debug = *debug_ptr;
    lmax_out   = *lmax_ptr;
    locp_out   = *locp_ptr;
    rlocal_out = *rlocal_ptr;

    (void) strncpy(infile, sdir_name, m9);
    infile[m9] = '\0';
    strcat(infile,"/");
    infile[m9+1] = '\0';
    strncat(infile,in_filename,m1);
    infile[m9+m1+1] = '\0';

    (void) strncpy(outfile, dir_name, m0);
    outfile[m0] = '\0';
    (void) strcat(outfile,"/");
    outfile[m0+1] = '\0';
    (void) strncat(outfile,out_filename,m2);
    outfile[m0+m2+1] = '\0';

    (void) strncpy(atom_out, atom, m3);
    atom_out[m3] = '\0';

    over_fourpi = 1.0/(16.0*atan(1.0));


    /* find the comment */
    strcpy(comment,"CPI formatted  pseudopotential");
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<comment>",w)!=0))
        w = get_word(fp);

    if (w!=NIL)
    {
        w = get_word(fp);
        p  = 0;
        tc = comment;
        while ((w!=NIL)&&(strcmp("<end>",w) != 0))
        {
            p = (strlen(w));
            strcpy(tc, w);
            for (p1=0;p1<p; ++p1) ++tc;
            strcpy(tc, " ");
            ++tc;

            w = get_word(fp);
        }
    }
    fclose(fp);


    /* define linear grid */
    nrl  = 2001;
    rmax = 40.0;
    drl  = rmax/((double) (nrl-1));

    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w != ((char *) EOF)) && (strcmp("<linear>",w) != 0))
        w = get_word(fp);
    if (w!=((char *) EOF))
    {
        fscanf(fp,"%d %lf",&nrl,&drl);
        rmax = ((double) (nrl-1))*drl;
    }
    fclose(fp);




    /* Read CPI psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<CPI>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <CPI> section not found\n");
        fclose(fp);
        exit(99);
    }

    argc = to_eoln(fp);


    fscanf(fp,"%lf %d",&Zion,&lmaxp);
    lmax = lmaxp-1;

    fscanf(fp,"%lf %lf %lf %lf",&dum1,&dum2,&dum3,&dum4);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);
    fscanf(fp,"%lf %lf %lf ",&dum1,&dum2,&dum3);

    fscanf(fp,"%d %lf",&Ngrid,&amesh);
    al = log(amesh);

    psi     = (double *) malloc(Ngrid*sizeof(double));
    psp     = (double *) malloc(Ngrid*sizeof(double));
    rgrid   = (double *) malloc(Ngrid*sizeof(double));
    tmp     = (double *) malloc(Ngrid*sizeof(double));
    tmp2    = (double *) malloc(Ngrid*sizeof(double));
    sc_rho  = (double *) malloc(Ngrid*sizeof(double));
    sc_drho = (double *) malloc(Ngrid*sizeof(double));

    for (i=0; i<Ngrid; ++i)
    {
        fscanf(fp,"%d %lf %lf %lf",&j, &r,&ul,&vl);
        rgrid[i]  = r;
        psi[i] = ul;
        psp[i] = vl;
    }


    /* check linear grid and redefine if necessary */
    if (rmax > rgrid[Ngrid-5])
    {
        rmax = rgrid[Ngrid-5];
        drl = rmax/((double) (nrl-1));
    }



    /* generate linear meshes */
    rl       = (double *) malloc(nrl*sizeof(double));
    nl       = (int *)    malloc(nrl*sizeof(int));
    psil     = (double **) malloc(lmaxp*sizeof(double*));
    pspl     = (double **) malloc(lmaxp*sizeof(double*));
    sc_rhol  = (double *) malloc(nrl*sizeof(double));
    sc_drhol = (double *) malloc(nrl*sizeof(double));

    r0    = rgrid[0];
    rl[0] = rgrid[0];
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((double) i);
        nl[i] = rint(log(rl[i]/r0)/al -0.5);
    }


    psil[0] = (double *) malloc(nrl*sizeof(double));
    pspl[0] = (double *) malloc(nrl*sizeof(double));

    cpi_Spline(rgrid,psp,Ngrid-4,0.0,0.0,tmp,tmp2);
    pspl[0][0] = psp[0];
    for (i=1; i<nrl; ++i)
    {
        pspl[0][i] = cpi_Splint(rgrid,psp,tmp,Ngrid-4,nl[i],rl[i]);
    }

    cpi_Spline(rgrid,psi,Ngrid-4,0.0,0.0,tmp,tmp2);
    psil[0][0] = psi[0];
    for (i=1; i<nrl; ++i)
    {
        psil[0][i] = cpi_Splint(rgrid,psi,tmp,Ngrid-4,nl[i],rl[i]);
    }

    for (p=1; p<lmaxp; ++p)
    {
        fscanf(fp,"%d %lf",&idum,&dum1);

        for (i=0; i<Ngrid; ++i)
        {
            fscanf(fp,"%d %lf %lf %lf",&j, &r,&ul,&vl);
            rgrid[i]  = r;
            psi[i] = ul;
            psp[i] = vl;
        }

        psil[p] = (double *) malloc(nrl*sizeof(double));
        pspl[p] = (double *) malloc(nrl*sizeof(double));

        cpi_Spline(rgrid,psp,Ngrid-4,0.0,0.0,tmp,tmp2);
        pspl[p][0] = psp[0];
        for (i=1; i<nrl; ++i)
        {
            pspl[p][i] = cpi_Splint(rgrid,psp,tmp,Ngrid-4,nl[i],rl[i]);
        }

        cpi_Spline(rgrid,psi,Ngrid-4,0.0,0.0,tmp,tmp2);
        psil[p][0] = psi[0];
        for (i=1; i<nrl; ++i)
        {
            psil[p][i] = cpi_Splint(rgrid,psi,tmp,Ngrid-4,nl[i],rl[i]);
        }


    }

    /* read semi-core */
    r_semicore = 0.0;
    value = fscanf(fp,"%lf   %lf %lf %lf", &r,&ul,&vl,&dum1);
    if (value!=EOF)
    {
        r_semicore =  99.99; /* not known?? */
        rgrid[0]  = r; sc_rho[0] = ul; sc_drho[0] = vl;
        for (i=1; i<Ngrid; ++i)
        {
            fscanf(fp,"%lf   %lf %lf %lf", &r,&ul,&vl,&dum1);
            rgrid[i]   = r;
            sc_rho[i]  = ul;
            sc_drho[i] = vl;
        }

        cpi_Spline(rgrid,sc_rho,Ngrid-4,0.0,0.0,tmp,tmp2);
        sc_rhol[0] = sc_rho[0];
        for (i=1; i<nrl; ++i)
        {
            sc_rhol[i] = cpi_Splint(rgrid,sc_rho,tmp,Ngrid-4,nl[i],rl[i]);
        }

        cpi_Spline(rgrid,sc_drho,Ngrid-4,0.0,0.0,tmp,tmp2);
        sc_drhol[0] = sc_drho[0];
        for (i=1; i<nrl; ++i)
        {
            sc_drhol[i] = cpi_Splint(rgrid,sc_drho,tmp,Ngrid-4,nl[i],rl[i]);
        }


    }
    free(nl);
    free(rgrid);
    free(psi);
    free(psp);
    free(tmp);
    free(tmp2);
    free(sc_rho);
    free(sc_drho);



    fclose(fp);


    /* write outfile */
    fp = fopen(outfile,"w+");
    fprintf(fp,"%s\n",atom_out);
    fprintf(fp,"%lf %lf %d   %d %d %lf\n",Zion,0.0,lmax,lmax_out,locp_out,rlocal_out);
    for (p=0; p<=lmax; ++p)
        fprintf(fp,"%lf ", -1.0);
    fprintf(fp,"\n");
    fprintf(fp,"%d %lf\n",nrl,drl);
    fprintf(fp,"%s\n",comment);

    /* appending pseudopotentials */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%12.8lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", pspl[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(pspl[p]);
    free(pspl);

    /* appending pseudowavefunctions */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%12.8lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", psil[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(psil[p]);
    free(psil);


    /* append semicore corrections */
    if (r_semicore != 0.0)
    {
        fprintf(fp,"%lf\n",r_semicore);
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%12.8lf %12.8lf\n", rl[k],
                    fabs(sc_rhol[k]*over_fourpi));
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%12.8lf %12.8lf\n", rl[k],
                    (sc_drhol[k]*over_fourpi));
    }
    free(sc_drhol);
    free(sc_rhol);
    free(rl);


    fclose(fp);


    if (debug)
    {
        printf("CPI pseudopotential Parameters\n\n");
        printf("atom : %s\n",atom);
        printf("Zion : %lf\n",Zion);
        printf(" lmax: %d\n",lmax);
        printf(" locp: %d\n",locp_out);
        printf(" rlocal: %lf\n\n",rlocal_out);
        printf(" r_semicore: %lf\n",r_semicore);

    }

    /* free malloc memory */
    free(infile);
    free(outfile);
    free(full_filename);
    free(atom_out);

    fflush(stdout);
    return;

} /* main */



/********************************
 *				*
 *	     cpi_Spline		*
 *				*
 ********************************/

void cpi_Spline(x,y,n,yp1,ypn,y2,u)
double 	x[],
y[];
int	n;
double	yp1;
double	ypn;
double	y2[];
double	u[];
{
    int	i,k;
    double sig,qn,un,p;

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


} /* cpi_Spline */

/********************************
 *				*
 *	     cpi_Splint		*
 *				*
 ********************************/


double	cpi_Splint(xa,ya,y2a,n,nx,x)
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

} /* cpi_Splint */


