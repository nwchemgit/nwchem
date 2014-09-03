/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "typesf2c.h"
#include "get_word.h"

extern double tetercc();
extern double cpi_Splint();
extern void   cpi_Spline();


#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define teter_parse_ TETER_PARSE
#endif

void FATR teter_parse_
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


    double   zatom,zion;      /* local psp parameters          */
    double over_fourpi;

    int      *nl;
    int      i,k,l,p,p1;
    int      Ngrid,nrl;
    double   *rgrid,*psi,*psp;
    double       *rl, *tmp, *tmp2, *sc_rho, *sc_rhol, *sc_drho, *sc_drhol,
    **psil,
    **pspl;
    double   drl,rmax;

    int      lmax,locp,lmaxp;
    double   r0,xx;
    int      n[10];
    int      pspdat,pspcode,pspxc;
    double   r2well,rcore[10],e99,e999;
    double   rchrg,fchrg,qchrg,pi;


    char   *w,*tc;
    FILE   *fp;

    char     comment[255];
    int      argc;



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

    pi          = 4.0*atan(1.0);
    over_fourpi = 1.0/(4.0*pi);



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







    /* Read TETER psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<TETER>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <TETER> section not found\n");
        fclose(fp);
        exit(99);
    }

    argc = to_eoln(fp);
    argc= get_line(fp,comment,255);

    fscanf(fp,"%lf %lf %d",&zatom,&zion,&pspdat);
    argc=to_eoln(fp);
    fscanf(fp,"%d %d %d %d %d %lf",&pspcode,&pspxc,&lmax,&locp,&Ngrid,&r2well);
    lmaxp = lmax+1;
    argc=to_eoln(fp);


    for (p=0; p<=lmax; ++p)
    {
        fscanf(fp,"%d %lf %lf %d %lf",&l,&e99,&e999,&(n[p]),&(rcore[p]));
        to_eoln(fp);
        to_eoln(fp);
    }
    fscanf(fp,"%lf %lf %lf",&rchrg,&fchrg,&qchrg);




    psi     = (double *) malloc(Ngrid*sizeof(double));
    psp     = (double *) malloc(Ngrid*sizeof(double));
    rgrid   = (double *) malloc(Ngrid*sizeof(double));
    tmp     = (double *) malloc(Ngrid*sizeof(double));
    tmp2    = (double *) malloc(Ngrid*sizeof(double));
    sc_rho  = (double *) malloc(Ngrid*sizeof(double));
    sc_drho = (double *) malloc(Ngrid*sizeof(double));


    /* define Teter grid */
    for (i=0; i<Ngrid; ++i)
    {
        xx = ((double) i);
        xx=xx/((double) (Ngrid-1));
        xx = (xx+0.01);
        xx = xx*xx*xx*xx*xx;
        rgrid[i]=100.0*xx-1.0e-8;
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

    r0    = rgrid[280];
    rl[0] = rgrid[280];
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((double) i);
        xx = (rl[i] + 1.0e-8)/100.0;
        xx = pow(xx,0.2);
        xx = xx-0.01;
        xx = (Ngrid-1)*xx;
        nl[i] = rint(xx-0.5);
    }


    /* read in pseudopotentials */
    for (p=0; p<=lmax; ++p)
    {
        pspl[p] = (double *) malloc(nrl*sizeof(double));

        to_eoln(fp);
        to_eoln(fp);
        for (i=0; i<Ngrid; ++i)  fscanf(fp,"%lf",&(psp[i]));


        cpi_Spline(rgrid,psp,Ngrid-4,0.0,0.0,tmp,tmp2);
        pspl[p][0] = psp[280];
        for (i=1; i<nrl; ++i)
        {
            pspl[p][i] = cpi_Splint(rgrid,psp,tmp,Ngrid-4,nl[i],rl[i]);
        }
    }

    /* read in wavefunctions */
    for (p=0; p<=lmax; ++p)
    {
        psil[p] = (double *) malloc(nrl*sizeof(double));

        to_eoln(fp);
        to_eoln(fp);
        for (i=0; i<Ngrid; ++i)  fscanf(fp,"%lf",&(psi[i]));


        cpi_Spline(rgrid,psi,Ngrid-4,0.0,0.0,tmp,tmp2);
        psil[p][0] = psi[280];
        for (i=1; i<nrl; ++i)
        {
            psil[p][i] = cpi_Splint(rgrid,psi,tmp,Ngrid-4,nl[i],rl[i]);
        }
    }

    /* define semicore */
    if (rchrg>0.0)
    {
        for (i=0; i<Ngrid; ++i)
        {
            /*
              xx = rgrid[i]/(rchrg);
              gg=sin(2.0*pi*xx)/( (2.0*pi*xx)*(1.0-4.0*xx*xx)*(1.0-xx*xx) );
              gg=gg*gg;
            */
            xx = rgrid[i]/(rchrg);
            sc_rho[i] = 4*pi*fchrg*tetercc(xx);
        }

        cpi_Spline(rgrid,sc_rho,Ngrid-4,0.0,0.0,tmp,tmp2);
        sc_rhol[0] = sc_rho[280];
        for (i=1; i<nrl; ++i)
        {
            sc_rhol[i] = cpi_Splint(rgrid,sc_rho,tmp,Ngrid-4,nl[i],rl[i]);
        }

        /* define to be zero for now since it is not used */
        for (i=1; i<nrl; ++i)
        {
            sc_drhol[i] = 0.0;
        }

    }

    fclose(fp);

    free(tmp2);
    free(tmp);
    free(rgrid);
    free(psp);
    free(psi);
    free(sc_rho);
    free(sc_drho);
    free(nl);


    /* find the comment */
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



    /* write outfile */
    fp = fopen(outfile,"w+");
    fprintf(fp,"%s\n",atom_out);
    if (locp_out!=-1) locp=locp_out;
    fprintf(fp,"%lf %lf %d   %d %d %lf\n",zion,0.0,lmax,lmax_out,locp,rlocal_out);
    for (p=0; p<=lmax; ++p)
        fprintf(fp,"%lf ", rcore[p]);
    fprintf(fp,"\n");
    fprintf(fp,"%d %lf\n",nrl,drl);
    fprintf(fp,"%s",comment);


    /* appending pseudopotentials */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%14.9lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", pspl[p][k]);
        fprintf(fp,"\n");
    }

    for (p=0; p<=lmax; ++p) free(pspl[p]);
    free(pspl);

    /* appending pseudowavefunctions */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%14.9lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", psil[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(psil[p]);
    free(psil);

    /* append semicore corrections */
    if (rchrg != 0.0)
    {
        fprintf(fp,"%lf\n",rchrg);
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%14.9lf %12.8lf\n", rl[k],
                    fabs(sc_rhol[k]*over_fourpi));
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%14.9lf %12.8lf\n", rl[k],
                    (sc_drhol[k]*over_fourpi));
    }
    free(sc_rhol);
    free(sc_drhol);
    free(rl);


    fclose(fp);


    if (debug)
    {
        printf("TETER pseudopotential Parameters\n\n");
        printf("atom : %s\n",atom_out);
        printf("Zatom= %lf\n",zatom);
        printf("Zion = %lf\n",zion);
        printf(" lmax= %d\n",lmax);
        printf(" locp= %d\n",locp);
        printf(" rlocal= %lf\n\n",rlocal_out);
        printf(" rcrhg=%lf  fchrg=%lf  qchrg=%lf\n",rchrg,fchrg,qchrg);
        printf("rcore: ");
        for (p=0; p<=lmax; ++p)
            printf("%lf ", rcore[p]);
        printf("\n");
        printf(" nrl=%d drl=%lf\n",nrl,drl);
        printf("comment:%s\n",comment);

        fflush(stdout);
    }

    /* free malloc memory */
    free(infile);
    free(outfile);
    free(full_filename);
    free(atom_out);

    return;

} /* main */


double tetercc(double xx)
{

    /*The c s are coefficients for Taylor expansion of the analytic form near xx=0, 1/2, and 1. */
    double   c21,c22,c23,c24;
    double   c31,c32,c33,c34;

    /*local variables */
    double pi,gg1cc,yy;

    pi = 4.0*atan(1.0);
    c21= 4.00/9.00;
    c22= -40.00/27.00;
    c23= 20.00/3.00-16.00*pi*pi/27.00;
    c24= -4160.00/243.00+160.00*pi*pi/81.00;
    c31= 1.00/36.00;
    c32= -25.00/108.00;
    c33= 485.00/432.00-pi*pi/27.00;
    c34=-4055.00/972.00+25.00*pi*pi/81.00;


    /* Cut off beyond 3/gcut=xcccrc */
    if (xx>3.000)
        gg1cc=0.000;

    /* Take care of difficult limits near x=0, 1/2, and 1 */
    else if (fabs(xx)<=1.e-9)
        gg1cc=1.00;

    else if (fabs(xx-0.500)<=1.0e-4)
        gg1cc=c21+(xx-0.500)*(c22+(xx-0.500)*(c23+(xx-0.500)*c24));

    else if (fabs(xx-1.00)<=1.e-04)
        gg1cc=c31+(xx-1.000)*(c32+(xx-1.000)*(c33+(xx-1.000)*c34));
    else
    {
        yy=sin(2.0*pi*xx)/( (2.0*pi*xx)*(1.0-4.0*xx*xx)*(1.0-xx*xx) );
        gg1cc=yy*yy;
    }

    return gg1cc;
}

