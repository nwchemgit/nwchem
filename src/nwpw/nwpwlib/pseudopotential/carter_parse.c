/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "typesf2c.h"
#include "get_word.h"


#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define carter_parse_ CARTER_PARSE
#endif

void FATR carter_parse_
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
    int      Ngrid,nrl=0;
    double   *rgrid,*psi,*psp;
    double       *rl, *tmp, *tmp2, *sc_rho, *sc_rhol, *sc_drho, *sc_drhol,
    **psil,
    **pspl;
    double   drl=0.0,rmax;

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





    /* Read CARTER psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<CARTER>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <CARTER> section not found\n");
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


    fscanf(fp,"%lf %lf %lf",&rchrg,&fchrg,&qchrg);
    argc=to_eoln(fp);
    argc=to_eoln(fp);
    argc=to_eoln(fp);
    argc=to_eoln(fp);
    rcore[0] = 0.0;

    psi     = (double *) malloc(Ngrid*sizeof(double));
    psp     = (double *) malloc(Ngrid*sizeof(double));
    rgrid   = (double *) malloc(Ngrid*sizeof(double));

    /* read in rgrid and pseudopotentials */
    for (i=0; i<Ngrid; ++i)  fscanf(fp,"%d %lf %lf",&k,&(rgrid[i]),&(psp[i]));

    /* define psi=0 */
    for (i=0; i<Ngrid; ++i)  psi[i] = 0.0;


    /* define semicore */

    fclose(fp);


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
    fprintf(fp,"%d %lf\n",Ngrid,rgrid[1]);
    fprintf(fp,"%s",comment);

    /* appending pseudopotentials */
    for (k=0; k<Ngrid; ++k)
       fprintf(fp,"%14.9lf %12.8lf\n", rgrid[k],psp[k]);

    /* appending pseudowavefunctions */
    for (k=0; k<Ngrid; ++k)
       fprintf(fp,"%14.9lf %12.8lf\n", rgrid[k],psi[k]);


    /* append semicore corrections */
    free(rgrid);
    free(psp);
    free(psi);


    fclose(fp);


    if (debug)
    {
        printf("CARTER pseudopotential Parameters\n\n");
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


