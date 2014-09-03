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
#define qmmm_parse_ QMMM_PARSE
#endif

void FATR qmmm_parse_
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

    int      i,k,p,p1;
    int      nrl;
    double       *rl,
    **psil,
    **pspl;
    double   drl,rmax;

    int      lmaxp,n_sigma;
    double   rc,rc1,rc2,rr1,rr2,ttt,sss,s_sigma;


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



    /* find the comment */
    strcpy(comment,"QMMM formatted  pseudopotential");
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
    drl  = rmax/((double)(nrl-1));

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




    /* Read QMMM psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<QMMM>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <QMMM> section not found\n");
        fclose(fp);
        exit(99);
    }

    argc = to_eoln(fp);

    if (!get_string(fp,atom))  printf("NO Atom name\n");
    if (!get_float(fp,&Zion))  printf("NO Zion\n");
    if (!get_int(fp,&n_sigma)) printf("NO n_sigma\n");
    if (!get_float(fp,&rc))    printf("NO rc\n");
    fclose(fp);

    lmax  = 0;
    lmaxp = lmax+1;


    /* generate linear meshes */
    rl       = (double *) malloc(nrl*sizeof(double));
    psil     = (double **) malloc(lmaxp*sizeof(double*));
    pspl     = (double **) malloc(lmaxp*sizeof(double*));

    rl[0] = 0.00004167;
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((double) i);
    }

    /* generate potential */
    rc1 = 1.0;
    for (i=0; i<n_sigma; ++i) rc1 *= rc;
    rc2 = rc1*rc;

    pspl[0] = (double *) malloc(nrl*sizeof(double));
    psil[0] = (double *) malloc(nrl*sizeof(double));
    if (Zion>0.0)
    {
        for (i=0; i<nrl; ++i)
        {
            rr1 = 1.0; for (p=0; p<n_sigma; ++p) rr1 *= rl[i];
            rr2 = rr1*rl[i];
            ttt = (rc1 - rr1);
            sss = (-rc2 - rr2);
            pspl[0][i] = -Zion*(ttt/sss);
            psil[0][i] = 0.0;
        }
    }
    else
    {
        for (i=0; i<nrl; ++i)
        {
            rr1 = 1.0; for (p=0; p<n_sigma; ++p) rr1 *= rl[i];
            rr2 = rr1*rl[i];
            ttt = (rc1 - rr1);
            sss = (rc2 - rr2);
            /* l'Hopital */
            if (fabs(sss)<1.0e-9)
                pspl[0][i] = (-Zion/rc) * ((double) n_sigma)/((double) (n_sigma+1));
            else
                pspl[0][i] = -Zion*(ttt/sss);
            psil[0][i] = 0.0;
        }

    }


    /* write outfile */
    fp = fopen(outfile,"w+");
    fprintf(fp,"%s\n",atom_out);
    fprintf(fp,"%lf %lf %d   %d %d %lf\n",Zion,0.0,lmax,lmax_out,locp_out,rlocal_out);
    fprintf(fp,"%lf\n", rc);
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
            fprintf(fp," %12.8lf %12.8lf", psil[p][k],psil[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(psil[p]);
    free(psil);

    fclose(fp);


    if (debug)
    {
        printf("QMMM pseudopotential Parameters\n\n");
        printf("atom : %s\n",atom);
        printf("Zion : %lf\n",Zion);
        printf(" lmax: %d\n",lmax);
        printf(" locp: %d\n",locp_out);
        printf(" rlocal: %lf\n\n",rlocal_out);

    }

    /* free malloc memory */
    free(infile);
    free(outfile);
    free(full_filename);
    free(atom_out);

    fflush(stdout);
    return;

} /* main */




