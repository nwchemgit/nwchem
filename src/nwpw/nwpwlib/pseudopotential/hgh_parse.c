/*
 $Id$
*/

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
#define hgh_parse_ HGH_PARSE
#endif

void FATR hgh_parse_
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
 Integer *n2)
{
    char *sdir_name    = _fcdtocp(fcd_sdir_name);
    char *dir_name     = _fcdtocp(fcd_dir_name);
    char *in_filename  = _fcdtocp(fcd_in_filename);
    char *out_filename = _fcdtocp(fcd_out_filename);

#else
(debug_ptr,lmax_ptr,locp_ptr,rlocal_ptr,
 sdir_name,n9,dir_name,n0,in_filename,n1,out_filename,n2)
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
{

#endif

    int      debug,done;
    int      lmax_out,locp_out;
    double   rlocal_out;

    int      lmax;
    int      nmax_l[4];

    char     atom[10];

    int      Zion;                  /* local psp parameters          */
    double   rloc,C1,C2,C3,C4;

    double   r[4];                  /* projector radii               */

    double   H[3][4];               /* diagonal overlap coefficients */
    double   K[3][4];               /* diagonal overlap coefficients */

    int    i,p,p1;
    char   *w,*tc;
    FILE   *fp;

    char     comment[255];
    char     line[255];
    int      argc,value;
    char     words[20][80];



    int		m9 = ((int) (*n9));
    int		m0 = ((int) (*n0));
    int		m1 = ((int) (*n1));
    int		m2 = ((int) (*n2));
    char *infile  = (char *) malloc(m9+m1+5);
    char *outfile = (char *) malloc(m0+m2+5);

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
    strcat(outfile,"/");
    outfile[m0+1] = '\0';
    strncat(outfile,out_filename,m2);
    outfile[m0+m2+1] = '\0';

    Zion = 0;
    rloc = 0.0;
    C1   = 0.0;
    C2   = 0.0;
    C3   = 0.0;
    C4   = 0.0;


    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    r[3] = 0.0;

    H[0][0] = 0.0;  /*s*/
    H[1][0] = 0.0;
    H[2][0] = 0.0;

    H[0][1] = 0.0;  /*p*/
    H[1][1] = 0.0;
    H[2][1] = 0.0;

    H[0][2] = 0.0;  /*d*/
    H[1][2] = 0.0;
    H[2][2] = 0.0;

    H[0][3] = 0.0;  /*f*/
    H[1][3] = 0.0;
    H[2][3] = 0.0;

    K[0][0] = 0.0;  /*s*/
    K[1][0] = 0.0;
    K[2][0] = 0.0;

    K[0][1] = 0.0;  /*p*/
    K[1][1] = 0.0;
    K[2][1] = 0.0;

    K[0][2] = 0.0; /*d*/
    K[1][2] = 0.0;
    K[2][2] = 0.0;

    K[0][3] = 0.0;  /*f*/
    K[1][3] = 0.0;
    K[2][3] = 0.0;

    /* find the comment */
    strcpy(comment,"Hartwigsen, Goedecker and Hutter pseudopotential");
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

    /* Read HGH psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<HGH>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <HGH> section not found\n");
        fclose(fp);
        exit(99);
    }

    argc = to_eoln(fp);

    if (!get_string(fp,atom)) printf("NO ATOM NAME\n");
    if (!get_int(fp,&Zion)) printf("NO ZION\n");

    if (!get_float(fp,&rloc)) printf("NO rlocN\n");
    if (!get_float(fp,&C1)) C1 = 0.0;
    if (!get_float(fp,&C2)) C2 = 0.0;
    if (!get_float(fp,&C3)) C3 = 0.0;
    if (!get_float(fp,&C4)) C4 = 0.0;

    argc = to_eoln(fp);
    done = get_end(fp);
    i = 0;
    if (!done)
    {
        if (!get_float(fp,&(r[i])))    r[i]    = 0.0;
        if (!get_float(fp,&(H[0][i]))) H[0][i] = 0.0;
        if (!get_float(fp,&(H[1][i]))) H[1][i] = 0.0;
        if (!get_float(fp,&(H[2][i]))) H[2][i] = 0.0;
        argc = to_eoln(fp);
        done = get_end(fp);
        ++i;
        while (!done)
        {
            if (!get_float(fp,&(r[i])))    r[i]    = 0.0;
            if (!get_float(fp,&(H[0][i]))) H[0][i] = 0.0;
            if (!get_float(fp,&(H[1][i]))) H[1][i] = 0.0;
            if (!get_float(fp,&(H[2][i]))) H[2][i] = 0.0;
            argc = to_eoln(fp);
            done = get_end(fp);
            if (!done)
            {
                if (!get_float(fp,&(K[0][i]))) K[0][i] = 0.0;
                if (!get_float(fp,&(K[1][i]))) K[1][i] = 0.0;
                if (!get_float(fp,&(K[2][i]))) K[2][i] = 0.0;
                argc = to_eoln(fp);
                done = get_end(fp);
            }
            ++i;
        }
    }
    lmax = i-1;
    fclose(fp);


    /* write outfile */
    fp = fopen(outfile,"w+");
    fprintf(fp,"%d\n",1);  /* set psp_type */
    fprintf(fp,"%s\n",atom);
    fprintf(fp,"%d\n",Zion);
    if (lmax<0) fprintf(fp,"%d",0);
    else fprintf(fp,"%d\n",lmax);
    fprintf(fp,"%lf  %lf %lf %lf %lf\n",rloc,C1,C2,C3,C4);

    if (lmax>=0)
    {
        fprintf(fp,"%lf  %lf %lf %lf\n",r[0],H[0][0],H[1][0],H[2][0]);
        for (i=1; i<=lmax; ++i)
        {
            fprintf(fp,"%lf %lf %lf %lf\n",r[i],H[0][i],H[1][i],H[2][i]);
            fprintf(fp,"%lf %lf %lf\n",         K[0][i],K[1][i],K[2][i]);
        }
    }
    fprintf(fp,"%s\n",comment);
    fclose(fp);


    if (debug)
    {
        printf("HGH pseudopotential Parameters\n\n");
        printf("atom : %s\n",atom);
        printf("Zion : %d\n",Zion);
        printf(" lmax: %d\n\n",lmax);
        printf(" vloc: %lf    %lf %lf %lf %lf\n\n",rloc,C1,C2,C3,C4);

        if (lmax>=0)
        {
            printf("l=%d   r=%lf \t H= %lf %lf %lf\n\n",0,r[0],H[0][0],H[1][0],H[2][0]);
            for (i=1; i<=lmax; ++i)
            {
                printf("l=%d   r=%lf \t H= %lf %lf %lf\n",i,r[i],H[0][i],H[1][i],H[2][i]);
                printf("           \t\t K= %lf %lf %lf\n\n",     K[0][i],K[1][i],K[2][i]);
            }
        }
    }

    /* free malloc memory */
    free(infile);
    free(outfile);
    free(full_filename);

    fflush(stdout);
    return;

} /* main */
