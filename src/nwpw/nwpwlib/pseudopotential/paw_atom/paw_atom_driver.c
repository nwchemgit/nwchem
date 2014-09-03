/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "typesf2c.h"


#include "paw_atom.h"
#include "paw_sdir.h"
#include "paw_output.h"
#include "paw_basis.h"
#include "paw_scattering.h"
#include "paw_loggrid.h"

#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define teter_parse_ PAW_ATOM_DRIVER
#endif

void FATR paw_atom_driver_
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
Integer *debug_ptr;
Integer *lmax_ptr;
Integer *locp_ptr;
double  *rlocal_ptr;
char    sdir_name[];
Integer *n9;
char    dir_name[];
Integer *n0;
char    in_filename[];
Integer *n1;
char    out_filename[];
Integer *n2;
char    atom[];
Integer *n3;
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



    int          m9 = ((int) (*n9));
    int          m0 = ((int) (*n0));
    int          m1 = ((int) (*n1));
    int          m2 = ((int) (*n2));
    int          m3 = ((int) (*n3));
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


    paw_set_sdir(sdir_name,m9);
    paw_set_debug(debug);
    paw_init_paw_scattering_set();

    if (debug) printf("\ninitializing atom parameters\n");
    paw_init_atom(atom_out,infile);


    if (debug) printf("\nentering the selfconsistent loop\n");
    paw_solve_atom();
    paw_print_atom();

    paw_init_paw_atom(infile);
    paw_solve_paw_atom(infile);
    paw_generate_basis_file(outfile);

    paw_end_paw_scattering();
    paw_end_paw_basis();
    paw_end_LogGrid();

}

