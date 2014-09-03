/*
 $Id$
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "typesf2c.h"
#include "name.h"
#include "loggrid.h"
#include "spline.h"
#include "atom.h"
#include "psp.h"
#include "rpsp.h"
#include "debug.h"

/** just include the source file here,it makes no sense to have a 
 **  separate object file. 
 **/
#include "rpspsolve.c"

#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define pspsolve_ PSPSOLVE
#endif

void FATR pspsolve_
#if defined(USE_FCD)
 
  (Integer * print_ptr,
   Integer * debug_ptr,
   Integer * lmax_ptr,
   Integer * locp_ptr,
   double *rlocal_ptr,
   Integer * efg_ptr,
   const _fcd fcd_sdir_name,
   Integer * n9,
   const _fcd fcd_dir_name,
   Integer * n0,
   const _fcd fcd_in_filename,
   Integer * n1, const _fcd fcd_out_filename, Integer * n2)
{
  char *sdir_name = _fcdtocp (fcd_sdir_name);
  char *dir_name = _fcdtocp (fcd_dir_name);
  char *in_filename = _fcdtocp (fcd_in_filename);
  char *out_filename = _fcdtocp (fcd_out_filename);

#else
 
  (print_ptr, debug_ptr, lmax_ptr, locp_ptr, rlocal_ptr,efg_ptr,
   sdir_name, n9, dir_name, n0, in_filename, n1, out_filename, n2)
     Integer *print_ptr;
     Integer *debug_ptr;
     Integer *lmax_ptr;
     Integer *locp_ptr;
     double *rlocal_ptr;
     Integer *efg_ptr;
     char sdir_name[];
     Integer *n9;
     char dir_name[];
     Integer *n0;
     char in_filename[];
     Integer *n1;
     char out_filename[];
     Integer *n2;
{

#endif

  int i, j, k, l, p, Nlinear, Nvalence,Ncore,istate,mch,kb_extra;
  int debug, print;
  double *rl, *rhol, **psil_ae,**psil, **psil_extra, **pspl;
  double over_fourpi, c, x, y,nu0;
  int Ngrid;
  double *vall, *rgrid;
  char name[255];
  int lmax_out, locp_out,efg_type;
  double rlocal_out,rmax;

  FILE *fp;

  int m9 = ((int) (*n9));
  int m0 = ((int) (*n0));
  int m1 = ((int) (*n1));
  int m2 = ((int) (*n2));
  char *infile = (char *) malloc (m9 + m1 + 5);
  char *outfile = (char *) malloc (m0 + m2 + 5);
  char *full_filename = (char *) malloc (m9 + 25 + 5);


  print = *print_ptr;
  debug = *debug_ptr;
  lmax_out = *lmax_ptr;
  locp_out = *locp_ptr;
  rlocal_out = *rlocal_ptr;
  efg_type = *efg_ptr;

  (void) strncpy (infile, sdir_name, m9);
  infile[m9] = '\0';
  strcat (infile, "/");
  infile[m9 + 1] = '\0';
  strncat (infile, in_filename, m1);
  infile[m9 + m1 + 1] = '\0';

  (void) strncpy (outfile, dir_name, m0);
  outfile[m0] = '\0';
  strcat (outfile, "/");
  outfile[m0 + 1] = '\0';
  strncat (outfile, out_filename, m2);
  outfile[m0 + m2 + 1] = '\0';
  
  over_fourpi = 1.0 / (16.0 * atan (1.0));

  set_debug_print(debug);
  init_Atom(infile);
  solve_Atom();
  if (debug)
    print_Atom (stdout);

  if (isRelativistic_Atom())
    {
      free(infile);
      free(outfile);
      free(full_filename);
      fprintf(stderr,"Relativstic Atom!\n");
      rpspsolve_(print_ptr, debug_ptr, lmax_ptr, locp_ptr, rlocal_ptr,
		     sdir_name, n9, dir_name, n0, in_filename, n1,
		     out_filename, n2);
      return;
    }
  init_Psp(infile);
  solve_Psp();
  if (debug)
    print_Psp(stdout);
  init_Linear(infile);
  /* allocate linear meshes */
  Nvalence = Nvalence_Psp();
  Ncore    = Ncore_Atom();
  Nlinear = nrl_Linear();
  psil    = (double **) malloc(Nvalence * sizeof (double *));
  psil_ae = (double **) malloc(Nvalence * sizeof (double *));
  pspl    = (double **) malloc(Nvalence * sizeof (double *));
  for (p = 0; p < Nvalence; ++p)
    {
      psil[p]    = (double *) malloc(Nlinear * sizeof (double));
      psil_ae[p] = (double *) malloc(Nlinear * sizeof (double));
      pspl[p]    = (double *) malloc(Nlinear * sizeof (double));
    }
  rl = (double *)   malloc(Nlinear * sizeof (double));
  rhol = (double *) malloc(Nlinear * sizeof (double));

  /* Norm-conserving output */
  if (NormConserving_Psp())
    {
      for (p = 0; p < Nvalence; ++p)
	{
	  Log_to_Linear(r_psi_Psp(p),  rl, psil[p]);
	  Log_to_Linear_zero(V_Psp(p), rl, pspl[p]);

	  /* normalize scattering state */
	  if (fill_Psp(p) == 0.0)
	  {
             rgrid = r_LogGrid();
             Ngrid = N_LogGrid();
             istate = Nvalence_Atom() + Ncore_Atom();
             rmax = 20.0;
             solve_Scattering_State_Atom(n_Psp(p),l_Psp(p),eigenvalue_Psp(p),rmax);
     
             rmax = 2.5*rcut_Psp(p);
             mch = rint(log(rmax/rgrid[0])/log_amesh_LogGrid());
             nu0 = r_psi_Psp(p)[mch]/r_psi_Atom(istate)[mch];
             for (i=0; i<Ngrid; ++i)
                r_psi_Atom(istate)[i] *= nu0;
	     Log_to_Linear(r_psi_Atom(istate), rl, psil_ae[p]);
	     normalize_Linear2(psil[p],psil_ae[p]);
	  }
          else
	     Log_to_Linear(r_psi_Atom(p+Ncore), rl, psil_ae[p]);
          
	}
      Log_to_Linear(rho_Psp(), rl, rhol);

      kb_extra = kb_extra_Psp();
      if (kb_extra>0)
      {
         psil_extra = (double **) malloc(kb_extra*sizeof(double*));
         for (p=0; p<kb_extra; ++p)
         {
            psil_extra[p] = (double *) malloc(Nlinear*sizeof(double));
            Log_to_Linear(r_psi_extra_Psp(p),rl,psil_extra[p]);
         }
      }


      if (debug)
	{
	  /* output pseudowavefunctions argv[1].psw */
	  strcpy(name, name_Atom());
	  strcat(name, ".psw.plt");
	  full_filename[0] = '\0';
	  strncpy(full_filename, sdir_name, m9);
	  full_filename[m9] = '\0';
	  strcat(full_filename, "/");
	  full_filename[m9 + 1] = '\0';
	  strcat(full_filename, name);

	  printf("Outputing pseudowavefunctions: %s\n", full_filename);
	  fp = fopen(full_filename, "w+");
	  for (k=0; k<Nlinear; ++k)
	    {
	      fprintf(fp, "%12.8lf", rl[k]);
	      for (p = 0; p < Nvalence; ++p)
		fprintf(fp, " %12.8lf", psil[p][k]);
              for (p=0; p<kb_extra; ++p)
                 fprintf(fp," %12.8lf", psil_extra[p][k]);

	      fprintf(fp, "\n");
	    }
	  fclose(fp);


	  /* output pseudopotentials argv[1].psp */
	  strcpy(name, name_Atom());
	  strcat(name, ".psp.plt");
	  full_filename[0] = '\0';
	  strncpy (full_filename, sdir_name, m9);
	  full_filename[m9] = '\0';
	  strcat (full_filename, "/");
	  full_filename[m9 + 1] = '\0';
	  strcat (full_filename, name);

	  printf("Outputing pseudopotentials: %s\n", full_filename);
	  fp = fopen(full_filename, "w+");
	  for (k = 0; k < Nlinear; ++k)
	    {
	      fprintf(fp, "%12.8lf", rl[k]);
	      for (p = 0; p < Nvalence; ++p)
		fprintf(fp, " %12.8lf", pspl[p][k]);
	      fprintf(fp, "\n");
	    }
	  fclose(fp);


	  /* output pseudodensity infile.psd */
	  strcpy(name, name_Atom());
	  strcat(name, ".psd.plt");
	  full_filename[0] = '\0';
	  strncpy(full_filename, sdir_name, m9);
	  full_filename[m9] = '\0';
	  strcat(full_filename, "/");
	  full_filename[m9 + 1] = '\0';
	  strcat(full_filename, name);

	  fp = fopen(full_filename, "w+");
	  for (k = 0; k < Nlinear; ++k)
	    fprintf(fp, "%12.8lf %12.8lf\n", rl[k], rhol[k]);
	  fclose(fp);
	}



      /* output datafile to be used for Kleinman-Bylander input, argv[1].psp */
      if (print)
	{
	  printf(" Creating datafile for Kleinman-Bylander input: %s\n",
		  outfile);
	}
      fp = fopen(outfile, "w+");
      if (efg_type && kb_extra) 
         fprintf(fp, "92\n");
      else if (efg_type) fprintf(fp, "9\n");
      else if (kb_extra)  fprintf(fp, "2\n");
      fprintf(fp, "%s\n", name_Atom());

      fprintf(fp, "%lf %lf %d   %d %d %lf\n", Zion_Psp(), Amass_Atom(),
	       lmax_Psp(), lmax_out, locp_out, rlocal_out);
      for (p = 0; p <= lmax_Psp(); ++p)
	fprintf(fp, "%lf ", rcut_Psp(p));
      fprintf(fp, "\n");
      if (kb_extra)
      {
         fprintf(fp, "%d\n", kb_extra);
         for (p = 0; p <= lmax_Psp(); ++p)
            fprintf(fp, "%d ", kb_expansion_Psp(p));
         fprintf(fp, "\n");
      }
      fprintf(fp, "%d %lf\n", nrl_Linear(), drl_Linear());
      fprintf(fp, "%s\n", comment_Psp());
       

      if (print)
	{
	  printf("  + Appending pseudopotentials:    %s thru %s\n",
		  spd_Name (0), spd_Name (lmax_Psp()));
	}
      for (k = 0; k < Nlinear; ++k)
	{
	  fprintf(fp, "%12.8lf", rl[k]);
	  for (p = 0; p <= lmax_Psp(); ++p)
	    fprintf (fp, " %12.8lf", pspl[p][k]);
	  fprintf(fp, "\n");
	}
      if (print)
	{
	  printf("  + Appending pseudowavefunctions and aewavefunctions: %s thru %s\n",
		  spd_Name(0), spd_Name(lmax_Psp()));
	}
      for (k = 0; k < Nlinear; ++k)
	{
	  fprintf(fp, "%12.8lf", rl[k]);
	  for (p = 0; p <= lmax_Psp(); ++p)
	    fprintf(fp, " %12.8lf", psil[p][k]);
	  for (p = 0; p<kb_extra; ++p)
	    fprintf(fp, " %12.8lf", psil_extra[p][k]);

          if (efg_type)
          {
	     for (p = 0; p <= lmax_Psp(); ++p)
	        fprintf(fp, " %12.8lf", psil_ae[p][k]);
          }

	  fprintf(fp, "\n");
	}

      /* append semicore corrections */
      if (r_semicore_Psp() != 0.0)
	{
	  if(print)
	    {
	      printf ("  + Appending semicore density\n");
	    }
	  Log_to_Linear(rho_semicore_Psp(), rl, rhol);
	  fprintf (fp, "%lf\n", r_semicore_Psp());
	  for (k = 0; k < Nlinear; ++k)
	    fprintf (fp, "%12.8lf %12.8lf\n", rl[k],
		     fabs(rhol[k] * over_fourpi));

	  if (print)
	    {
	      printf ("  + Appending semicore density gradient\n");
	    }
	  Log_to_Linear(drho_semicore_Psp(), rl, rhol);
	  for (k = 0; k < Nlinear; ++k)
	    fprintf (fp, "%12.8lf %12.8lf\n", rl[k], (rhol[k] * over_fourpi));
	}


      fclose(fp);
    }				/* Norm-conserving PSP output */


  /* output datafile to be used for Vanderbilt input, argv[1].vbt */
  if (Vanderbilt_Psp())
    {
      printf("Creating Vanderbilt pseudopotential input: %s\n", outfile);
      fp = fopen(outfile, "w+");
      fprintf(fp, "0 %s\n", name_Atom());
      fprintf(fp, "%lf %lf %d\n", Zion_Psp(), Amass_Atom(), lmax_Psp());

      /* output grid parameters */
      fprintf (fp, "%d %20.15le \n", Nlinear, drl_Linear());

      /* output local potential */
      Log_to_Linear(Vlocal_Psp(), rl, rhol);
      for (k = 0; k < Nlinear; ++k)
	fprintf(fp, "%20.15le  %20.15le\n", rl[k], rhol[k]);

      for (l = 0; l <= lmax_Psp(); ++l)
	{
	  /* output ns */
	  fprintf(fp, "%d\n", ns_Psp(l));

	  /* output rcut */
	  for (i = 0; i < ns_Psp (l); ++i)
	    fprintf(fp, "%le ", rcut_il_Psp(i, l));
	  fprintf(fp, "\n");

	  /* output D0   */
	  for (j = 0; j < ns_Psp(l); ++j)
	    for (i = 0; i < ns_Psp(l); ++i)
	      fprintf(fp, "%20.15le ",
		       0.5 * (D0_Psp (i, j, l) + D0_Psp (j, i, l)));
	  fprintf(fp, "\n");

	  /* output q    */
	  for (j = 0; j < ns_Psp(l); ++j)
	    for (i = 0; i < ns_Psp(l); ++i)
	      fprintf(fp, "%20.15le ",
		       0.5 * (q_Psp(i, j, l) + q_Psp(j, i, l)));
	  fprintf(fp, "\n");

	  /* output Beta */
	  for (i = 0; i < ns_Psp(l); ++i)
	    Log_to_Linear(Beta_Psp(i, l), rl, pspl[i]);
	  for (k = 0; k < Nlinear; ++k)
	    {
	      fprintf(fp, "%20.15le ", rl[k]);
	      for (i = 0; i < ns_Psp (l); ++i)
		fprintf(fp, " %20.15le", pspl[i][k]);
	      fprintf(fp, "\n");
	    }

	  /* output u    */
	  for (i = 0; i < ns_Psp(l); ++i)
	    Log_to_Linear(r_hard_psi_il_Psp(i, l), rl, pspl[i]);
	  for (k = 0; k < Nlinear; ++k)
	    {
	      fprintf(fp, "%20.15le ", rl[k]);
	      for (i = 0; i < ns_Psp(l); ++i)
		fprintf(fp, " %20.15le", pspl[i][k]);
	      fprintf(fp, "\n");
	    }

	  /* output w    */
	  for (i = 0; i < ns_Psp(l); ++i)
	    Log_to_Linear(r_psi_il_Psp(i, l), rl, pspl[i]);
	  for (k = 0; k < Nlinear; ++k)
	    {
	      fprintf(fp, "%20.15le ", rl[k]);
	      for (i = 0; i < ns_Psp(l); ++i)
		fprintf(fp, " %20.15le", pspl[i][k]);
	      fprintf(fp, "\n");
	    }
	}
      fclose(fp);

    }				/* Vanderbilt PSP output */

   /******************************************************************/
   /******************* output AE information ************************/
   /******************************************************************/


  if (debug)
    {

      /* output all-electron wavefunctions */
      printf("Outputing all-electron wavefunctions:");
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      for (p = 0; p <= (Ncore_Atom() + Nvalence_Atom()); ++p)
	{
	  sprintf(name, "%s.%1d%s", name_Atom(), n_Atom(p),
		   spd_Name(l_Atom(p)));
	  full_filename[0] = '\0';
	  strncpy (full_filename, sdir_name, m9);
	  full_filename[m9] = '\0';
	  strcat (full_filename, "/");
	  full_filename[m9 + 1] = '\0';
	  strcat  (full_filename, name);

	  printf(" %s", full_filename);
	  fp = fopen(full_filename, "w+");
	  for (k = 0; k < Ngrid; ++k)
	    fprintf(fp, "%12.8lf %12.8lf\n", rgrid[k], r_psi_Atom(p)[k]);
	  fclose (fp);
	}
      printf("\n");

      /* output density argv[1].dns */
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      strcpy(name, name_Atom());
      strcat(name, ".dns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf("Outputing atom density: %s\n", full_filename);
      fp = fopen(full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf(fp, "%12.8lf %12.8lf\n", rgrid[k], rho_Atom()[k]);
      fclose(fp);

      /* output core density argv[1].cdns */
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      strcpy (name, name_Atom());
      strcat (name, ".cdns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf("Outputing core density: %s\n", full_filename);
      fp = fopen(full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf(fp, "%12.8lf %12.8lf\n", rgrid[k], rho_core_Atom()[k]);
      fclose(fp);

      /* output core density gradient argv[1].cddns */
      vall = alloc_LogGrid();
      Derivative_LogGrid(rho_core_Atom(), vall);
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      strcpy(name, name_Atom());
      strcat(name, ".cddns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat(full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat(full_filename, name);

      printf("Outputing core density gradient: %s\n", full_filename);
      fp = fopen(full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf(fp, "%12.8lf %12.8lf\n", rgrid[k], vall[k]);
      fclose(fp);
      dealloc_LogGrid(vall);

      /* output semicore density infile.sdns */
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      strcpy (name, name_Atom());
      strcat (name, ".sdns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf("Outputing semicore density: %s\n", full_filename);
      fp = fopen(full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf(fp, "%12.8lf %12.8lf\n", rgrid[k], rho_semicore_Psp()[k]);
      fclose(fp);

      /* output semicore density gradient infile.sddns */
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      strcpy(name, name_Atom());
      strcat(name, ".sddns.plt");
      full_filename[0] = '\0';
      strncpy(full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat(full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat(full_filename, name);

      printf("Outputing semicore density gradient: %s\n", full_filename);
      fp = fopen(full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf(fp, "%12.8lf %12.8lf\n", rgrid[k], (drho_semicore_Psp())[k]);
      fclose(fp);

      /* output all-electron potential infile.pot */
      vall = Vall_Atom();
      Ngrid = N_LogGrid();
      rgrid = r_LogGrid();
      strcpy (name, name_Atom());
      strcat (name, ".pot.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf("Outputing all-electron potential(non-screened): %s\n",
	      full_filename);
      fp = fopen(full_filename, "w+");

      c = 0.0;
      for (k = 0; k < Ngrid; ++k)
	{
	  x = rgrid[k] / rcut_Psp(0);
	  x = pow (x, 3.5);
	  x = exp (-x);
	  y = 1.0 - x;
	  fprintf(fp, "%12.8lf %12.8lf %12.8lf\n", rgrid[k], vall[k],
		   y * vall[k] + c * x);
	}
      fclose(fp);
    }


  /* free malloc memory */
  free(infile);
  free(outfile);
  free(full_filename);
  for (p = 0; p < Nvalence; ++p)
    {
      free(psil[p]);
      free(psil_ae[p]);
      free(pspl[p]);
    }
  free(psil);
  free(psil_ae);
  free(pspl);
  free(rl);
  free(rhol);
  end_Linear();

  end_Psp();
  end_Atom();

}				/* main */
