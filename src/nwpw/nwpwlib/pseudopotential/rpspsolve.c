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

#if defined(CRAY) || defined(CRAY_T3D)
#include <fortran.h>
#if !defined(__crayx1)
#define USE_FCD
#endif
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define rpspsolve_ RPSPSOLVE
#endif

void FATR rpspsolve_
#if defined(USE_FCD)
  (Integer * print_ptr,
   Integer * debug_ptr,
   Integer * lmax_ptr,
   Integer * locp_ptr,
   double *rlocal_ptr,
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
 
  (print_ptr, debug_ptr, lmax_ptr, locp_ptr, rlocal_ptr,
   sdir_name, n9, dir_name, n0, in_filename, n1, out_filename, n2)
     Integer *print_ptr;
     Integer *debug_ptr;
     Integer *lmax_ptr;
     Integer *locp_ptr;
     double *rlocal_ptr;
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

  int i, j, k, l, p, Nlinear, Nvalence;
  int debug, print;
  double *rl, *rhol, **psil, **pspl;
  double over_fourpi, c, x, y;
  int Ngrid;
  double *vall, *rgrid;
  char name[255];
  int lmax_out, locp_out, nvh;
  double rlocal_out, vx;

  FILE *fp;

  int m9 = ((int) (*n9));
  int m0 = ((int) (*n0));
  int m1 = ((int) (*n1));
  int m2 = ((int) (*n2));
  char *infile = (char *) malloc (m9 + m1 + 5);
  char *outfile = (char *) malloc (m0 + m2 + 5);
  char *soutfile = (char *) malloc (m0 + m2 + 9);
  char *full_filename = (char *) malloc (m9 + 25 + 5);

  print = *print_ptr;
  debug = *debug_ptr;
  lmax_out = *lmax_ptr;
  locp_out = *locp_ptr;
  rlocal_out = *rlocal_ptr;

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

  over_fourpi = 0.25/M_PI;

/********************
 *   we have already solved for the Relativsitic AE wavefunctions using atom 
 *  in the call to pspsolve
 *******************/
  set_debug_print (debug);


  init_RelPsp (infile);


  solve_RelPsp ();


  if (debug)
    print_RelPsp (stdout);


  init_Linear (infile);

  /* allocate linear meshes */
  Nvalence = Nvalence_RelPsp ();
  Nlinear = nrl_Linear ();
  psil = (double **) malloc (Nvalence * sizeof (double *));
  pspl = (double **) malloc (Nvalence * sizeof (double *));
  for (p = 0; p < Nvalence; ++p)
    {
      psil[p] = (double *) malloc (Nlinear * sizeof (double));
      pspl[p] = (double *) malloc (Nlinear * sizeof (double));
    }
  rl = (double *) malloc (Nlinear * sizeof (double));
  rhol = (double *) malloc (Nlinear * sizeof (double));


  /* Norm-conserving output */
  for (p = 0; p < Nvalence; ++p)
    {
      Log_to_Linear (r_psi_RelPsp (p), rl, psil[p]);
      Log_to_Linear_zero (V_RelPsp (p), rl, pspl[p]);

      /* normalize scattering state */
      if (fill_RelPsp (p) == 0.0)
	{
	  normalize_Linear (psil[p]);
	}
    }

  Log_to_Linear (rho_RelPsp (), rl, rhol);


  if (debug)
    {
      /* output pseudowavefunctions argv[1].psw */
      strcpy (name, name_Atom ());
      strcat (name, ".psw.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing pseudowavefunctions: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Nlinear; ++k)
	{
	  fprintf (fp, "%12.8lf", rl[k]);
	  for (p = 0; p < Nvalence; ++p)
	    fprintf (fp, " %12.8lf", psil[p][k]);
	  fprintf (fp, "\n");
	}
      fclose (fp);

      /* output pseudopotentials argv[1].psp.plt */
      strcpy (name, name_Atom ());
      strcat (name, ".psp.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing pseudopotentials: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Nlinear; ++k)
	{
	  fprintf (fp, "%12.8lf", rl[k]);
	  for (p = 0; p < Nvalence; ++p)
	    fprintf (fp, " %12.8lf", pspl[p][k]);
	  fprintf (fp, "\n");
	}
      fclose (fp);

      /* output pseudodensity infile.psd.plt */
      strcpy (name, name_Atom ());
      strcat (name, ".psd.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      fp = fopen (full_filename, "w+");
      for (k = 0; k < Nlinear; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rl[k], rhol[k]);
      fclose (fp);
    }

  /* output datafile to be used for Kleinman-Bylander input, argv[1].psp */
  if (print)
    {
      printf (" Creating datafile for Kleinman-Bylander input: %s\n",
	      outfile);
    }


  fp = fopen (outfile, "w+");
  fprintf (fp, "7 %s\n", name_Atom ());
  fprintf (fp, "%lf %lf %d   %d %d %lf\n", Zion_RelPsp (), Amass_Atom (),
	   lmax_RelPsp (), lmax_out, locp_out, rlocal_out);
  for (p = 0; p <= lmax_RelPsp (); ++p)
    fprintf (fp, "%lf ", rcut_RelPsp (p));
  fprintf (fp, "\n");
  fprintf (fp, "%d %lf\n", nrl_Linear (), drl_Linear ());
  fprintf (fp, "%s\n", comment_RelPsp ());

  if (print)
    {
      printf ("  + Appending pseudopotentials:    %s thru %s\n",
	      spd_Name (0), spd_Name (lmax_RelPsp ()));
    }

  nvh=Nvalence/2;
  for (k = 0; k < Nlinear; ++k)
    {
      fprintf (fp, "%12.8lf", rl[k]);
      for (p = 0; p < nvh; ++p)
	{
/**************************************************************
 *  Here we output the v average
 *     v_avg(l,r)=((l*v(l-1/2) + (l+1)*v(l+1/2))/(2*l+1)
 *************************************************************/
	  vx = ((p) * pspl[2 * p][k] + (p+1) * pspl[2 * p + 1][k])/(p+p+1);
	  fprintf (fp, " %12.8lf", vx);
	}
      fprintf (fp, "\n");
    }
  if (print)
    {
      printf ("  + Appending pseudowavefunctions: %s thru %s\n",
	      spd_Name (0), spd_Name (lmax_RelPsp ()));
    }


  for (k = 0; k < Nlinear; ++k)
    {
      fprintf (fp, "%12.8lf ", rl[k]);
      for (p = 0; p < Nvalence; ++p)
	fprintf (fp, " %12.8lf", psil[p][k]);
      fprintf (fp, "\n");
    }

  /* output spin-orbit datafile to be used for Kleinman-Bylander input, argv[1].dat */
  if (print)
    {
      printf
	(" Creating spin-orbit datafile for Kleinman-Bylander input: %s\n",
	 outfile);
    }
  if (print)
    {
      printf ("  + Appending Spin-Orbit pseudopotentials:    %s thru %s\n",
	      spd_Name (0), spd_Name (lmax_RelPsp ()));
    }
  for (k = 0; k < Nlinear; ++k)
    {
      fprintf (fp, "%12.8lf", rl[k]);
/**************************************************************
 *  Here we output the v spin orbit
 *     v_spin_orbit(l,r)=2*(v(l+1/2) - v(l-1/2))/(2*l+1)
 *
 *   so that V_spin_orbit|Psi>= -1/2(1+kappa)*v_spin_orbit|Psi>
 *   its confusing because L*S -> - (1+kappa)/2
 *************************************************************/
      for (p = 1; p < nvh;++p)
	{
	  vx = (pspl[2*p+1][k] - pspl[2*p][k]);
	  vx *= 2. / (2. * p + 1.);
	  fprintf (fp, " %15.8lf", vx);
	}
      fprintf (fp, "\n");
    }


  /* append semicore corrections */
  if (r_semicore_RelPsp () != 0.0)
    {
      if (print)
	{
	  printf ("  + Appending semicore density\n");
	}
      Log_to_Linear (rho_semicore_RelPsp (), rl, rhol);
      fprintf (fp, "%lf\n", r_semicore_RelPsp ());
      for (k = 0; k < Nlinear; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rl[k],
		 fabs (rhol[k] * over_fourpi));

      if (print)
	{
	  printf ("  + Appending semicore density gradient\n");
	}
      Log_to_Linear (drho_semicore_RelPsp (), rl, rhol);
      for (k = 0; k < Nlinear; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rl[k], (rhol[k] * over_fourpi));
    }
  fclose (fp);

   /******************************************************************/
   /******************* output AE information ************************/
   /******************************************************************/


  if (debug)
    {

      /* output all-electron wavefunctions */
      printf ("Outputing all-electron wavefunctions:");
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      for (p = 0; p < (Ncore_Atom () + Nvalence_Atom ()); ++p)
	{
	  sprintf (name, "%s.%1d%s%s", name_Atom (), n_Atom (p),
		   spd_Name (l_Atom (p)), spin_Name (p));
	  full_filename[0] = '\0';
	  strncpy (full_filename, sdir_name, m9);
	  full_filename[m9] = '\0';
	  strcat (full_filename, "/");
	  full_filename[m9 + 1] = '\0';
	  strcat (full_filename, name);

	  printf (" %s", full_filename);
	  fp = fopen (full_filename, "w+");
	  for (k = 0; k < Ngrid; ++k)
	    fprintf (fp, "%12.8lf %12.8lf\n", rgrid[k], r_psi_Atom (p)[k]);
	  fclose (fp);
	}
      printf ("\n");

      /* output density argv[1].dns */
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      strcpy (name, name_Atom ());
      strcat (name, ".dns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing atom density: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rgrid[k], rho_Atom ()[k]);
      fclose (fp);

      /* output core density argv[1].cdns */
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      strcpy (name, name_Atom ());
      strcat (name, ".cdns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing core density: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rgrid[k], rho_core_Atom ()[k]);
      fclose (fp);

      /* output core density gradient argv[1].cddns */
      vall = alloc_LogGrid ();
      Derivative_LogGrid (rho_core_Atom (), vall);
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      strcpy (name, name_Atom ());
      strcat (name, ".cddns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);
      printf ("Outputing core density gradient: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rgrid[k], vall[k]);
      fclose (fp);
      dealloc_LogGrid (vall);

      /* output semicore density infile.sdns */
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      strcpy (name, name_Atom ());
      strcat (name, ".sdns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing semicore density: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
	fprintf (fp, "%12.8lf %12.8lf\n", rgrid[k], rho_semicore_RelPsp ()[k]);
      fclose (fp);

      /* output semicore density gradient infile.sddns */
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      strcpy (name, name_Atom ());
      strcat (name, ".sddns.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing semicore density gradient: %s\n", full_filename);
      fp = fopen (full_filename, "w+");
      for (k = 0; k < Ngrid; ++k)
      {
        vx=drho_semicore_RelPsp()[k];
	fprintf (fp, "%12.8lf %12.8lf\n", rgrid[k], vx);
      }
      fclose (fp);

      /* output all-electron potential infile.pot */
      vall = Vall_Atom ();
      Ngrid = N_LogGrid ();
      rgrid = r_LogGrid ();
      strcpy (name, name_Atom ());
      strcat (name, ".pot.plt");
      full_filename[0] = '\0';
      strncpy (full_filename, sdir_name, m9);
      full_filename[m9] = '\0';
      strcat (full_filename, "/");
      full_filename[m9 + 1] = '\0';
      strcat (full_filename, name);

      printf ("Outputing all-electron potential(non-screened): %s\n",
	      full_filename);
      fp = fopen (full_filename, "w+");

      c = 0.0;
      for (k = 0; k < Ngrid; ++k)
	{
	  x = rgrid[k] / rcut_RelPsp (0);
	  x = pow (x, 3.5);
	  x = exp (-x);
	  y = 1.0 - x;
	  fprintf (fp, "%12.8lf %12.8lf %12.8lf\n", rgrid[k], vall[k],
		   y * vall[k] + c * x);
	}
      fclose (fp);
    }

  /* free malloc memory */
  free (infile);
  free (outfile);
  free (full_filename);
  for (p = 0; p < Nvalence; ++p)
    {
      free (psil[p]);
      free (pspl[p]);
    }
  free (psil);
  free (pspl);
  free (rl);
  free (rhol);
  end_Linear ();
  fflush (stdout);
  return;
}
