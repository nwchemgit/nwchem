c     Check density
c
      logical function chk_dens(rtdb, g_dens)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c      
      integer rtdb
      integer g_dens(2)
c
      double precision denmax,denmax1,denmax2
      double precision tol_rho
c
c     get user defined tol_rho
      tol_rho=1.D-10
      if (.not. rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho))
     &   tol_rho=1.D-10

      denmax1 = 0.0d0
      denmax2 = 0.0d0
      denmax=0.0d0
      call ga_maxelt(g_dens(1), denmax1)
      call ga_maxelt(g_dens(2), denmax2)
      denmax=max(denmax1,denmax2)
c
      chk_dens = .true.
      if (denmax.lt.tol_rho) chk_dens = .false.
c
      return
      end
c $Id$
