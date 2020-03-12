c
      subroutine hess_energyinit(rtdb)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
c  This routine does a few checks to make sure that the energy will
c  be calculated without symmetry and with reasonable thresholds so
c  the rest of the Hessian calculation will work properly.  Note that
c  these values should be reset after the calculation is done.
c
      integer rtdb
      integer geomnew, filesize, memsize
      double precision ithresh,tol_rho,tol_rho_rd,e_conv,
     .     job_grid_acc
      integer itol2e
      logical xc_got2nd      
      external xc_got2nd
      logical grid_clinit,gstatus
      external grid_clinit
      character*15 gridtype
      character*24 rtdb_errmsg
      logical int_disablesim
      external int_disablesim
c     
c  Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hess_energyinit: theory not specified',555,
     &       INPUT_ERR)
      if (theory.eq.'hyb') theory = 'dft'
c
c  Check to see if the current geometry has symmetry.  If so, create
c  a new geometry without symmetry and put onto the rtdb.  When symmetry
c  is implemented into the code, this can go away.
c
      if (.not. geom_create(geomold, 'geometry'))
     *  call errquit('hess_init: geom_create problem', 555, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geomold, 'geometry'))
     *  call errquit('hess_init: no geometry ', 555, GEOM_ERR)
      if (.not. geom_makec1(geomold, geomnew))
     $  call errquit('hess_init: failure to decrease symmetry',555,
     &       GEOM_ERR)
      if (.not.geom_rtdb_store(rtdb,geomnew,'geometry'))
     &      call errquit(' hess_init: geom_rtdb_store failed',555,
     &       RTDB_ERR)
      if (.not.geom_destroy(geomnew))
     &      call errquit(' hess_init: geom_destroy failed',555,
     &       GEOM_ERR)
c
c  Need to check to make sure that the convergence tolerances for
c  the wavefunction are OK.  Use the user defined one if available,
c  otherwise check to make sure that the threshold gets set to a
c  "reasonable" value.
c
      if (rtdb_get(rtdb, 'hessian:thresh', MT_DBL,1, ithresh)) then
          if (.not. rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (.not. rtdb_put(rtdb, 'scf:thresh', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        if (theory.eq.'scf') then
        else
          if (.not. rtdb_get(rtdb, 'dft:g_conv', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (.not. rtdb_put(rtdb, 'dft:g_conv', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        endif
      else
        ithresh = 1.0d-6
          if (.not. rtdb_get(rtdb, 'scf:thresh', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (tthresh.lt.ithresh) ithresh = tthresh
          if (.not. rtdb_put(rtdb, 'scf:thresh', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        if (theory.eq.'scf') then
        else
        ithresh = 1.0d-6
          if (.not. rtdb_get(rtdb, 'dft:g_conv', MT_DBL, 1, tthresh))
     $       tthresh = 1.0d-4
          if (tthresh.lt.ithresh) ithresh = tthresh
          if (.not. rtdb_put(rtdb, 'dft:g_conv', MT_DBL, 1, ithresh))
     $      call errquit('hess_energyinit: rtdb_put failed',ithresh,
     &       RTDB_ERR)
        endif
      endif
      if(theory.eq.'dft') then
         if(.not.xc_got2nd()) call errquit(
     A        'analytic 2nds not ready for these XC functionals',0,
     &       CAPMIS_ERR)
c
c     set grid to fine or leave the user to set it to xfine
c
         if (rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype)) then
            if(gridtype.eq.'coarse'.or.gridtype.eq.'xcoarse'.
     $           or.gridtype.eq.'medium') gridtype='fine'
         else
            gridtype='fine'
         endif
 
         if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     $        call errquit('hessinit:rtdbcput failed',0, RTDB_ERR)
c
c     tighten itol2e to avoid cphf  residual woes to lt 1d-12
c
      if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, 
     &   itol2e)) then
          itol2e =12
      else
         itol2e = max(itol2e,12)
      endif
      if (.not. rtdb_put(rtdb, 'dft:itol2e', mt_int, 1, 
     &   itol2e))call errquit('hessint:rtdbput failed',1,
     &       RTDB_ERR)
c     set grid_acc consistent with fine (at least), ie 1d-7
c     make it even tighter for cases with diffuse basis sets and
c     low frequencies ... > 1d-14. (Edo 04/08/2015)
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc)) then
          job_grid_acc = 1d-14
      else
         job_grid_acc = min(job_grid_acc,1d-14)
      endif
      if (.not. rtdb_put(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))call errquit('hessint:rtdbput failed',1,
     &       RTDB_ERR)
c
c     weight derivatives bit turned off
c
      if (.not. rtdb_put(rtdb, 'dft:weight derivatives', mt_log, 1, 
     & .false.))call errquit('hessint:rtdbput failed',2, RTDB_ERR)
      
c
c     set tol_rho to at least 1d-11 (it is 1d-10 otherwise)
c
         tol_rho=1d-11
         if (rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho_rd))
     .        tol_rho=min(tol_rho,tol_rho_rd)
         if (.not.rtdb_put(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho)) 
     .        call errquit(' hessinit: rtdbput failed ',11,0)
c     
c     kill symmetry for XC grid build
c     
         if (.not. rtdb_put(rtdb, 'dft:gridsym',
     $        MT_LOG, 1, .false.)) call errquit(
     $        ' task_hessian: rtdb_put failed',0, RTDB_ERR)
c
c     get xfac for doing xc_gotxc right on restart
c     but ignore the status since there are not xfac 
c     on rtdb when we are in start mode
c
         gstatus=grid_clinit(rtdb)
         call grid_cleanup(.false.)
c
c     set deltaE at least to 1d-7 
c
      if (.not.rtdb_get(rtdb,'dft:e_conv',mt_dbl,1,e_conv))then
        e_conv=1d-7
      else
        e_conv=min(e_conv,1d-7)
      endif
      if (.not. rtdb_put(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))
     &   call errquit('hess_init: rtdb_put failed', 100,
     &       RTDB_ERR)

       call dft_inpanae(rtdb)
      endif
#if 0
c
c     make the Fock build go direct
c
      filesize = -1
      memsize  = -1
      if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, filesize))
     $  call errquit('hess_energyinit: failed to store filesize',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $  call errquit('hess_energyinit: failed to store filesize',0,
     &       RTDB_ERR)
#endif
c
c     skip Simint integrals since they do not have 2nds yet
c
      if(.not.int_disablesim(rtdb,rtdb_errmsg))
     &        call errquit(rtdb_errmsg,0, RTDB_ERR)
      return
      end
c $Id$
