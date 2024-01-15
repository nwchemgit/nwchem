      subroutine x2c_scf_init(rtdb)
      implicit none
      integer rtdb
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "x2c_scf.fh"
c
c     check for additional x2c keywords and set
c     lx2c_keepU to T if we need it for any of the properties:
      lx2c_keepU = .false.
      lx2c_efg = .false.
      if (rtdb_get(rtdb, 'x2c:efg', mt_log, 1,
     &   lx2c_efg))    lx2c_keepU = .true.
      lx2c_hfc = .false.
      if (rtdb_get(rtdb, 'x2c:hfc', mt_log, 1,
     &   lx2c_hfc))    lx2c_keepU = .true.
      lx2c_moment = .false.
      if (rtdb_get(rtdb, 'x2c:moment', mt_log, 1,
     &   lx2c_moment)) lx2c_keepU = .true.
      lx2c_noatm = .false.
      if (rtdb_get(rtdb, 'x2c:noatm', mt_log, 1,
     &  lx2c_noatm)) continue

      if (lx2c_efg .and. ga_nodeid().eq.0) write (luout,*)
     &  'x2c: will calculate EFGs, too'
      
      return
      end
      subroutine x2c_scf_terminate(rtdb,
     i     ipol,
     g     g_rel,g_x2c_u,g_x2c_u2c,g_x2c_c2u)
      implicit none
      integer rtdb
      integer ipol
      integer g_rel(*)
      integer g_x2c_u(*),g_x2c_u2c,g_x2c_c2u
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "x2c_scf.fh"
c
c     deallocate x2c arrays if necessary
      if (.not. ga_destroy(g_rel(1))) call errquit(
     &     'dft_scf: ga_destroy failed g_rel',0, GA_ERR)
      if (ipol.gt.1) then
         if (.not. ga_destroy(g_rel(2))) call errquit(
     &        'dft_scf: ga_destroy failed g_rel',0, GA_ERR)
      end if
      if (lx2c_keepU) then
         if(.not.ga_destroy(g_x2c_u(1))) call errquit(
     &        'dft_scf: ga_destroy failed x2c_u(1)',1, GA_ERR)
         if(.not.ga_destroy(g_x2c_u(2))) call errquit(
     &        'dft_scf: ga_destroy failed x2c_u(2)',2, GA_ERR)
         if (lx2c_u2c) then
            if(.not.ga_destroy(g_x2c_u2c)) call errquit(
     &           'dft_scf: ga_destroy failed x2c_u2c',51, GA_ERR)
         end if
         if (lx2c_c2u) then
            if(.not.ga_destroy(g_x2c_c2u)) call errquit(
     &           'dft_scf: ga_destroy failed x2c_c2u',54, GA_ERR)
         end if
      endif
      return
      end
      subroutine x2c_prop(rtdb,
     &     nexc, g_dens, g_dens_at, g_x2c_u, g_x2c_u2c, g_x2c_c2u)
      implicit none
      integer rtdb
      integer nexc, g_dens, g_dens_at, g_x2c_u, g_x2c_u2c, g_x2c_c2u
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "x2c_scf.fh"
c
      if (lx2c_efg) then
c     x2c EFG 
        if (ga_nodeid().eq.0)
     &    write (luout,*) 'calling calc_x2c_efg_scalar' 
        call calc_x2c_efg_scalar(rtdb,
     &    nexc, g_dens, g_dens_at, g_x2c_u)
      end if
      
      if (lx2c_hfc) then
c     x2c hyperfine coupling
         if (ga_nodeid().eq.0)
     &        write (luout,*) 'calling calc_x2c_hfc_scalar' 
         call calc_x2c_hfc_scalar(rtdb, lx2c_u2c,
     &     nexc, g_dens, g_dens_at, g_x2c_u, g_x2c_u2c, g_x2c_c2u)

      end if       
      
      if (lx2c_moment) then
c     x2c dipole moment
         if (ga_nodeid().eq.0)
     &        write (luout,*) 'calling calc_x2c_moment_scalar' 
         call calc_x2c_moment_scalar(rtdb,
     &        nexc, g_dens, g_dens_at, g_x2c_u, 'dipole', 3)
      end if
      return
      end
      logical function x2c_noatm()
      implicit none
#include "x2c_scf.fh"      
      x2c_noatm=lx2c_noatm
      return
      end
