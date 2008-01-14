      subroutine tce_aod1(rtdb,g_aod1,dipole)
c
c $Id: tce_aod1.F,v 1.3 2008-01-14 15:20:33 jhammond Exp $
c
c     Returns a GA handle for AO Fock matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb               ! Run-time database
      integer g_aod1(3)          ! GA handle for AO dipole matrices
      integer g_dens(2)          ! GA handle for density matrices
      double precision dipole(3) ! SCF dipole moments
      integer i
c
c     ===============
c     Allocate g_dens
c     ===============
c
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO density',
     1    -1,-1,g_dens(i))) call errquit('tce_ao1e: GA problem',2,
     2    GA_ERR)
      enddo
c
c     =====================
c     Make density matrices
c     =====================
c
      do i=1,ipol
        call ga_dgemm('N','T',nbf,nbf,nocc(i),1.0d0,
     1    g_movecs(i),g_movecs(i),0.0d0,g_dens(i))
        if (ipol.eq.1) call ga_dscal(g_dens(i),2.0d0)
        if (util_print('ao1e',
     1    print_debug)) call ga_print(g_dens(i))
      enddo
c
c     ===============
c     Allocate g_aod1
c     ===============
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole X',
     1  -1,-1,g_aod1(1))) call errquit('tce_aod1: GA problem',0,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole Y',
     1  -1,-1,g_aod1(2))) call errquit('tce_aod1: GA problem',1,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole Z',
     1  -1,-1,g_aod1(3))) call errquit('tce_aod1: GA problem',2,
     2  GA_ERR)
c
c     ====================
c     Make dipole matrices
c     ====================
c
      call int_dip_ga(ao_bas_han,ao_bas_han,
     1  g_aod1(1),g_aod1(2),g_aod1(3))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(1))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(2))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(3))
c
c     ===================
c     Make dipole moments
c     ===================
c
c     Note: do not attempt to compute these in MO basis;
c     frozen core approximation can give you wrong results.
c
      do i=1,ipol
        dipole(1) = dipole(1) + ga_ddot(g_dens(i),g_aod1(1))
        dipole(2) = dipole(2) + ga_ddot(g_dens(i),g_aod1(2))
        dipole(3) = dipole(3) + ga_ddot(g_dens(i),g_aod1(3))
      enddo
      do i=1,ipol
        if (.not.ga_destroy(g_dens(i)))
     1  call errquit('tce_aod1: failed to destroy g_dens',0,GA_ERR)
      enddo
      return
      end
