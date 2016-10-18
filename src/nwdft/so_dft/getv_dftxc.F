c
c     Exchange-correlation part
c
      subroutine getv_dftxc(rtdb, nexc, edftxc, rho_n, g_densso, 
     &              g_fockso, g_v, which)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "oep.fh"
#include "dftpara.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
c
      integer rtdb, g_densso(2), g_fockso(2), g_v(4)
      double precision edftxc(2)
      character*4 which
c
      integer me,type
      double precision tol2e,denmax
      double precision rho_n
      integer g_tmp(2),g_dens(2)
      integer nao, nexc
      double precision tmat
c
      character*32 pname
c
      logical xc_gotxc, chk_dens
      external xc_gotxc, chk_dens
c
      me=ga_nodeid()
      tol2e=10.d0**(-itol2e)
      pname = 'getv_dftxc :'
c
c     check if any dft xc is present
      if (.not.xc_gotxc()) return
c
      call ga_sync
c
c     Query the full density matrix
      call ga_inquire(g_densso(1),type,nao,nao)
      nao=nao/2
c
c     Allocate some scratch space
      if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_tmp(1)))
     &        call errquit(pname//'error creating ga', 0, GA_ERR)
      if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_tmp(2)))
     &        call errquit(pname//'error creating ga', 0, GA_ERR)
      if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(1)))
     &        call errquit(pname//'error creating ga', 0, GA_ERR)
      if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(2)))
     &        call errquit(pname//'error creating ga', 0, GA_ERR)
c
c     Get real part of the density matrix 
c     g_dens(1) = alpha, g_dens(2)= beta
      call ga_dens_aabbr(g_dens, g_densso, nao)
c
c     Check density
      if (.not. chk_dens(rtdb, g_dens)) goto 1000
c
      if (which .eq. 'summ') then       ! summation density
c        Form a total density matrix
         call ga_dadd(0.5d0, g_dens(2), 0.5d0, g_dens(1), g_dens(1))
         call ga_copy(g_dens(1),g_dens(2))
      else if (which .eq. 'diff') then  ! difference density
c        Form a difference density matrix
         call ga_dadd(0.5d0, g_dens(2),-0.5d0, g_dens(1), g_dens(1))
         call ga_copy(g_dens(1),g_dens(2))
      else
c        Treat the densities g_dens(1) and g_dens(2) as is
      end if
c
      edftxc(1) = 0.d0
      edftxc(2) = 0.d0
      call ga_zero(g_tmp(1))
      call ga_zero(g_tmp(2))
      rho_n = 0.d0
      tmat = 0.d0
      call grid_quadv0(rtdb, g_dens, g_tmp, nexc, rho_n, edftxc, tmat)
c
c     Accumulate contribution 
      call ga_zero(g_v(1))
      call ga_zero(g_v(2))
      call ga_dadd(1.0d0, g_v(1), 1.0d0, g_tmp(1), g_v(1))
      call ga_dadd(1.0d0, g_v(2), 1.0d0, g_tmp(2), g_v(2))
c
c     == patch in the exchange-correlation potential into the fock matrix ==
      call ga_sync()
      call ga_dadd_patch( 1.d0, g_fockso(1), 1, nao, 1, nao,
     &     1.0d0, g_v(1),  1, nao, 1, nao,
     &     g_fockso(1), 1, nao, 1, nao)
      call ga_dadd_patch( 1.d0, g_fockso(1), 1+nao, 2*nao,
     &     1+nao, 2*nao, 1.0d0, g_v(2),  1, nao,
     &     1, nao, g_fockso(1), 1+nao, 2*nao, 1+nao, 2*nao)
c
 1000 continue
c
c     Deallocate scratch
      if (.not. ga_destroy(g_tmp(1))) call errquit
     $          (pname//'ga corrupt?',0, GA_ERR)
      if (.not. ga_destroy(g_tmp(2))) call errquit
     $          (pname//'ga corrupt?',1, GA_ERR)
      if (.not. ga_destroy(g_dens(1))) call errquit
     $          (pname//'ga corrupt?',0, GA_ERR)
      if (.not. ga_destroy(g_dens(2))) call errquit
     $          (pname//'ga corrupt?',1, GA_ERR)
c
      return
      end
c $Id$
