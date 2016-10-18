c
c     Coulomb part
c
      subroutine getv_coul(rtdb, ecoul, g_densso, g_fockso, g_v)
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
      double precision ecoul
c
      double precision jfac(4),kfac(4)
c
      integer me,type
      double precision tol2e
      integer g_tmp(2),g_dens(2)
      integer nao
c
      character*32 pname
c
      logical chk_dens
      external chk_dens
c
      me=ga_nodeid()
      tol2e=10.d0**(-itol2e)
c
      pname = 'getv_coul: '
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
      call ga_dens_aabbr(g_dens, g_densso, nao)
c
c     Form a total density matrix
      call ga_dadd(1.d0, g_dens(1), 1.d0, g_dens(2), g_dens(1))
      call ga_copy(g_dens(1),g_dens(2))
c
c     Check density
      if (.not. chk_dens(rtdb, g_dens)) goto 1000
c
c     Calculate Coulomb
      jfac(1) = 1.0d0
      jfac(2) = 1.0d0
      kfac(1) = 0.0d0
      kfac(2) = 0.0d0
      call ga_zero(g_tmp(1))
      call ga_zero(g_tmp(2))
      call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_tmp, .false.)
c
c     Accumulate contribution 
      call ga_zero(g_v(1))
      call ga_zero(g_v(2))
      call ga_dadd(1.0d0, g_v(1), 1.0d0, g_tmp(1), g_v(1))
      call ga_dadd(1.0d0, g_v(2), 1.0d0, g_tmp(2), g_v(2))
      call ga_copy(g_v(1), g_v(2))
      ecoul = 0.5d0*ga_ddot(g_dens(1),g_v(1))
c
c     == patch in the Coulomb into the fock matrix ==
      call ga_sync()
      call ga_dadd_patch( 1.d0, g_fockso(1), 1, nao,
     &     1, nao,
     &     1.0d0, g_v(1),  1, nao,
     &     1, nao,
     &     g_fockso(1), 1, nao,
     &     1, nao)
      call ga_dadd_patch( 1.d0, g_fockso(1), 1+nao, 2*nao,
     &     1+nao, 2*nao,
     &     1.0d0, g_v(2),  1, nao,
     &     1, nao,
     &     g_fockso(1), 1+nao, 2*nao,
     &     1+nao, 2*nao)
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
