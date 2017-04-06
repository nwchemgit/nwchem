c
c     Wrapper for Coulomb part using charge density fitting
c
      subroutine getv_coul_cd(rtdb, ecoul, g_densso, g_fockso, g_v,
     &     k_3cERI, nTotEl, n_batch, n3c_int,
     &     k_3cwhat, n3c_dbl, IOLGC, natoms, iVcoul_opt)
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
      integer nTotEl, n_batch, n3c_int
      integer n3c_dbl, IOLGC, natoms, iVcoul_opt
      double precision k_3cERI(n3c_dbl)
      integer k_3cwhat(n3c_int)
c
      integer me,type
      double precision tol2e
      integer g_tmp(2),g_dens(2)
      integer nao
c
      integer lcd_coef, icd_coef
c
      character*32 pname
c
      logical chk_dens
      external chk_dens
c
      me=ga_nodeid()
      tol2e=10.d0**(-itol2e)
c
      pname = 'getv_coul_cd: '
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
c     Calculate Coulomb using CD fitting
c
c     Fit the electron charge density.
c
      icd_coef = 0
      if (.not.ma_push_get(MT_Dbl,nbf_cd,'cd_coef',lcd_coef, icd_coef))
     &   call errquit(pname//'cannot allocate cd_coef',0, MA_ERR)
c
      call dft_fitcd(1,Dbl_MB(icd_coef),k_3cERI,Ecoul,
     &               g_dens, nTotEl, n_batch, n3c_int,
     &               k_3cwhat, n3c_dbl, IOLGC,
     &               natoms, .false., 0d0, .false.)
c
      call dft_getvc(Dbl_MB(icd_coef),k_3cERI,Ecoul,
     &               g_tmp, iVcoul_opt, n_batch,
     &               n3c_int, k_3cwhat, n3c_dbl,
     &               IOLGC, .false., 1)
c
c     == patch in the Coulomb into the fock matrix ==
      call ga_sync()
      call ga_dadd_patch( 1.d0, g_fockso(1), 1, nao,
     &     1, nao,
     &     1.0d0, g_tmp,  1, nao,
     &     1, nao,
     &     g_fockso(1), 1, nao,
     &     1, nao)
      call ga_dadd_patch( 1.d0, g_fockso(1), 1+nao, 2*nao,
     &     1+nao, 2*nao,
     &     1.0d0, g_tmp,  1, nao,
     &     1, nao,
     &     g_fockso(1), 1+nao, 2*nao,
     &     1+nao, 2*nao)
c
      if (.not.ma_pop_stack(lcd_coef))
     &   call errquit(pname//'cannot pop stack cd_coef',0, MA_ERR)
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
c $Id: getv_coul.F 28689 2016-10-18 15:10:10Z niri $
