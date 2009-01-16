c
      subroutine shell_fock_build_cam_multiple(geom, basis, 
     &     nder, nvectors, nbf_ao,
     &     jfac, kfac, tol2e, 
     &     g_dens, g_fock, g_temp, ipol, 
     &     doexch, docoul, doxc, osym, oprint)
c
c     == wrapper routine to handle exchange attenuated functionals ==
c     == for TDDFT algorithms 1 & 3 ==
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "case.fh"
c
      integer geom, basis
      integer nder, nvectors, ndu, nfu, nbf_ao
      integer g_dens, g_fock, g_temp, ipol
      double precision jfac, kfac, tol2e
      logical doexch, docoul, doxc, osym, oprint
c
c     == local variables ==
      integer n
      integer l_jfloc,k_jfloc        ! Prefactor for HF Coulomb as input for shell_fock_build
      integer l_kfloc,k_kfloc        ! Prefactor for HF exchange as input for shell_fock_build
c
c     == preliminaries ==
      if ((.not.doexch).and.(.not.docoul).and.(.not.doxc)) return
c
      if (ipol.eq.1) then
      if (.not.ma_push_get(mt_dbl,nvectors,'J factors',
     &  l_jfloc,k_jfloc)) 
     &call errquit('shell_fock_build_att_alg13: failed to allocate jf',
     &   0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nvectors,'K factors',
     &  l_kfloc,k_kfloc)) 
     &call errquit('shell_fock_build_att_alg13: failed to allocate kf',
     &   0, MA_ERR)
      else
      if (.not.ma_push_get(mt_dbl,4*nvectors,'J factors',
     & l_jfloc,k_jfloc)) 
     &call errquit('shell_fock_build_att_alg13: failed to allocate jf',
     &   0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,4*nvectors,'K factors',
     & l_kfloc,k_kfloc)) 
     &call errquit('shell_fock_build_att_alg13: failed to allocate kf',
     &   0, MA_ERR)
      endif
c
c     == build matrices for polarized & unpolarized calculations ==
      if (ipol.eq.1) then  ! unpolarized calculations
c
c       == exchange part ==
        if (doexch) then
          call case_setflags(.true.)    ! set attenuation for exchange
          do n=1,nvectors
             dbl_mb(k_jfloc+n-1)=0.d0
             dbl_mb(k_kfloc+n-1)=-kfac*0.5d0
          enddo
          call ga_zero(g_temp)
          call shell_fock_buildjk(geom,basis,nder,nvectors,
     &       dbl_mb(k_jfloc),dbl_mb(k_kfloc),tol2e,g_dens,g_temp,osym)
          call ga_add(1d0, g_temp, 1.0d0, g_fock, g_fock)
        end if  ! doexch
c
c       == coulomb part ==
        if (docoul) then
          call case_setflags(.false.)  ! unset attenuation for coulomb
          do n=1,nvectors
             dbl_mb(k_jfloc+n-1)=jfac
             dbl_mb(k_kfloc+n-1)=0.d0
          enddo
          call ga_zero(g_temp)
          call shell_fock_buildjk(geom,basis,nder,nvectors,
     &       dbl_mb(k_jfloc),dbl_mb(k_kfloc),tol2e,g_dens,g_temp,osym)
          call ga_add(1d0, g_temp, 1.0d0, g_fock, g_fock)
        end if ! docoul
c
c       == exchange-correlation part ==
        if (doxc) then
          call case_setflags(.false.)  ! unset attenuation for exchange-correlation
          do n=1,nvectors
             dbl_mb(k_jfloc+n-1)=0.d0
             dbl_mb(k_kfloc+n-1)=0.d0
          enddo
          call ga_zero(g_temp)
          call shell_fock_buildxc(geom,basis,nder,nvectors,
     &       dbl_mb(k_jfloc),dbl_mb(k_kfloc),tol2e,g_dens,g_temp,osym)
          call ga_add(1d0, g_temp, 1.0d0, g_fock, g_fock)
        end if ! doxc
c
      else  ! polarized calculations
c
        ndu=2*nvectors
        nfu=2*nvectors
c
c       == exchange part ==
        if (doexch) then
         call case_setflags(.true.)  ! set attenuation for exchange
         do n=1,2*nvectors
            dbl_mb(k_jfloc+n-1)=0.0d0
            dbl_mb(k_kfloc+n-1)=kfac
         enddo
         call ga_zero(g_temp)
         call uhf_fock_setup(g_dens,g_temp,
     &     dbl_mb(k_jfloc),dbl_mb(k_kfloc),ndu,nfu,nbf_ao,oprint)
         call shell_fock_buildjk(geom,basis,nder,ndu,
     &     dbl_mb(k_jfloc),dbl_mb(k_kfloc),tol2e,g_dens,g_temp,osym)
         call uhf_fock_finish(g_dens,g_temp,
     &     dbl_mb(k_jfloc),dbl_mb(k_kfloc),ndu,nfu,nbf_ao,oprint)
         call ga_add(1d0, g_temp, 1.0d0, g_fock, g_fock)
        end if  ! doexch
c
c       == coulomb part ==
        if (docoul) then
          call case_setflags(.false.) ! unset attenuation for coulomb
          do n=1,2*nvectors
             dbl_mb(k_jfloc+n-1)=jfac
             dbl_mb(k_kfloc+n-1)=0.d0
          enddo
          call ga_zero(g_temp)
          call uhf_fock_setup(g_dens,g_temp,
     &      dbl_mb(k_jfloc),dbl_mb(k_kfloc),ndu,nfu,nbf_ao,oprint)
          call shell_fock_buildjk(geom,basis,nder,ndu,
     &      dbl_mb(k_jfloc),dbl_mb(k_kfloc),tol2e,g_dens,g_temp,osym)
          call uhf_fock_finish(g_dens,g_temp,
     &      dbl_mb(k_jfloc),dbl_mb(k_kfloc),ndu,nfu,nbf_ao,oprint)
          call ga_add(1d0, g_temp, 1.0d0, g_fock, g_fock)
        end if  ! docoul
c
c       == exchange-correlation part ==
        if (doxc) then
          call case_setflags(.false.) ! unset attenuation for exchange-correlation
          do n=1,2*nvectors
             dbl_mb(k_jfloc+n-1)=0.d0
             dbl_mb(k_kfloc+n-1)=0.d0
          enddo
          call ga_zero(g_temp)
          call uhf_fock_setup(g_dens,g_temp,
     &      dbl_mb(k_jfloc),dbl_mb(k_kfloc),ndu,nfu,nbf_ao,oprint)
          call shell_fock_buildxc(geom,basis,nder,ndu,
     &      dbl_mb(k_jfloc),dbl_mb(k_kfloc),tol2e,g_dens,g_temp,osym)
          call uhf_fock_finish(g_dens,g_temp,
     &      dbl_mb(k_jfloc),dbl_mb(k_kfloc),ndu,nfu,nbf_ao,oprint)
          call ga_add(1d0, g_temp, 1.0d0, g_fock, g_fock)
        end if  ! doxc
c
      end if  ! ipol 
c
c     == deallocate ==
      if (.not.ma_pop_stack(l_kfloc)) 
     & call errquit('shell_fock_build_att_alg13: failed to deallocate 
     & kfloc',0, MA_ERR)
      if (.not.ma_pop_stack(l_jfloc)) 
     & call errquit('shell_fock_build_att_alg13: failed to deallocate 
     & jfloc',0, MA_ERR)
c
      return
      end
