c
      subroutine shell_fock_build_cam_single(geom, basis, 
     &     nder, ndens, nbf_ao, jfac, kfac, tol2e, 
     &     g_dens, g_fock, g_work, ipol, 
     &     doexch, docoul, doxc, osym, oprint)
c
c     == wrapper routine to handle exchange attenuated functionals ==
c     == for TDDFT algorithms 2 & 4 ==
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
      integer nder, ndens, ndu, nfu, nbf_ao
      integer g_dens, g_fock, g_work, ipol
      double precision jfac, kfac, tol2e
      logical doexch, docoul, doxc, osym, oprint
c
c     == local variables ==
      double precision jf(4)   ! Prefactor for HF Coulomb as input for shell_fock_build
      double precision kf(4)   ! Prefactor for HF exchange as input for shell_fock_build
c
c     == preliminaries ==
      if ((.not.doexch).and.(.not.docoul).and.(.not.doxc)) return
c
c     == build matrices for polarized & unpolarized calculations ==
      if (ipol.eq.1) then  ! unpolarized calculations
c
c       == exchange part ==
        if (doexch) then
          call case_setflags(.true.)    ! set attenuation for exchange
          jf(1)= 0.0
          kf(1)=-kfac*0.5d0
          call ga_zero(g_work)
          call shell_fock_buildjk(geom,basis,nder,ndens,
     &       jf,kf,tol2e,g_dens,g_work,osym)
          call ga_add(1d0, g_work, 1.0d0, g_fock, g_fock)
        end if  ! doexch
c
c       == coulomb part ==
        if (docoul) then
          call case_setflags(.false.)  ! unset attenuation for coulomb
          jf(1)= jfac
          kf(1)= 0.0
          call ga_zero(g_work)
          call shell_fock_buildjk(geom,basis,nder,ndens,
     &       jf,kf,tol2e,g_dens,g_work,osym)
          call ga_add(1d0, g_work, 1.0d0, g_fock, g_fock)
        end if ! docoul
c
c       == exchange-correlation part ==
        if (doxc) then
          call case_setflags(.false.)  ! unset attenuation for exchange-correlation
          jf(1)= 0.0
          kf(1)= 0.0
          call ga_zero(g_work)
          call shell_fock_buildxc(geom,basis,nder,ndens,
     &       jf,kf,tol2e,g_dens,g_work,osym)
          call ga_add(1d0, g_work, 1.0d0, g_fock, g_fock)
        end if ! doxc
c
      else  ! polarized calculations
c
        ndu = 2
        nfu = 2
c
c       == exchange part ==
        if (doexch) then
          call case_setflags(.true.)  ! set attenuation for exchange
          jf(1)=0.0d0
          kf(1)=kfac
          jf(2)=0.0d0
          kf(2)=kfac
          call ga_zero(g_work)
          call uhf_fock_setup(g_dens,g_work,jf,kf,ndu,
     &       nfu,nbf_ao,oprint)
          call shell_fock_buildjk(geom,basis,nder,ndu,
     &      jf,kf,tol2e,g_dens,g_work,osym)
          call uhf_fock_finish(g_dens,g_work,jf,kf,ndu,
     &       nfu,nbf_ao,oprint)
          call ga_add(1d0, g_work, 1.0d0, g_fock, g_fock)
        end if  ! doexch
c
c       == coulomb part ==
        if (docoul) then
          call case_setflags(.false.) ! unset attenuation for coulomb
          jf(1)=jfac
          kf(1)=0.0
          jf(2)=jfac
          kf(2)=0.0
          call ga_zero(g_work)
          call uhf_fock_setup(g_dens,g_work,jf,kf,ndu,
     &       nfu,nbf_ao,oprint)
          call shell_fock_buildjk(geom,basis,nder,ndu,
     &      jf,kf,tol2e,g_dens,g_work,osym)
          call uhf_fock_finish(g_dens,g_work,jf,kf,ndu,
     &       nfu,nbf_ao,oprint)
          call ga_add(1d0, g_work, 1.0d0, g_fock, g_fock)
        end if ! docoul
c
c       == exchange-correlation part ==
        if (doxc) then
          call case_setflags(.false.) ! unset attenuation for exchange-correlation
          jf(1)=0.0
          kf(1)=0.0
          jf(2)=0.0
          kf(2)=0.0
          call ga_zero(g_work)
          call uhf_fock_setup(g_dens,g_work,jf,kf,ndu,
     &       nfu,nbf_ao,oprint)
          call shell_fock_buildxc(geom,basis,nder,ndu,
     &      jf,kf,tol2e,g_dens,g_work,osym)
          call uhf_fock_finish(g_dens,g_work,jf,kf,ndu,
     &       nfu,nbf_ao,oprint)
          call ga_add(1d0, g_work, 1.0d0, g_fock, g_fock)
        end if ! doxc
c
      end if  ! ipol 
c
      return
      end
