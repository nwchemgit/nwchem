      Subroutine xc_getv(rtdb, Exc, ecoul,nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, g_wght, g_xyz,g_nq, 
     &                   wght_GA, rho_n, rdens_atom,
     &                   cetobfr, natoms)
c
      implicit none
#include "errquit.fh"
c      
      integer nExc
      integer iVxc_opt
      integer g_xcinv, g_dens(2), g_vxc(4), g_wght, g_xyz,g_nq
      integer natoms
      logical IOLGC, wght_GA
      integer rtdb
c
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
#include "case.fh"
c
      integer cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
      integer g_jk(4), g_d(4)
      logical havehfxc
c     
      integer  ga_create_atom_blocked
      logical xc_gotxc
      external ga_create_atom_blocked,xc_gotxc
c
c--> XC Energy
c
      double precision Exc(2)
      double precision ecoul ! [output]
c
c This driver routine solves for the XC energy and potential (Vxc) via 
c numerical quadrature methods. The results are obtained either by direct 
c numerical integration or by means of a LSQ fit of the Vxc to a set of 
c Gaussian functions. This fitted function can be used to evaluate Vxc
c via a summation of a series of 3-center overlap integrals (3OIs). The
c algorithms are formulated in terms of matrix products. See subsequent
c subroutines for further explanation.
c 
c              XC Energy and Potential Index Key, Vxc(pq,i)
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      
c              2        1   |   Vxc^up   Vxc^dw   
c              1        2   |   Vxc      
c              2        2   |   Vxc^up   Vxc^dw   
c
c           nTcols = ipol 
c
      integer me,nTrows,nTcols
      integer lTmat,iTmat,g_truevxc(2)
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time, 
     ,     oprint_oep
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e
      integer g_tmp(2)
c
c     timings
c
      double precision time1_2e,time2_2e
      double precision time1_xc,time2_xc
c
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
      oprint_oep = util_print('oep', print_high)
      Exc(1)=0.d0
      Exc(2)=0.d0
c
      me=ga_nodeid()
      havehfxc=abs(xfac(1)).gt.1d-8
c      nproc=ga_nnodes()
c
      if (oprint_intermediate_xc)then
         write(luout,*)' rtdb, Exc, nExc, iVxc_opt, g_xcinv: ',
     &               rtdb, Exc, nExc, iVxc_opt, g_xcinv
         write(luout,*)'g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA:',
     &               g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA
         write(luout,*)' Fock XC matrix entering xc_getv: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
         call ga_print(g_dens(1))
         if(ipol.eq.2)call ga_print(g_dens(2))
      endif
c
      if(util_print('dft timings', print_high))
     &      time1_2e=util_cpusec()   ! start 2e build time
      if (havehfxc .or. (.not. CDFIT))then
c       
c        Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
         tol2e=10.d0**(-itol2e)
         call ga_sync
         if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
         if (ipol.eq.1) then
            if (.not. CDFIT) then
              if (.not.cam_exch) then  ! for regular calculations
c
c               set up prefactors
                kfac(1) = -0.5d0*xfac(1)
                jfac(1) = 0.0d0
                jfac(2) = 1.0d0
                kfac(2) = 0.0d0
c
c               get some work space
                g_vxc(2) = ga_create_atom_blocked(geom,ao_bas_han,'jk')
c
c               calculate the exchange and coulomb parts
                call ga_zero(g_vxc(2))
                g_dens(2)=g_dens(1)
                call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_dens, g_vxc, .false.)
                Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
                ecoul = 0.5d0*ga_ddot(g_dens(1),g_vxc(2)) 
                call ga_dadd(1d0,g_vxc(1),1d0,g_vxc(2),g_vxc(1))
                if (.not. ga_destroy(g_vxc(2))) call errquit
     $             ('xc_getv: ga corrupt?',0, GA_ERR)
              else  ! CAM calculations
c
c               get some work space
                g_tmp(1)=ga_create_atom_blocked(geom,ao_bas_han,'work')
                call ga_zero(g_tmp(1))
c
                g_tmp(2)=ga_create_atom_blocked(geom,ao_bas_han,'work')
                call ga_zero(g_tmp(2))
c
c               set up prefactors for exchange
                kfac(1) = -0.5d0*xfac(1)
                jfac(1) = 0.0d0
                kfac(2) = 0.0d0
                jfac(2) = 0.0d0
                g_dens(2)=g_dens(1)
                call case_setflags(.true.)
                call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_dens, g_tmp, .false.)
                Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_tmp(1))
                call ga_dadd(1d0,g_vxc(1),1d0,g_tmp(1),g_vxc(1))
                call case_setflags(.false.)
c
c               calculate the full Coulomb
                call ga_zero(g_tmp(1))
                call ga_zero(g_tmp(2))
c
c               set up prefactors for coulomb
                kfac(1) = 0.0d0
                jfac(1) = 1.0d0
                kfac(2) = 0.0d0
                jfac(2) = 0.0d0
                g_dens(2)=g_dens(1)
                call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_dens, g_tmp, .false.)
                ecoul = 0.5d0*ga_ddot(g_dens(1),g_tmp(1)) 
                call ga_dadd(1d0,g_vxc(1),1d0,g_tmp(1),g_vxc(1))
c
c               destroy work space
                if (.not. ga_destroy(g_tmp(1))) call errquit
     $             ('xc_getv: ga corrupt?',0, GA_ERR)
                if (.not. ga_destroy(g_tmp(2))) call errquit
     $             ('xc_getv: ga corrupt?',0, GA_ERR)
              end if
            else  ! with CDFIT
c
c             set up prefactors
              kfac(1) = -0.5d0*xfac(1)
              jfac(1) = 0.0d0
c
c             calculate the non-CAM exchange
              if (.not.cam_exch) then ! for regular calculations with CDFIT
                call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
              else
c
c               calculate the CAM exchange
                call case_setflags(.true.)  ! turn on attenuation
                call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
                call case_setflags(.false.) ! turn off attenuation
              end if
              Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
            endif
         else  ! spin-polarized calculations
            if (CDFIT) then
              jfac(1)=0.d0
              jfac(2)=0.d0
              kfac(1)=-1.0d0*xfac(1)
              kfac(2)=-1.0d0*xfac(1)
              if (.not.cam_exch) then ! for regular calculations
               call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_vxc, .false.)
              else
c 
c              calculate the CAM exchange
               call case_setflags(.true.)  ! turn on attenuation
               call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_vxc, .false.)
               call case_setflags(.false.) ! turn off attenuation
              end if
              Exc(1) = Exc(1)+0.5d0*(ga_ddot(g_dens(1),g_vxc(1)) +
     &              ga_ddot(g_dens(2),g_vxc(2)))
            else
              if (.not.cam_exch) then ! for regular calculations
               jfac(1) = 1.0d0
               jfac(2) = 0.0d0
               jfac(3) = 1.0d0
               jfac(4) = 0.0d0
               kfac(1) = 0.0d0
               kfac(2) = 1.0d0
               kfac(3) = 0.0d0
               kfac(4) = 1.0d0
               g_jk(1) = g_vxc(1) ! This assignment is assumed
               g_jk(2) = g_vxc(2)
               g_jk(3) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
               g_jk(4) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
               call ga_zero(g_jk(3))
               call ga_zero(g_jk(4))
               g_d(1)  = g_dens(1)
               g_d(2)  = g_dens(1)
               g_d(3)  = g_dens(2)
               g_d(4)  = g_dens(2)
               call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &              tol2e, oskel, g_d(1), g_jk(1), .false.)
               ecoul = 0.5d0*( ! Alpha coulomb energy
     $              ga_ddot(g_dens(1),g_jk(1)) + 
     $              ga_ddot(g_dens(1),g_jk(3)))
               ecoul = ecoul + 0.5d0*( ! Beta coulomb energy
     $              ga_ddot(g_dens(2),g_jk(1)) + 
     $              ga_ddot(g_dens(2),g_jk(3)))
               exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy
     $              ga_ddot(g_dens(1),g_jk(2)) +
     $              ga_ddot(g_dens(2),g_jk(4)))
               call ga_dadd(1.0d0, g_jk(1), 1.0d0, g_jk(3), g_jk(1))
               call ga_copy(g_jk(1), g_jk(3))
               call ga_dadd(1.0d0, g_jk(1), -xfac(1), g_jk(2),
     $              g_jk(1))
               call ga_dadd(1.0d0, g_jk(3), -xfac(1), g_jk(4),
     $              g_jk(2))
               if (.not. ga_destroy(g_jk(3))) call errquit
     $              ('xc_getv: ga corrupt?',0, GA_ERR)
               if (.not. ga_destroy(g_jk(4))) call errquit
     $              ('xc_getv: ga corrupt?',1, GA_ERR)
              else
c
c              Allocate some scratch space
               g_tmp(1)=ga_create_atom_blocked(geom, ao_bas_han,'tmp1')
               g_tmp(2)=ga_create_atom_blocked(geom, ao_bas_han,'tmp2')
c
c              Calculate Coulomb
               jfac(1) = 1.0d0
               jfac(2) = 1.0d0
               kfac(1) = 0.0d0
               kfac(2) = 0.0d0
               call ga_zero(g_tmp(1))
               call ga_zero(g_tmp(2))
               call case_setflags(.false.)
               call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_tmp, .false.)
c
c              Accumulate contribution 
               call ga_dadd(1.0d0, g_vxc(1), 1.0d0, g_tmp(1), g_vxc(1))
               call ga_dadd(1.0d0, g_vxc(2), 1.0d0, g_tmp(2), g_vxc(2))
               call ga_dadd(1.0d0, g_vxc(1), 1.0d0, g_vxc(2), g_vxc(1))
               call ga_copy(g_vxc(1), g_vxc(2))
               ecoul = 0.5d0*( ! Alpha coulomb energy
     $              ga_ddot(g_dens(1),g_tmp(1)) + 
     $              ga_ddot(g_dens(1),g_tmp(2)))
               ecoul = ecoul + 0.5d0*( ! Beta coulomb energy
     $              ga_ddot(g_dens(2),g_tmp(1)) + 
     $              ga_ddot(g_dens(2),g_tmp(2)))
c
c              Calculate Exchange
               jfac(1) = 0.0d0
               jfac(2) = 0.0d0
               kfac(1) =-1.0d0*xfac(1)
               kfac(2) =-1.0d0*xfac(1)
               call ga_zero(g_tmp(1))
               call ga_zero(g_tmp(2))
               call case_setflags(.true.) ! turn on attenuation
               call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_tmp, .false.)
               call case_setflags(.false.) ! turn off attenuation
c
c              Accumulate contribution 
               call ga_dadd(1.0d0, g_vxc(1), 1.0d0, g_tmp(1), g_vxc(1))
               call ga_dadd(1.0d0, g_vxc(2), 1.0d0, g_tmp(2), g_vxc(2))
               exc(1) = exc(1) + 0.5d0*( ! Exchange energy
     $              ga_ddot(g_dens(1),g_tmp(1)) +
     $              ga_ddot(g_dens(2),g_tmp(2)))
c
c              Deallocate scratch
               if (.not. ga_destroy(g_tmp(1))) call errquit
     $              ('xc_getv: ga corrupt?',0, GA_ERR)
               if (.not. ga_destroy(g_tmp(2))) call errquit
     $              ('xc_getv: ga corrupt?',1, GA_ERR)
c
              end if
            endif
         endif
         if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
         call ga_sync
      endif
      if(util_print('dft timings', print_high))
     &      time2_2e=util_cpusec()   ! end 2e build time
c
c     print fock_2e build time
c
      if(util_print('dft timings', print_high)) then
       if (me.eq.0) then
         write(*,"(4x,'Fock_2e Build Time:',F13.1,'s')") 
     &              time2_2e-time1_2e
       endif
      end if
c
c     Get the DFT exchange-correlation contribution
      if(util_print('dft timings', print_high)) 
     &      time1_xc=util_cpusec()   ! start xc build time
      if (xc_gotxc()) then  
         if(xcfit) then       
            nTrows = nbf_xc
            nTcols = ipol
            if (.not.ma_push_get(MT_Dbl,nTrows*nTcols,'Tmat',lTmat,
     &           iTmat))call errquit('xc_getv: cannot allocate Tmat',0,
     &       MA_ERR)
            call dfill(nTrows*nTcols,0.D0,dbl_mb(iTmat),1)
         endif
         
         if(havehfxc.or.(.not.cdfit)) then
               if(.not.ga_duplicate(g_vxc(1),g_truevxc(1),'g vxc 1'))
     .         call errquit('xcgetv: gaduplicate failed',1, GA_ERR)
               call ga_zero(g_truevxc(1))
               if(ipol.eq.2) then
                  if(.not.ga_duplicate(g_vxc(2),g_truevxc(2),'gv21'))
     .         call errquit('xcgetv: gaduplicate failed',1, GA_ERR)
                  call ga_zero(g_truevxc(2))
               endif
         else
               g_truevxc(1)=g_vxc(1)
               g_truevxc(2)=g_vxc(2)
         endif
c
         if(xc_gotxc())            !! Qin, to skip the call when only HF EX used
     &        call grid_quadv0(rtdb, g_dens, g_truevxc, 
     &           nexc,rho_n,  Exc, dbl_mb(itmat))
c
              if(havehfxc.or.(.not.cdfit)) then
                call ga_dadd(1d0,g_vxc(1),1d0,g_truevxc(1),g_vxc(1))
                if (.not. ga_destroy(g_truevxc(1))) call errquit(
     &              ' xc_getv: ga_destroy failed ',0, GA_ERR)
                if(ipol.eq.2) then
                 call ga_dadd(1d0,g_vxc(2),1d0,g_truevxc(2),g_vxc(2))
                 if (.not. ga_destroy(g_truevxc(2))) call errquit(
     &              ' xc_getv: ga_destroy failed ',0, GA_ERR)
                endif
              endif
         if(util_print('dft timings', print_high)) 
     &         time2_xc=util_cpusec()   ! end xc build time
c
c        print fock_xc build time
         if(util_print('dft timings', print_high)) then
          if (me.eq.0) then
           write(*,"(4x,'Fock_xc Build Time:',F13.1,'s')") 
     &                 time2_xc-time1_xc
          endif
         end if
c
c        In case we are performing an xc fit calculation
         if(xcfit) then
c     
c     symmetrize the "T" vector
c     
            if (oskel)then
               call sym_vec_symmetrize(
     .              geom,xc_bas_han,Dbl_MB(iTmat))
               if (ipol.gt.1)then
                  call sym_vec_symmetrize(geom, xc_bas_han, 
     &                    Dbl_MB(iTmat+nbf_xc))
               endif
            endif
            call xc_fitv(rtdb,Dbl_MB(iTmat), nTrows, nTcols,
     &           g_vxc, g_xcinv, IOLGC)  
            if (.not.ma_pop_stack(lTmat))
     &           call errquit('xc_getv: cannot pop stack',0, MA_ERR)
c     
         endif
      endif
c
      if (oprint_intermediate_xc)then
         write(*,*)' Fock XC matrix leaving xc_getv: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
      endif
c
      return
      end
