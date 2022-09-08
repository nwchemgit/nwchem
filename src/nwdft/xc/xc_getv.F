c $Id$
c
C> \ingroup nwdft
C> @{
C>
C> \file xc_getv.F
C> Calculate exchange-correlation energy
C> 
C> \brief Calculate the exchange-correlation energy and Fock matrix
C> contributions
C>
C> This driver routine solves for the XC energy and potential (Vxc) via 
C> numerical quadrature methods. The results are obtained either by
C> direct numerical integration or by means of a LSQ fit of the Vxc to
C> a set of Gaussian functions. This fitted function can be used to
C> evaluate Vxc via a summation of a series of 3-center overlap
C> integrals (3OIs). The algorithms are formulated in terms of matrix
C> products. See subsequent subroutines for further explanation.
C>
      Subroutine xc_getv(rtdb, Exc, ecoul,nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, g_wght, g_xyz,g_nq, 
     &                   wght_GA, rho_n, rdens_atom,
     &                   cetobfr, natoms, cd_coef)
c
      implicit none
#include "errquit.fh"
c      
      integer nExc      !< [Input] The number of energy terms
                        !< - nExc=1: Exc(1) = exchange + correlation 
                        !< - nExc=2: Exc(1) = exchange,
                        !<           Exc(2) = correlation 
      integer iVxc_opt  !< [Input] If 1 then do density fitting for 
                        !< exchange
      integer g_xcinv   !< [Work] GA for the inversion of the fitting
                        !< matrix
      integer g_dens(2) !< [Input] The density matrices, if ipol=1
                        !< g_dens(1)=\f$D^\alpha+D^\beta\f$, else
                        !< g_dens(1)=\f$D^\alpha\f$ and 
                        !< g_dens(2)=\f$D^\beta\f$.
      integer g_vxc(4)  !< [Output] DFT Fock matrix contributions, if
                        !< ipol=1 g_vxc(1)=\f$F^\alpha+F^\beta\f$, else
                        !< g_vxc(1)=\f$F^\alpha\f$ and
                        !< g_vxc(2)=\f$F^\beta\f$.
      integer g_wght    !< [Work] The grid point weights if wght_GA
      integer g_xyz     !< [Work] The grid point coordinates if wght_GA
      integer g_nq      !< [Unused]
      integer natoms    !< [Input] The number of atoms
      logical IOLGC     !< [Input] .TRUE. do not use disk for exchange
                        !< fitting, .FALSE. store data on disk
      logical wght_GA   !< [Input] .TRUE. store grid points in GA,
                        !< .FALSE. store grid points on file
      integer rtdb      !< [Input] The RTDB handle
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
#include "dftps.fh"
c
      integer cetobfr(2,natoms) !< [Unused]
      double precision rho_n    !< [Output] The number of electrons
                                !< obtained by integrating the density
      double precision rdens_atom(ipol*natoms*natoms) !< [Unused]
      double precision jfac(4),kfac(4)
      integer g_jk(4), g_d(4)
      logical havehfxc
c     
      integer  ga_create_atom_blocked
cc AJL/Begin/FDE
c      logical xc_gotxc
c      external ga_create_atom_blocked,xc_gotxc
      external ga_create_atom_blocked
cc AJL/End
c
c--> XC Energy
c
      double precision Exc(2) !< [Output] The energy terms
                              !< - nExc=1: Exc(1) = exchange +
                              !<   correlation 
                              !< - nExc=2: Exc(1) = exchange,
                              !<           Exc(2) = correlation 
      double precision ecoul  !< [Output] The Coulomb energy
      double precision cd_coef(nbf_cd,*)
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
cc AJL/Begin/FDE
c      integer me,nTrows,nTcols
      integer me
c      integer lTmat,iTmat,g_truevxc(2)
c AJL: These parameters are unused?
c      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time 
cc AJL/Unused
c     ,     oprint_oep
c      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
cc AJL/End
      double precision tol2e
      integer g_tmp(2)
c
c     timings
c
      double precision time1_2e,time2_2e
cc AJL/Begin/FDE
c      double precision time1_xc,time2_xc
cc AJL/End
c
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
cc AJL/Unused
c      oprint_oep = util_print('oep', print_high) 
      Exc(1)=0.d0
      Exc(2)=0.d0
cc AJL/Begin/FDE
c      iTmat=0
cc AJL/End
c
      me=ga_nodeid()
      havehfxc=abs(xfac(1)).gt.1d-8
c
      if (oprint_intermediate_xc)then
c         write(luout,*)' rtdb, Exc, nExc, iVxc_opt, g_xcinv: ',
c     &               rtdb, Exc, nExc, iVxc_opt, g_xcinv
c         write(luout,*)'g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA:',
c     &               g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA
         write(luout,*)' Fock XC matrix entering xc_getv: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
c         call ga_print(g_dens(1))
c         if(ipol.eq.2)call ga_print(g_dens(2))
      endif
c
      if(oprint_time)
     &      time1_2e=util_cpusec()   ! start 2e build time
      if (havehfxc .or. (.not. CDFIT))then
c       
c        Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
         tol2e=10.d0**(-itol2e)
         call ga_sync
         if (odftps) call pstat_on(ps_f2e)
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
              if (cam_exch) call case_setflags(.true.)
                call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
              if (cam_exch)  call case_setflags(.false.) ! turn off attenuation
              Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
            endif
         else  ! spin-polarized calculations
            if (CDFIT) then
              jfac(1)=0.d0
              jfac(2)=0.d0
              kfac(1)=-1.0d0*xfac(1)
              kfac(2)=-1.0d0*xfac(1)
              if (cam_exch) call case_setflags(.true.)
              call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_vxc, .false.)
              if (cam_exch) call case_setflags(.false.) ! turn off attenuation
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
         if (odftps) call pstat_off(ps_f2e)
         if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
         call ga_sync
      endif
      if(oprint_time)
     &      time2_2e=util_cpusec()   ! end 2e build time
c
c     print fock_2e build time
c
      if(oprint_time) then
       if (me.eq.0) then
         write(luout,"(4x,'Fock_2e Build Time:',F13.1,'s')") 
     &              time2_2e-time1_2e
       endif
      end if
c
c     Get the DFT exchange-correlation contribution
c
cc AJl/Begin/FDE
c  I have moved this to a new subroutine so the XC evaluation can be
c  called multiple times, as is needed for the FDE evaluation of the
c  non-additive XC energy
c
c      if(util_print('dft timings', print_high)) 
c     &      time1_xc=util_cpusec()   ! start xc build time
c      if (xc_gotxc()) then  
c         if(xcfit) then       
c            nTrows = nbf_xc
c            nTcols = ipol
c            if (.not.ma_push_get(MT_Dbl,nTrows*nTcols,'Tmat',lTmat,
c     &           iTmat))call errquit('xc_getv: cannot allocate Tmat',0,
c     &       MA_ERR)
c            call dfill(nTrows*nTcols,0.D0,dbl_mb(iTmat),1)
c         endif
c         
c         if(havehfxc.or.(.not.cdfit)) then
c               if(.not.ga_duplicate(g_vxc(1),g_truevxc(1),'g vxc 1'))
c     .         call errquit('xcgetv: gaduplicate failed',1, GA_ERR)
c               call ga_zero(g_truevxc(1))
c               if(ipol.eq.2) then
c                  if(.not.ga_duplicate(g_vxc(2),g_truevxc(2),'gv21'))
c     .         call errquit('xcgetv: gaduplicate failed',1, GA_ERR)
c                  call ga_zero(g_truevxc(2))
c               endif
c         else
c               g_truevxc(1)=g_vxc(1)
c               g_truevxc(2)=g_vxc(2)
c         endif
cc
c         call grid_quadv0(rtdb, g_dens, g_truevxc, 
c     &                    nexc,rho_n,  Exc, dbl_mb(itmat))
cc
c         if(havehfxc.or.(.not.cdfit)) then
c             call ga_dadd(1d0,g_vxc(1),1d0,g_truevxc(1),g_vxc(1))
c             if (.not. ga_destroy(g_truevxc(1))) call errquit(
c     &           ' xc_getv: ga_destroy failed ',0, GA_ERR)
c             if(ipol.eq.2) then
c                 call ga_dadd(1d0,g_vxc(2),1d0,g_truevxc(2),g_vxc(2))
c                 if (.not. ga_destroy(g_truevxc(2))) call errquit(
c     &               ' xc_getv: ga_destroy failed ',0, GA_ERR)
c             endif
c         endif
c         if(util_print('dft timings', print_high)) 
c     &         time2_xc=util_cpusec()   ! end xc build time
cc
cc        print fock_xc build time
c         if(util_print('dft timings', print_high)) then
c          if (me.eq.0) then
c           write(*,"(4x,'Fock_xc Build Time:',F13.1,'s')") 
c     &                 time2_xc-time1_xc
c          endif
c         end if
cc
cc        In case we are performing an xc fit calculation
c         if(xcfit) then
cc     
cc     symmetrize the "T" vector
cc     
c            if (oskel)then
c               call sym_vec_symmetrize(
c     .              geom,xc_bas_han,Dbl_MB(iTmat))
c               if (ipol.gt.1)then
c                  call sym_vec_symmetrize(geom, xc_bas_han, 
c     &                    Dbl_MB(iTmat+nbf_xc))
c               endif
c            endif
c            call xc_fitv(rtdb,Dbl_MB(iTmat), nTrows, nTcols,
c     &           g_vxc, g_xcinv, IOLGC)  
c            if (.not.ma_pop_stack(lTmat))
c     &           call errquit('xc_getv: cannot pop stack',0, MA_ERR)
cc     
c         endif
c      endif
cc
      if (odftps) call pstat_on(ps_getvxc)
      call xc_getvxc(rtdb, Exc, nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, rho_n, 0, 0, cd_coef)
      if (odftps) call pstat_off(ps_getvxc)
c last two inputs -------------------------------->.fde_option., g_dens_fde
c
cc AJL/End
c
      if (oprint_intermediate_xc)then
         write(luout,*)' Fock XC matrix leaving xc_getv: '
         call util_flush(6)
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
         call util_flush(6)
      endif
c
      return
      end
C>
C> @}
