      Subroutine xc_getv(rtdb, Exc, ecoul,nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, g_wght, g_xyz,g_nq, 
     &                   wght_GA, rho_n, rdens_atom,
     &                   cetobfr, natoms)
c
C$Id: xc_getv.F,v 1.38 2000-06-09 17:08:05 edo Exp $
c
      implicit none
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
#include "util.fh"
#include "sym.fh"
c
      integer cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
      integer g_jk(4), g_d(4)
      logical havehfxc
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
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
      integer me,i,nTrows,nTcols
      integer lTmat,iTmat,g_truevxc(2)
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time, grid_on_file
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e,tot
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
      Exc(1)=0.d0
      Exc(2)=0.d0
      me=ga_nodeid()
      havehfxc=abs(xfac(1)).gt.1d-8
c      nproc=ga_nnodes()
c
      if (oprint_intermediate_xc)then
         write(*,*)' rtdb, Exc, nExc, iVxc_opt, g_xcinv: ',
     &               rtdb, Exc, nExc, iVxc_opt, g_xcinv
         write(*,*)' g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA: ',
     &               g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA
         write(*,*)' Fock XC matrix entering xc_getv: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
         call ga_print(g_dens(1))
         if(ipol.eq.2)call ga_print(g_dens(2))
      endif
c
      if (havehfxc .or. (.not. CDFIT))then
c       
c        Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
         tol2e=10.d0**(-itol2e)
         call ga_sync
         if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
         if (ipol.eq.1) then
            kfac(1) = -0.5d0*xfac(1)
            jfac(1)=0.0d0
            if (.not. CDFIT) then
             jfac(2) = 1.0d0
             kfac(2) = 0d0
              g_vxc(2) = ga_create_atom_blocked(geom,ao_bas_han,'jk')
              g_dens(2)=g_dens(1)
              call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
              Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
              ecoul = 0.5d0*ga_ddot(g_dens(1),g_vxc(2)) 
              call ga_dadd(1d0,g_vxc(1),1d0,g_vxc(2),g_vxc(1))
              if (.not. ga_destroy(g_vxc(2))) call errquit
     $             ('xc_getv: ga corrupt?',0)
            else
              call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
            Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
            endif
         else
            if (CDFIT) then
               jfac(1)=0.d0
               jfac(2)=0.d0
               kfac(1)=-1.0d0*xfac(1)
               kfac(2)=-1.0d0*xfac(1)
               call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_vxc, .false.)
               Exc(1) = Exc(1)+0.5d0*(ga_ddot(g_dens(1),g_vxc(1)) +
     &              ga_ddot(g_dens(2),g_vxc(2)))
            else
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
     $              ('xc_getv: ga corrupt?',0)
               if (.not. ga_destroy(g_jk(4))) call errquit
     $              ('xc_getv: ga corrupt?',1)
            endif
         endif
         if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
         call ga_sync
      endif
c
      tot=-xfac(1)
      do i=1,numfunc
        tot=tot+xfac(i)+cfac(i)
      enddo
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &     grid_on_file))then
         grid_on_file = .false.
      endif
c
      if (abs(tot).gt.1e-8) then
         if(xcfit) then       
            nTrows = nbf_xc
            nTcols = ipol
            if (.not.ma_push_get(MT_Dbl,nTrows*nTcols,'Tmat',lTmat,
     &           iTmat))call errquit('xc_getv: cannot allocate Tmat',0)
            call dfill(nTrows*nTcols,0.D0,dbl_mb(iTmat),1)
         endif
         
         if (grid_on_file)then
            if(havehfxc.or.(.not.cdfit)) then
               call ga_duplicate(g_vxc(1),g_truevxc(1),'g vxc 1')
               call ga_zero(g_truevxc(1))
               if(ipol.eq.2) then
                  call ga_duplicate(g_vxc(2),g_truevxc(2),'g vxc 21')
                  call ga_zero(g_truevxc(2))
               endif
            else
               g_truevxc(1)=g_vxc(1)
               g_truevxc(2)=g_vxc(2)
            endif
            call grid_quadv0(rtdb, g_dens, g_truevxc, 
     ,           nexc,rho_n,  Exc, dbl_mb(itmat))
            if(havehfxc.or.(.not.cdfit)) then
               call ga_dadd(1d0,g_vxc(1),1d0,g_truevxc(1),g_vxc(1))
               if (.not. ga_destroy(g_truevxc(1))) call errquit(
     .              ' xc_getv: ga_destroy failed ',0)
               if(ipol.eq.2) then
               call ga_dadd(1d0,g_vxc(2),1d0,g_truevxc(2),g_vxc(2))
               if (.not. ga_destroy(g_truevxc(2))) call errquit(
     .              ' xc_getv: ga_destroy failed ',0)
               endif
            endif
         else
            if(.not.xcfit) then
c     
c     Option A:  Compute via direct numerical quadrature.
c     
c     compute the integrated density (numerical)
c     
c     
               call xc_quadv0(rtdb, nExc, g_vxc,g_dens,g_wght,
     /              g_xyz,g_nq, 
     &              wght_GA, Exc, rho_n, rdens_atom, 
     &              cetobfr, natoms)
               
            else
c     
c     Option B:  Compute via a LSQ fitting procedure.
c     
c     Define the dimensions of the "T" matrix.
c     
               if (.not. bas_numbf(XC_bas_han,nbf_xc) )then
                  call errquit('Exiting in getvxc.',1)
               endif
c     
c     
               call xc_quadv1(rtdb, Dbl_MB(iTmat), nExc, g_dens, Exc, 
     &              rho_n, rdens_atom, cetobfr, 
     &              natoms)
            endif
         endif
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
            call xc_fitv(Dbl_MB(iTmat), nTrows, nTcols,
     &           g_vxc, g_xcinv, IOLGC)
            if (.not.ma_pop_stack(lTmat))
     &           call errquit('xc_getv: cannot pop stack',0)
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
c     
      return
      end

