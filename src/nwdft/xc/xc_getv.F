      Subroutine xc_getv(rtdb, Exc, nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, g_wght, g_xyz, 
     &                   wght_GA, rho_n, rdens_atom, cntoce,
     &                   cntobfr, cetobfr, natoms)
c
C$Id: xc_getv.F,v 1.17 1997-12-08 23:48:51 d3h449 Exp $
c
      implicit none
c      
      integer nExc
      integer iVxc_opt
      integer g_xcinv, g_dens(2), g_vxc(4), g_wght, g_xyz
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
      integer cntoce(nshells_ao), cntobfr(2,nshells_ao), 
     &        cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
      integer d(4),f(4)
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
c--> XC Energy
c
      double precision Exc(2)
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
c
c      integer g_scr,me,nproc,i,nTrows,nbf_xc,nbf_ao,nTcols
      integer g_scr,me,nproc,i,nTrows,nTcols
c      integer lTmat,iTmat,itol2e
      integer lTmat,iTmat
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time
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
      nproc=ga_nnodes()
c
      if (oprint_intermediate_xc)then
         write(*,*)' rtdb, Exc, nExc, iVxc_opt, g_xcinv: ',
     &               rtdb, Exc, nExc, iVxc_opt, g_xcinv
         write(*,*)' g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA: ',
     &               g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA
         write(*,*)' Fock XC matrix entering xc_getv: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
      endif
c
      if (abs(xfac(1)).gt.1e-8)then
c       
c        Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
         tol2e=10.d0**(-itol2e)
         if (ipol.eq.1) then
            jfac(1)=0.
            kfac(1) = -0.5d0*xfac(1)
            if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
            call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &                   tol2e, oskel, g_dens(1), g_vxc(1))
            if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
         else
c
c           bogus Ga
c
c            call gacrea(geom, AO_bas_han, nbf_ao, nbf_ao, 
c     &                  'ga scratch', g_scr, 'atom')
            g_scr = ga_create_atom_blocked(geom, AO_bas_han, 
     &                                     'ga scratch')
            call ga_zero(g_scr)
            jfac(1)=0.d0
            jfac(2)=0.d0
            jfac(3)=0.d0
            jfac(4)=0.d0
            kfac(1)=0.0d0
            kfac(2)=-1.0d0*xfac(1)
            kfac(3)=0.0d0
            kfac(4)=-1.0d0*xfac(1)
            d(1) = g_dens(1)
            d(2) = g_dens(1)
            d(3) = g_dens(2)
            d(4) = g_dens(2)
c            f(1) = g_a_coul
            f(1) = g_scr
            f(2) = g_vxc(1)
c            f(3) = g_b_coul
            f(3) = g_scr
            f(4) = g_vxc(2)
            if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
            call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &                   tol2e, oskel, d, f)
            if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
c            call gadest(g_scr)
            if (.not. ga_destroy(g_scr)) call errquit
     &         ('xc_getv: could not destroy g_scr', 0)
         endif
         call ga_sync
c
c        Symmetrize Vxc?
c
c         if (oskel)then
c            call sym_symmetrize(geom, AO_bas_han, .false., g_vxc(1))
c            if (ipol.gt.1)then
c               call sym_symmetrize(geom, AO_bas_han, .false., 
c     &              g_vxc(2))
c            endif
c         endif
c       
c        Compute the exact exchange energy.
c       
         if (ipol.eq.1)then
            Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
         else
            Exc(1) = Exc(1)+0.5d0*(ga_ddot(g_dens(1),g_vxc(1)) +
     &                             ga_ddot(g_dens(2),g_vxc(2)))
         end if
      endif
c
      tot=-xfac(1)
      do i=1,10
        tot=tot+xfac(i)+cfac(i)
      enddo
c
      if (abs(tot).gt.1e-8) then
         if (iVxc_opt.eq.0) then
c         
c           Option A:  Compute via direct numerical quadrature.
c         
            if (oprint_time)call dft_tstamp('Before call to xc_quadv0.')
            call xc_quadv0(rtdb, nExc, g_vxc, g_dens, g_wght, g_xyz, 
     &                     wght_GA, Exc, rho_n, rdens_atom, cntoce,
     &                     cntobfr, cetobfr, natoms)
            if (oprint_time)call dft_tstamp(' After call to xc_quadv0.')
c      write(6,*)'from xc_getv; Exc(1), Exc(2): ',Exc(1), Exc(2)
c
c           Symmetrize Vxc?
c
c            if (oskel)then
c               call sym_symmetrize(geom, AO_bas_han, .false., g_vxc(1))
c               if (ipol.gt.1)then
c                  call sym_symmetrize(geom, AO_bas_han, .false., 
c     &                 g_vxc(2))
c               endif
c            endif
c
         elseif (iVxc_opt.eq.1 )then
c         
c           Option B:  Compute via a LSQ fitting procedure.
c         
c           Define the dimensions of the "T" matrix.
c         
            if (.not. bas_numbf(XC_bas_han,nbf_xc) )then
               call errquit('Exiting in getvxc.',1)
            endif
c
            nTrows = nbf_xc
            nTcols = ipol
c         
c           Allocate scratch space for the "T" matrix.
c         
            if (.not.ma_push_get(MT_Dbl,nTrows*nTcols,'Tmat',lTmat,
     &         iTmat))call errquit('xc_getv: cannot allocate Tmat',0)
            call dfill(nTrows*nTcols,0.D0,dbl_mb(iTmat),1)
c
            if (oprint_time)call dft_tstamp('Before call to xc_quadv1.')
            call xc_quadv1(rtdb, Dbl_MB(iTmat), nExc, g_dens, Exc, 
     &                     rho_n, rdens_atom, cntoce, cntobfr, cetobfr, 
     &                     natoms)
            if (oprint_time)call dft_tstamp(' After call to xc_quadv1.')
c
c           symmetrize the "T" vector
c
            if (oskel)then
               call sym_vec_symmetrize(geom, xc_bas_han, Dbl_MB(iTmat))
               if (ipol.gt.1)then
                  call sym_vec_symmetrize(geom, xc_bas_han, 
     &                 Dbl_MB(iTmat+nbf_xc))
               endif
            endif
c      
            if (oprint_time)call dft_tstamp('Before call to xc_fitv.  ')
            call xc_fitv(Dbl_MB(iTmat), nTrows, nTcols,
     &                   g_vxc, g_xcinv, IOLGC)
            if (oprint_time)call dft_tstamp(' After call to xc_fitv.  ')
c
c           Symmetrize Vxc?
c
c            if (oskel)then
c               call sym_symmetrize(geom, AO_bas_han, .false., g_vxc(1))
c               if (ipol.gt.1)then
c                  call sym_symmetrize(geom, AO_bas_han, .false., 
c     &                 g_vxc(2))
c               endif
c            endif
c
            if (.not.ma_pop_stack(lTmat))
     &         call errquit('xc_getv: cannot pop stack',0)
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
