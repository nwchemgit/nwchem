      Subroutine xc_exso(rtdb,Exc,ecoul,nExc,g_densso,g_fockso)
c
C$Id: xc_exso.F,v 1.4 2000-04-26 23:29:11 edo Exp $
c
      implicit none
c      
      integer nExc
      integer g_densso(2),g_fockso(2)
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
      double precision jfac(4),kfac(4)
      integer g_jk(4), g_d(4), g_dens(2)
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c     
c--> XC Energy
c
      double precision Exc(2)
      double precision ecoul ! [output]
      integer me,nao,type
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e
c******************************************************************************
c
c     Compute the matrix elements for the exact exchange contribution: 
c     Kaa, Kab and Kbb. 
c     The contributions to Kaa and Kbb from Re(Daa) and Re(Dbb) are already 
c     included through xc_getv by way of g_dens(1) and g_dens(2). The imaginary 
c     parts of Kaa and Kbb will be calculated here. Also are the real and imaginary 
c     parts of Kab. Kba is obtained from Kab by the Hermitian property. 
c     
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
c      oprint_intermediate_xc = .true. 
      oprint_time = util_print('dft timings', print_high)
c      Exc(1)=0.d0
c      Exc(2)=0.d0
      me=ga_nodeid()
      if (oprint_intermediate_xc .and. me .eq. 0)then
         write(*,*)' rtdb, Exc, nExc: ',rtdb, Exc, nExc
         write(*,*)' Fock XC matrix entering xc_getv: '
      endif
!      call ga_print(g_fockso(1))
!      call ga_print(g_fockso(2))
!      call ga_print(g_densso(1))
!      call ga_print(g_densso(2))
c
      if (abs(xfac(1)).gt.1e-8)then
         call ga_inquire(g_densso(1),type,nao,nao)
         nao=nao/2
         if(me.eq.0)write(*,*)"nao=", nao
c     
c     Compute the exact exchange potential (as in Hartree-Fock calculations).
c     
         tol2e=10.d0**(-itol2e)
         if(me.eq.0)write(*,*)"xfac(1),tol2e",xfac(1),tol2e  
         call ga_sync
         if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
c
c     ???????? what to do with charge density fit? 
c
         if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(1)))
     &        call errquit('xc_exso: error creating ga', 0)
         if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(2)))
     &        call errquit('xc_exso: error creating ga', 0)
         g_jk(1) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         g_jk(2) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         g_jk(3) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         g_jk(4) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
c
c     Im(Kaa) and Im(Kbb) (g_dens(1) and g_dens(2))
c
         if(me.eq.0)write(*,*)"before aabbi" 
         call ga_dens_aabbi(g_dens, g_densso, nao)
         if(me.eq.0)write(*,*)"after aabbi" 
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call ga_print(g_dens(1))
         call ga_print(g_dens(2))
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1))
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy (i.i=-1)
     $        onem*ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(2),1,nao,1,nao, 
     &                   -xfac(1),g_jk(2),    1,nao,1,nao,
     &                            g_fockso(2),1,nao,1,nao)
         call ga_dadd_patch(1.0d0,g_fockso(2),nao+1,2*nao,nao+1,2*nao, 
     &                   -xfac(1),g_jk(4),    1,nao,1,nao,
     &                            g_fockso(2),nao+1,2*nao,nao+1,2*nao)
c
         if (oprint_intermediate_xc)then
            if(me.eq.0)write(*,*)'Fock XC matrix:Im(Kaa) and Im(Kbb):'
            call ga_print(g_jk(2))
            call ga_print(g_jk(4))
         endif
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
c
c     Re(Kab) and Im(Kab) (g_dens(1) and g_dens(2))
c
         call ga_dens_abri(g_dens, g_densso, nao)
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1))
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy 
     $        ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(1),1,nao,nao+1,2*nao, 
     &                   -xfac(1),g_jk(2),    1,nao,1,nao,
     &                            g_fockso(1),1,nao,nao+1,2*nao)
         call ga_dadd_patch(1.0d0,g_fockso(2),1,nao,nao+1,2*nao, 
     &                   -xfac(1),g_jk(4),    1,nao,1,nao,
     &                            g_fockso(2),1,nao,nao+1,2*nao)
         if (oprint_intermediate_xc)then
            if(me.eq.0)write(*,*)'Fock XC matrix:Re(Kab) and Im(Kab):'
            call ga_print(g_jk(2))
            call ga_print(g_jk(4))
         endif
c
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
c
c     Re(Kba) and Im(Kba) (g_dens(1) and g_dens(2))
c
         call ga_dens_bari(g_dens, g_densso, nao)
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1))
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy 
     $        ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(1),nao+1,2*nao,1,nao, 
     &                   -xfac(1),g_jk(2),    1,nao,1,nao,
     &                            g_fockso(1),nao+1,2*nao,1,nao)
         call ga_dadd_patch(1.0d0,g_fockso(2),nao+1,2*nao,1,nao, 
     &                   -xfac(1),g_jk(4),    1,nao,1,nao,
     &                            g_fockso(2),nao+1,2*nao,1,nao)
         if (oprint_intermediate_xc)then
            if(me.eq.0)write(*,*)'Fock XC matrix:Re(Kba) and Im(Kba):'
            call ga_print(g_jk(2))
            call ga_print(g_jk(4))
         endif
c
         if (.not. ga_destroy(g_dens(1))) call errquit
     $        ('xc_exso: ga corrupt?',0)
         if (.not. ga_destroy(g_dens(2))) call errquit
     $        ('xc_exso: ga corrupt?',1)
         if (.not. ga_destroy(g_jk(1))) call errquit
     $        ('xc_exso: ga corrupt?',0)
         if (.not. ga_destroy(g_jk(2))) call errquit
     $        ('xc_exso: ga corrupt?',1)
         if (.not. ga_destroy(g_jk(3))) call errquit
     $        ('xc_exso: ga corrupt?',0)
         if (.not. ga_destroy(g_jk(4))) call errquit
     $        ('xc_exso: ga corrupt?',1)
      endif
      if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
      if (oprint_intermediate_xc)then
         if(me.eq.0)write(*,*)'Fock XC matrix:Leaving xc_exso:'
         call ga_print(g_fockso(1))
         call ga_print(g_fockso(2))
      endif
      call ga_sync
c     
      return
      end

