      Subroutine xc_exso(rtdb,Exc,ecoul,nExc,g_densso,g_fockso)
c
C$Id: xc_exso.F,v 1.2 2000-04-25 15:37:40 edo Exp $
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
c--> XC Energy
c
      double precision Exc(2)
      double precision ecoul ! [output]
      integer me,type
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
      oprint_intermediate_xc = .true. 
      oprint_time = util_print('dft timings', print_high)
c      Exc(1)=0.d0
c      Exc(2)=0.d0
      me=ga_nodeid()
      if (oprint_intermediate_xc)then
         write(*,*)' rtdb, Exc, nExc: ',rtdb, Exc, nExc
         write(*,*)' Fock XC matrix entering xc_getv: '
         call ga_print(g_fockso(1))
         call ga_print(g_fockso(2))
      endif
c
      if (abs(xfac(1)).gt.1e-8)then
         call ga_inquire(g_densso(1),type,nbf,nbf)
         nbf=nbf/2
c     
c     Compute the exact exchange potential (as in Hartree-Fock calculations).
c     
         tol2e=10.d0**(-itol2e)
         call ga_sync
         if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
c
c     ???????? what to do with charge density fit? 
c
         if(.not.ga_create(mt_dbl,nbf,nbf,'temp den',0,0,g_dens(1)))
     &        call errquit('xc_exso: error creating ga', 0)
         if(.not.ga_create(mt_dbl,nbf,nbf,'temp den',0,0,g_dens(2)))
     &        call errquit('xc_exso: error creating ga', 0)
         if(.not.ga_create(mt_dbl,nbf,nbf,'temp den',0,0,g_jk(1)))
     &        call errquit('xc_exso: error creating ga', 0)
         if(.not.ga_create(mt_dbl,nbf,nbf,'temp den',0,0,g_jk(2)))
     &        call errquit('xc_exso: error creating ga', 0)
         if(.not.ga_create(mt_dbl,nbf,nbf,'temp den',0,0,g_jk(3)))
     &        call errquit('xc_exso: error creating ga', 0)
         if(.not.ga_create(mt_dbl,nbf,nbf,'temp den',0,0,g_jk(4)))
     &        call errquit('xc_exso: error creating ga', 0)
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
c
c     Im(Kaa) and Im(Kbb) (g_dens(1) and g_dens(2))
c
         call ga_dens_aabbi(g_dens, g_densso, nbf)
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1))
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy (i.i=-1)
     $        onem*ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(2),1,nbf,1,nbf, 
     &                   -xfac(1),g_jk(2),    1,nbf,1,nbf,
     &                            g_fockso(2),1,nbf,1,nbf)
         call ga_dadd_patch(1.0d0,g_fockso(2),nbf+1,2*nbf,nbf+1,2*nbf, 
     &                   -xfac(1),g_jk(4),    1,nbf,1,nbf,
     &                            g_fockso(2),nbf+1,2*nbf,nbf+1,2*nbf)
c
         if (oprint_intermediate_xc)then
            write(*,*)' Fock XC matrix: Im(Kaa) and Im(Kbb): '
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
         call ga_dens_abri(g_dens, g_densso, nbf)
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1))
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy 
     $        ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(1),1,nbf,nbf+1,2*nbf, 
     &                   -xfac(1),g_jk(2),    1,nbf,1,nbf,
     &                            g_fockso(1),1,nbf,nbf+1,2*nbf)
         call ga_dadd_patch(1.0d0,g_fockso(2),1,nbf,nbf+1,2*nbf, 
     &                   -xfac(1),g_jk(4),    1,nbf,1,nbf,
     &                            g_fockso(2),1,nbf,nbf+1,2*nbf)
         if (oprint_intermediate_xc)then
            write(*,*)' Fock XC matrix: Re(Kab) and Im(Kab): '
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
         call ga_dens_bari(g_dens, g_densso, nbf)
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1))
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy 
     $        ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(1),nbf+1,2*nbf,1,nbf, 
     &                   -xfac(1),g_jk(2),    1,nbf,1,nbf,
     &                            g_fockso(1),nbf+1,2*nbf,1,nbf)
         call ga_dadd_patch(1.0d0,g_fockso(2),nbf+1,2*nbf,1,nbf, 
     &                   -xfac(1),g_jk(4),    1,nbf,1,nbf,
     &                            g_fockso(2),nbf+1,2*nbf,1,nbf)
         if (oprint_intermediate_xc)then
            write(*,*)' Fock XC matrix: Re(Kba) and Im(Kba): '
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
         write(*,*)' Fock XC matrix: Leaving xc_exso: '
         call ga_print(g_fockso(1))
         call ga_print(g_fockso(2))
      endif
      call ga_sync
c     
      return
      end

