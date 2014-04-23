      Subroutine xc_exso(rtdb,Exc,ecoul,nExc,g_densso,g_fockso)
c
C$Id$
c
      implicit none
#include "errquit.fh"
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
#include "case.fh"
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
      double precision tol2e, denmax1, denmax2, denmax
c
      call ga_inquire(g_densso(1),type,nao,nao)
      nao=nao/2
c
c     Check if this is a cam type calculation
      call case_setflags(.false.)
      if (cam_exch) call case_setflags(.true.)
c         
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
      if(oprint_intermediate_xc)then
c         call ga_print(g_fockso(1))
c         call ga_print(g_fockso(2))
c         call ga_print(g_densso(1))
c         call ga_print(g_densso(2))
      endif
c
      if (abs(xfac(1)).gt.1e-8)then
c     
c     Compute the exact exchange potential (as in Hartree-Fock calculations).
c     
         tol2e=10.d0**(-itol2e)
         call ga_sync
         if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
c
c     ???????? what to do with charge density fit? 
c
         if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(1)))
     &        call errquit('xc_exso: error creating ga', 0, GA_ERR)
         if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(2)))
     &        call errquit('xc_exso: error creating ga', 0, GA_ERR)
         g_jk(1) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         g_jk(2) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         g_jk(3) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         g_jk(4) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
         jfac(1)=zero
         jfac(2)=zero
         jfac(3)=zero
         jfac(4)=zero
         kfac(1)=zero
         kfac(2)=one
         kfac(3)=zero
         kfac(4)=one 
c
c     Im(Kaa) and Im(Kbb) (g_dens(1) and g_dens(2))
c
         call ga_dens_aabbi(g_dens, g_densso, nao)
         denmax1 = 0.0d0
         denmax2 = 0.0d0
         denmax=0.0d0
         call ga_maxelt(g_dens(1), denmax1)
         call ga_maxelt(g_dens(2), denmax2)
         denmax=max(denmax1,denmax2)
         if(denmax.lt.1.0d-6)goto 1001
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         if(oprint_intermediate_xc)then 
c            call ga_print(g_dens(1))
c            call ga_print(g_dens(2))
         endif
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1), .true.)
         call ga_scale(g_jk(2),dble(-1.0d0))
         call ga_scale(g_jk(4),dble(-1.0d0))
         call asym_fock2e(g_jk(2))
         call asym_fock2e(g_jk(4))
c
         exc(1) = exc(1) + xfac(1)*0.5d0*( ! All exchange energy (i.i=-1)
     $        onem*ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(2),1,nao,1,nao, 
     &                    xfac(1),g_jk(2),    1,nao,1,nao,
     &                            g_fockso(2),1,nao,1,nao)
         call ga_dadd_patch(1.0d0,g_fockso(2),nao+1,2*nao,nao+1,2*nao, 
     &                    xfac(1),g_jk(4),    1,nao,1,nao,
     &                            g_fockso(2),nao+1,2*nao,nao+1,2*nao)
c
         if (oprint_intermediate_xc)then
            if(me.eq.0)write(*,*)'Fock XC matrix:Im(Kaa) and Im(Kbb):'
            call ga_print(g_jk(2))
            call ga_print(g_jk(4))
         endif
 1001    continue
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
c
c     Re(Kab) and Re(Kba) (g_dens(1) and g_dens(2)) ???
c     HvD: g_dens(1:2) are anti-symmetric and therefore Im(D) so we
c     should be calculating Im(Kab) and Im(Kba) unless I am very
c     mistaken.
c     
         call ga_dens_abbar(g_dens, g_densso, nao)
         denmax1 = 0.0d0
         denmax2 = 0.0d0
         denmax=0.0d0
         call ga_maxelt(g_dens(1), denmax1)
         call ga_maxelt(g_dens(2), denmax2)
         denmax=max(denmax1,denmax2)
         if(denmax.lt.1.0d-6)goto 1002
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1), .true.)
         call ga_scale(g_jk(2),dble(-1.0d0))
         call ga_scale(g_jk(4),dble(-1.0d0))
         call ga_transpose(g_jk(4),g_jk(1))
         call ga_add(dble(0.5d0),g_jk(2),dble(0.5d0),g_jk(1), 
     &        g_jk(2))          !g_jk(2) has Kab
         call ga_transpose(g_jk(2),g_jk(4))
c
         exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy 
     $        ga_ddot(g_dens(1),g_jk(2)) +
     $        ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(1),1,nao,nao+1,2*nao, 
     &                   -xfac(1),g_jk(2),    1,nao,1,nao,
     &                            g_fockso(1),1,nao,nao+1,2*nao)
         call ga_dadd_patch(1.0d0,g_fockso(1),nao+1,2*nao,1,nao, 
     &                   -xfac(1),g_jk(4),    1,nao,1,nao,
     &                            g_fockso(1),nao+1,2*nao,1,nao)
         if (oprint_intermediate_xc)then
            if(me.eq.0)write(*,*)'Fock XC matrix:Re(Kab) and Im(Kab):'
            call ga_print(g_jk(2))
            call ga_print(g_jk(4))
         endif
c
 1002    continue 
         call ga_zero(g_jk(1))
         call ga_zero(g_jk(2))
         call ga_zero(g_jk(3))
         call ga_zero(g_jk(4))
c
c     Im(Kab) and Im(Kba) (g_dens(1) and g_dens(2))
c
         call ga_dens_abbai(g_dens, g_densso, nao)
         denmax1 = 0.0d0
         denmax2 = 0.0d0
         denmax=0.0d0
         call ga_maxelt(g_dens(1), denmax1)
         call ga_maxelt(g_dens(2), denmax2)
         denmax=max(denmax1,denmax2)
         if(denmax.lt.1.0d-6)goto 1003
         g_d(1)  = g_dens(1)
         g_d(2)  = g_dens(1)
         g_d(3)  = g_dens(2)
         g_d(4)  = g_dens(2)
         call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &        tol2e, oskel, g_d(1), g_jk(1), .true.)
         call ga_scale(g_jk(2),dble(-1.0d0))
         call ga_scale(g_jk(4),dble(-1.0d0))
         call ga_transpose(g_jk(4),g_jk(1))
         call ga_add(dble(0.5d0),g_jk(2),dble(-0.5d0),g_jk(1),
     &        g_jk(2))
         call ga_transpose(g_jk(2),g_jk(4))
         call ga_scale(g_jk(4),dble(-1.0d0))
c
         exc(1) = exc(1) + xfac(1)*0.5d0*( ! All exchange energy 
     $        onem*ga_ddot(g_dens(1),g_jk(2)) +
     $        onem*ga_ddot(g_dens(2),g_jk(4)))
         call ga_dadd_patch(1.0d0,g_fockso(2),1,nao,nao+1,2*nao,
     &                    xfac(1),g_jk(2),    1,nao,1,nao,
     &                            g_fockso(2),1,nao,nao+1,2*nao)
         call ga_dadd_patch(1.0d0,g_fockso(2),nao+1,2*nao,1,nao, 
     &                    xfac(1),g_jk(4),    1,nao,1,nao,
     &                            g_fockso(2),nao+1,2*nao,1,nao)
         if (oprint_intermediate_xc)then
            if(me.eq.0)write(*,*)'Fock XC matrix:Re(Kba) and Im(Kba):'
            call ga_print(g_jk(2))
            call ga_print(g_jk(4))
         endif
c
 1003    continue
         if (.not. ga_destroy(g_dens(1))) call errquit
     $        ('xc_exso: ga corrupt?',0, GA_ERR)
         if (.not. ga_destroy(g_dens(2))) call errquit
     $        ('xc_exso: ga corrupt?',1, GA_ERR)
         if (.not. ga_destroy(g_jk(1))) call errquit
     $        ('xc_exso: ga corrupt?',0, GA_ERR)
         if (.not. ga_destroy(g_jk(2))) call errquit
     $        ('xc_exso: ga corrupt?',1, GA_ERR)
         if (.not. ga_destroy(g_jk(3))) call errquit
     $        ('xc_exso: ga corrupt?',0, GA_ERR)
         if (.not. ga_destroy(g_jk(4))) call errquit
     $        ('xc_exso: ga corrupt?',1, GA_ERR)
      endif
      if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
      if (oprint_intermediate_xc)then
         if(me.eq.0)write(*,*)'Fock XC matrix:Leaving xc_exso:'
c         call ga_print(g_fockso(1))
c         call ga_print(g_fockso(2))
      endif
      call ga_sync
c     
c     Unset attenuation before exiting
      if (cam_exch) call case_setflags(.false.)

      return
      end
c
      subroutine asym_fock2e(g_fock)
      implicit none
#include "errquit.fh"
c
#include "global.fh"
c
      integer g_fock, g_tmp
c
      if (.not. ga_duplicate(g_fock, g_tmp, 'tmp'))
     $     call errquit('asym_fock2e: ga duplicate?',0, GA_ERR)
c
      call ga_transpose(g_fock,g_tmp)
      call ga_add(0.5d0,g_fock,-0.5d0,g_tmp,g_fock)

      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('asym_fock2e: ga corrupt?',0, GA_ERR)
c
      return
      end
