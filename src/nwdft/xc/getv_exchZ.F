c
      subroutine getv_exchZ(rtdb, exchZ, g_densso, g_fockso, g_v, which)
c
      implicit none
#include "errquit.fh"
c      
      integer rtdb, g_densso(2), g_fockso(2), g_v(4)
      double precision exchZ
      character*4 which
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
c
      double precision jfac(4),kfac(4)
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer me,type
      double precision tol2e
      integer g_tmp(2),g_dens(2)
      integer nao
c
      me=ga_nodeid()
      tol2e=10.d0**(-itol2e)
      call ga_sync
c
c     Query the full density matrix
      call ga_inquire(g_densso(1),type,nao,nao)
      nao=nao/2  ! get the number of atomic orbitals
c
c     Allocate some scratch space
      g_tmp(1)=ga_create_atom_blocked(geom, ao_bas_han,'tmp1')
      g_tmp(2)=ga_create_atom_blocked(geom, ao_bas_han,'tmp2')
      if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(1)))
     &        call errquit('xc_exso: error creating ga', 0, GA_ERR)
      if(.not.ga_create(mt_dbl,nao,nao,'temp den',0,0,g_dens(2)))
     &        call errquit('xc_exso: error creating ga', 0, GA_ERR)
c
      jfac(1) = 0.0d0
      jfac(2) = 0.0d0
      kfac(1) =-0.5d0*xfac(1)
      kfac(2) =-0.5d0*xfac(1)
      call ga_zero(g_tmp(1))
      call ga_zero(g_tmp(2))
c
      exchZ = 0.d0
c
      if (which .eq. 'real') then
c
c        Get real part of the density matrix
         call ga_dens_aabbr(g_dens, g_densso, nao)
c
c        Form the difference density matrix
         call ga_dadd(1.d0, g_dens(1), -1.d0, g_dens(2), g_dens(1))
         call ga_copy(g_dens(1),g_dens(2))
c
         call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_tmp, .false.)
c
c        Accumulate contribution 
         call ga_zero(g_v(1))
         call ga_zero(g_v(2))
         call ga_dadd(1.0d0, g_v(1), 1.0d0, g_tmp(1), g_v(1))
         call ga_dadd(1.0d0, g_v(2), 1.0d0, g_tmp(2), g_v(2))
c        exchZ = 0.5d0*(ga_ddot(g_dens(1),g_v(1)) +
c     $              ga_ddot(g_dens(2),g_v(2)))
         exchZ = 0.5d0*ga_ddot(g_dens(1),g_v(1))
c
c        Patch in real part of exchange
         call ga_sync()
         call ga_dadd_patch(1.d0, g_fockso(1), 1, nao, 1, nao,
     &     1.0d0, g_v(1), 1, nao, 1, nao,
     &     g_fockso(1), 1, nao, 1, nao)
         call ga_dadd_patch(1.d0, g_fockso(1), 1+nao, 2*nao,
     &     1+nao, 2*nao,  -1.0d0, g_v(2),  1, nao,   ! note the -ve sign (Pauli-Z matrix)
     &     1, nao, g_fockso(1), 1+nao, 2*nao, 1+nao, 2*nao)
c
      else  ! imaginary part
c
c        Get imag part of the density matrix
         call ga_dens_aabbi(g_dens, g_densso, nao)
c
c        Form the difference density matrix
         call ga_dadd(1.d0, g_dens(1), -1.d0, g_dens(2), g_dens(1))
         call ga_copy(g_dens(1),g_dens(2))
c
         call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_dens, g_tmp, .true.)  ! true implies no symmetrization
         call asym_fock2e(g_tmp(1))  ! asymmetrize outside
         call asym_fock2e(g_tmp(2))  ! asymmetrize outside
c
c        Accumulate contribution
         call ga_zero(g_v(1))
         call ga_zero(g_v(2))
         call ga_dadd(1.0d0, g_v(1), 1.0d0, g_tmp(1), g_v(1))
         call ga_dadd(1.0d0, g_v(2), 1.0d0, g_tmp(2), g_v(2))
         exchZ = -0.5d0*ga_ddot(g_dens(1),g_v(1))  ! -ve sign from i*i = -1
c
c        == patch in the imag part ==
         call ga_sync()
         call ga_dadd_patch(1.0d0, g_fockso(2), 1, nao, 1, nao,
     &    1.d0, g_v(1), 1, nao, 1, nao,
     &    g_fockso(2), 1, nao, 1, nao)
         call ga_dadd_patch(1.0d0,g_fockso(2),nao+1,2*nao,nao+1,2*nao,
     &    -1.d0, g_v(2), 1, nao, 1, nao,   ! note the -ve signe (Pauli-Z matrix)
     &    g_fockso(2), nao+1, 2*nao, nao+1, 2*nao)
c
      end if  ! which
c
c     Deallocate scratch
      if (.not. ga_destroy(g_tmp(1))) call errquit
     $          ('xc_getv: ga corrupt?',0, GA_ERR)
      if (.not. ga_destroy(g_tmp(2))) call errquit
     $          ('xc_getv: ga corrupt?',1, GA_ERR)
      if (.not. ga_destroy(g_dens(1))) call errquit
     $          ('xc_getv: ga corrupt?',0, GA_ERR)
      if (.not. ga_destroy(g_dens(2))) call errquit
     $          ('xc_getv: ga corrupt?',1, GA_ERR)
c
      return
      end
