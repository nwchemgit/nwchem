      subroutine tce_ao1e_shellfockbuild(rtdb,g_ao1e)
c
c $Id$
c
c     Returns a GA handle for AO Fock matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb               ! Run-time database
      integer g_ao1e(2)          ! GA handle for AO Fock matrices
      integer g_dens(2)          ! GA handle for AO density matrices
      integer g_hcore            ! GA handle for AO hcore (T+V) matrices
c     logical symmetry           ! True if symmetry is used in integral evaluation
c
c 7/1/03 Remember that there is a global variable with the same name (tce.fh)
c
      integer idim(3)            ! Dimension array for NGA create
      integer ichunk(3)          ! Chunk array for NGA create
      integer alo(2),ahi(2)      ! Regions for patch
      integer blo(3),bhi(3)      ! Regions for patch
      integer g_ndens,g_nfock    ! NGA handles for work density and Fock matrices
      double precision jf(4)     ! Prefactor for HF Coulomb as input for shell_fock_build
      double precision kf(4)     ! Prefactor for HF exchange as input for shell_fock_build
      integer i
      integer ndu,nfu
c
c     ====================
c     Initialize int_API's
c     ====================
c
c     symmetry=(sym_number_ops(geom).gt.0).and.sym_abelian_group(geom)
      call scf_get_fock_param(rtdb,tol2e)
c
c     ================================
c     Allocate g_ao1e, g_dens, g_hcore
c     ================================
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'alpha AO Fock',
     1  -1,-1,g_ao1e(1))) call errquit('tce_ao1e: GA problem',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'beta AO Fock',
     1  -1,-1,g_ao1e(2))) call errquit('tce_ao1e: GA problem',1,GA_ERR)
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO density',
     1    -1,-1,g_dens(i))) call errquit('tce_ao1e: GA problem',2,
     2    GA_ERR)
      enddo
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO H core',
     1  -1,-1,g_hcore)) call errquit('tce_ao1e: GA problem',3,GA_ERR)
c
c     =====================
c     Make density matrices
c     =====================
c
      do i=1,ipol
        call ga_dgemm('N','T',nbf,nbf,nocc(i),1.0d0,
     1    g_movecs(i),g_movecs(i),0.0d0,g_dens(i))
        if (ipol.eq.1) call ga_dscal(g_dens(i),2.0d0)
        if (util_print('ao1e',
     1    print_debug)) call ga_print(g_dens(i))
      enddo
c
c     ===========================
c     Make core Hamiltonian (T+V)
c     ===========================
c
      call ga_zero(g_hcore)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'kinetic',.false.)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'potential',.false.)
      if (util_print('ao1e',
     1  print_debug)) call ga_print(g_hcore)
c
c     =====================
c     Call shell_fock_build
c     =====================
c
c     Note: even when a KS reference wave function is used,
c     the Fock operator is defined in a conventional HF way,
c     and is never the KS Hamiltonian.
c
      idim(1)=ipol
      idim(2)=nbf
      idim(3)=nbf
      ichunk(1)=ipol
      ichunk(2)=-1
      ichunk(3)=-1
      if (.not.nga_create(mt_dbl,3,idim,'Work density matrix',
     1  ichunk,g_ndens)) call errquit
     2  ('tce_ao1e: failed to create g_ndens',0,GA_ERR)
      if (.not.nga_create(mt_dbl,3,idim,'Work Fock matrix',
     1  ichunk,g_nfock)) call errquit
     2  ('tce_ao1e: failed to create g_nfock',0,GA_ERR)
      do i=1,ipol
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=nbf
        blo(1)=i
        bhi(1)=i
        blo(2)=1
        bhi(2)=nbf
        blo(3)=1
        bhi(3)=nbf
        call nga_copy_patch('N',g_dens(i),alo,ahi,g_ndens,blo,bhi)
      enddo
      call ga_zero(g_nfock)
      if (ipol.eq.1) then
        jf(1)=1.0d0
        kf(1)=-0.5d0
        call shell_fock_build(geom,ao_bas_han,0,1,
     1    jf,kf,tol2e,g_ndens,g_nfock,.true.)
      else if (ipol.eq.2) then
        jf(1)=1.0d0
        kf(1)=1.0d0
        jf(2)=1.0d0
        kf(2)=1.0d0
        ndu=2
        nfu=2
        call uhf_fock_setup(g_ndens,g_nfock,jf,kf,ndu,
     1    nfu,nbf,.false.)
        call shell_fock_build(geom,ao_bas_han,0,ndu,
     1    jf,kf,tol2e,g_ndens,g_nfock,.true.)
        call uhf_fock_finish(g_ndens,g_nfock,jf,kf,ndu,
     1    nfu,nbf,.false.)
      endif
      do i=1,ipol
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=nbf
        blo(1)=i
        bhi(1)=i
        blo(2)=1
        bhi(2)=nbf
        blo(3)=1
        bhi(3)=nbf
        call nga_copy_patch('N',g_nfock,blo,bhi,g_ao1e(i),alo,ahi)
      enddo
      do i=1,ipol
        call ga_add(1.0d0,g_hcore,1.0d0,g_ao1e(i),g_ao1e(i))
      enddo
      if (ipol.eq.1) call ga_copy(g_ao1e(1),g_ao1e(2))
      if (util_print('ao1e',print_debug)) then
        call ga_print(g_ao1e(1))
        call ga_print(g_ao1e(2))
      endif
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ga_destroy(g_hcore))
     1  call errquit('tce_ao1e: failed to destroy g_hcore',0,GA_ERR)
      do i=1,ipol
        if (.not.ga_destroy(g_dens(i)))
     1  call errquit('tce_ao1e: failed to destroy g_aodens',0,GA_ERR)
      enddo
      if (.not.ga_destroy(g_ndens))
     1  call errquit('tce_ao1e: failed to destroy g_ndens',0,GA_ERR)
      if (.not.ga_destroy(g_nfock))
     1  call errquit('tce_ao1e: failed to destroy g_nfock',0,GA_ERR)
c
c     ======
c     Return
c     ======
c
      do i = 1, ipol
        if (.not.ga_destroy(g_movecs(i)))
     1    call errquit('tce_ao1e: GA problem',0,GA_ERR)
      enddo
      return
      end
      subroutine tce_ao1e_fock2e(rtdb,g_ao1e)
c
c     Returns a GA handle for AO Fock matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb               ! Run-time database
      integer g_ao1e(2)          ! GA handle for AO Fock matrices
      integer g_dens(4)          ! GA handle for AO density matrices
      integer g_fock(4)          ! GA handle for AO density matrices
      integer g_hcore            ! GA handle for AO hcore (T+V) matrices
      double precision jf(4)     ! Prefactor for HF Coulomb
      double precision kf(4)     ! Prefactor for HF exchange
      integer i
c
c     ====================
c     Initialize int_API's
c     ====================
c
      call scf_get_fock_param(rtdb,tol2e)
c
c     ================================
c     Allocate g_ao1e, g_dens, g_hcore
c     ================================
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'alpha AO Fock',
     1  -1,-1,g_ao1e(1))) call errquit('tce_ao1e: GA problem',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'beta AO Fock',
     1  -1,-1,g_ao1e(2))) call errquit('tce_ao1e: GA problem',1,GA_ERR)
      do i=1,4
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO density',
     1    -1,-1,g_dens(i))) call errquit('tce_ao1e: GA problem',2,
     2    GA_ERR)
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO Fock',
     1    -1,-1,g_fock(i))) call errquit('tce_ao1e: GA problem',3,
     2    GA_ERR)
        call ga_zero(g_fock(i))
      enddo
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO H core',
     1  -1,-1,g_hcore)) call errquit('tce_ao1e: GA problem',4,GA_ERR)
c
c     =====================
c     Make density matrices
c     =====================
c
      do i=1,ipol
        call ga_dgemm('N','T',nbf,nbf,nocc(i),1.0d0,
     1    g_movecs(i),g_movecs(i),0.0d0,g_dens(i))
        if (ipol.eq.1) call ga_dscal(g_dens(i),2.0d0)
        if (util_print('ao1e',
     1    print_debug)) call ga_print(g_dens(i))
      enddo
      if (ipol.eq.2) then
        call ga_copy(g_dens(1),g_dens(3))
        call ga_copy(g_dens(2),g_dens(4))
      endif
c
c     ===========================
c     Make core Hamiltonian (T+V)
c     ===========================
c
      call ga_zero(g_hcore)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'kinetic',.false.)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'potential',.false.)
      if (util_print('ao1e',
     1  print_debug)) call ga_print(g_hcore)
c
c     ============
c     Call fock_2e
c     ============
c
c     Note: even when a KS reference wave function is used,
c     the Fock operator is defined in a conventional HF way,
c     and is never the KS Hamiltonian.
c
      if (ipol.eq.1) then
        jf(1)=1.0d0
        kf(1)=-0.5d0
        call fock_2e(geom,ao_bas_han,1,
     1    jf,kf,tol2e,.false.,g_dens,g_fock,.false.)
      else if (ipol.eq.2) then
        jf(1)=1.0d0
        kf(1)=0.0d0
        jf(2)=1.0d0
        kf(2)=0.0d0
        jf(3)=0.0d0
        kf(3)=1.0d0
        jf(4)=0.0d0
        kf(4)=1.0d0
        call fock_2e(geom,ao_bas_han,4,
     1    jf,kf,tol2e,.false.,g_dens,g_fock,.false.)
      endif
      if (ipol.eq.1) then
        call ga_add(1.0d0,g_hcore,1.0d0,g_fock(1),g_ao1e(1))
        call ga_copy(g_ao1e(1),g_ao1e(2))
      else
        call ga_add(1.0d0,g_hcore,1.0d0,g_fock(1),g_ao1e(1))
        call ga_add(1.0d0,g_ao1e(1),1.0d0,g_fock(2),g_ao1e(1))
        call ga_add(1.0d0,g_ao1e(1),-1.0d0,g_fock(3),g_ao1e(1))
        call ga_add(1.0d0,g_hcore,1.0d0,g_fock(1),g_ao1e(2))
        call ga_add(1.0d0,g_ao1e(2),1.0d0,g_fock(2),g_ao1e(2))
        call ga_add(1.0d0,g_ao1e(2),-1.0d0,g_fock(4),g_ao1e(2))
      endif
      if (util_print('ao1e',print_debug)) then
        call ga_print(g_ao1e(1))
        call ga_print(g_ao1e(2))
      endif
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ga_destroy(g_hcore))
     1  call errquit('tce_ao1e: failed to destroy g_hcore',0,GA_ERR)
      do i=1,4
        if (.not.ga_destroy(g_dens(i)))
     1  call errquit('tce_ao1e: failed to destroy g_dens',0,GA_ERR)
        if (.not.ga_destroy(g_fock(i)))
     1  call errquit('tce_ao1e: failed to destroy g_fock',0,GA_ERR)
      enddo
c
c     ======
c     Return
c     ======
c
      return
      end
