      subroutine tce_aod1(rtdb,g_aod1,dipole)
!
! $Id$
!
!     Returns a GA handle for AO dipole matrix.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce_main.fh"
      integer rtdb               ! Run-time database
      integer g_aod1(3)          ! GA handle for AO dipole matrices
      integer g_dens(2)          ! GA handle for density matrices
      double precision dipole(3) ! SCF dipole moments
      integer i
!
!     ===============
!     Allocate g_dens
!     ===============
!
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO density',
     1    -1,-1,g_dens(i))) call errquit('tce_aod1: GA problem',2,
     2    GA_ERR)
      enddo
!
!     =====================
!     Make density matrices
!     =====================
!
      do i=1,ipol
        call ga_dgemm('N','T',nbf,nbf,nocc(i),1.0d0,
     1    g_movecs(i),g_movecs(i),0.0d0,g_dens(i))
        if (ipol.eq.1) call ga_dscal(g_dens(i),2.0d0)
        if (util_print('aod1',
     1    print_debug)) call ga_print(g_dens(i))
      enddo
!
!     ===============
!     Allocate g_aod1
!     ===============
!
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole X',
     1  -1,-1,g_aod1(1))) call errquit('tce_aod1: GA problem',0,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole Y',
     1  -1,-1,g_aod1(2))) call errquit('tce_aod1: GA problem',1,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole Z',
     1  -1,-1,g_aod1(3))) call errquit('tce_aod1: GA problem',2,
     2  GA_ERR)
!
!     ====================
!     Make dipole matrices
!     ====================
!
      call int_dip_ga(ao_bas_han,ao_bas_han,
     1  g_aod1(1),g_aod1(2),g_aod1(3))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(1))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(2))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(3))
!
!     ===================
!     Make dipole moments
!     ===================
!
!     Note: do not attempt to compute these in MO basis;
!     frozen core approximation can give you wrong results.
!
      do i=1,ipol
        dipole(1) = dipole(1) + ga_ddot(g_dens(i),g_aod1(1))
        dipole(2) = dipole(2) + ga_ddot(g_dens(i),g_aod1(2))
        dipole(3) = dipole(3) + ga_ddot(g_dens(i),g_aod1(3))
      enddo
      do i=1,ipol
        if (.not.ga_destroy(g_dens(i)))
     1  call errquit('tce_aod1: failed to destroy g_dens',0,GA_ERR)
      enddo
      return
      end



      subroutine tce_aoq1(rtdb,g_aoq1,quadrupole)
!
! $Id$
!
!     Returns a GA handle for AO quadrupole matrix.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce_main.fh"
      integer rtdb                   ! Run-time database
      integer g_aoq1(6)              ! GA handle for AO quadrupole matrices
      integer g_dens(2)              ! GA handle for density matrices
      double precision quadrupole(6) ! SCF quadrupole moments
      integer i
!
!     ===============
!     Allocate g_dens
!     ===============
!
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO density',
     1    -1,-1,g_dens(i))) call errquit('tce_aoq1: GA problem',2,
     2    GA_ERR)
      enddo
!
!     =====================
!     Make density matrices
!     =====================
!
      do i=1,ipol
        call ga_dgemm('N','T',nbf,nbf,nocc(i),1.0d0,
     1    g_movecs(i),g_movecs(i),0.0d0,g_dens(i))
        if (ipol.eq.1) call ga_dscal(g_dens(i),2.0d0)
        if (util_print('aoq1',
     1    print_debug)) call ga_print(g_dens(i))
      enddo
!
!     ===============
!     Allocate g_aoq1
!     ===============
!
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO quadrupole XX',
     1  -1,-1,g_aoq1(1))) call errquit('tce_aoq1: GA problem',0,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO quadrupole XY',
     1  -1,-1,g_aoq1(2))) call errquit('tce_aoq1: GA problem',1,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO quadrupole XZ',
     1  -1,-1,g_aoq1(3))) call errquit('tce_aoq1: GA problem',2,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO quadrupole YY',
     1  -1,-1,g_aoq1(4))) call errquit('tce_aoq1: GA problem',3,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO quadrupole YZ',
     1  -1,-1,g_aoq1(5))) call errquit('tce_aoq1: GA problem',4,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO quadrupole ZZ',
     1  -1,-1,g_aoq1(6))) call errquit('tce_aoq1: GA problem',5,
     2  GA_ERR)
!
!     ========================
!     Make quadrupole matrices
!     ========================
!
      call int_qdr_ga(ao_bas_han,ao_bas_han,
     1                g_aoq1(1),g_aoq1(2),g_aoq1(3),
     2                g_aoq1(4),g_aoq1(5),g_aoq1(6))
      if (util_print('aoq1',print_debug)) then
        do i=1,6
            call ga_print(g_aoq1(i))
        enddo
      endif
!
!     =======================
!     Make quadrupole moments
!     =======================
!
!     Note: do not attempt to compute these in MO basis;
!     frozen core approximation can give you wrong results.
!
      do i=1,ipol
        quadrupole(1) = quadrupole(1) + ga_ddot(g_dens(i),g_aoq1(1))
        quadrupole(2) = quadrupole(2) + ga_ddot(g_dens(i),g_aoq1(2))
        quadrupole(3) = quadrupole(3) + ga_ddot(g_dens(i),g_aoq1(3))
        quadrupole(4) = quadrupole(4) + ga_ddot(g_dens(i),g_aoq1(4))
        quadrupole(5) = quadrupole(5) + ga_ddot(g_dens(i),g_aoq1(5))
        quadrupole(6) = quadrupole(6) + ga_ddot(g_dens(i),g_aoq1(6))
      enddo
      do i=1,ipol
        if (.not.ga_destroy(g_dens(i)))
     1  call errquit('tce_aoq1: failed to destroy g_dens',0,GA_ERR)
      enddo
      return
      end




      subroutine tce_aoo1(rtdb,g_aoo1,octpole)
!
! $Id$
!
!     Returns a GA handle for AO octupole matrix.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce_main.fh"
      integer rtdb                   ! Run-time database
      integer g_aoo1(10)             ! GA handle for AO octupole matrices
      integer g_dens(2)              ! GA handle for density matrices
      double precision octpole(10)   ! SCF octupole moments
      integer i
!
!     ===============
!     Allocate g_dens
!     ===============
!
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf,nbf,'AO density',
     1    -1,-1,g_dens(i))) call errquit('tce_aoo1: GA problem',2,
     2    GA_ERR)
      enddo
!
!     =====================
!     Make density matrices
!     =====================
!
      do i=1,ipol
        call ga_dgemm('N','T',nbf,nbf,nocc(i),1.0d0,
     1    g_movecs(i),g_movecs(i),0.0d0,g_dens(i))
        if (ipol.eq.1) call ga_dscal(g_dens(i),2.0d0)
        if (util_print('aoo1',
     1    print_debug)) call ga_print(g_dens(i))
      enddo
!
!     ===============
!     Allocate g_aoo1
!     ===============
!
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole XXX',
     1  -1,-1,g_aoo1(1))) call errquit('tce_aoo1: GA problem',0,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole XXY',
     1  -1,-1,g_aoo1(2))) call errquit('tce_aoo1: GA problem',1,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole XXZ',
     1  -1,-1,g_aoo1(3))) call errquit('tce_aoo1: GA problem',2,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole XYY',
     1  -1,-1,g_aoo1(4))) call errquit('tce_aoo1: GA problem',3,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole XYZ',
     1  -1,-1,g_aoo1(5))) call errquit('tce_aoo1: GA problem',4,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole XZZ',
     1  -1,-1,g_aoo1(6))) call errquit('tce_aoo1: GA problem',5,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole YYY',
     1  -1,-1,g_aoo1(7))) call errquit('tce_aoo1: GA problem',6,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole YYZ',
     1  -1,-1,g_aoo1(8))) call errquit('tce_aoo1: GA problem',7,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole YZZ',
     1  -1,-1,g_aoo1(9))) call errquit('tce_aoo1: GA problem',8,
     2  GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO octupole ZZZ',
     1  -1,-1,g_aoo1(10))) call errquit('tce_aoo1: GA problem',9,
     2  GA_ERR)
!
!     ======================
!     Make octupole matrices
!     ======================
!
!    octpole ao integer calculation order int_oct_ga: xxx, xxy, xxz, xyy, xyz,
!                                                     xzz, yyy, yyz, yzz, zzz
!
      call int_oct_ga(ao_bas_han,ao_bas_han,
     1          g_aoo1(1),g_aoo1(2),g_aoo1(3),g_aoo1(4),g_aoo1(5),
     2          g_aoo1(6),g_aoo1(7),g_aoo1(8),g_aoo1(9),g_aoo1(10))
      if (util_print('aoo1',print_debug)) then
        do i=1,10
            call ga_print(g_aoo1(i))
        enddo
      endif
!
!     =======================
!     Make octupole moments
!     =======================
!
!     Note: do not attempt to compute these in MO basis;
!     frozen core approximation can give you wrong results.
!
      do i=1,ipol
        octpole(1)  = octpole(1)  + ga_ddot(g_dens(i),g_aoo1(1))
        octpole(2)  = octpole(2)  + ga_ddot(g_dens(i),g_aoo1(2))
        octpole(3)  = octpole(3)  + ga_ddot(g_dens(i),g_aoo1(3))
        octpole(4)  = octpole(4)  + ga_ddot(g_dens(i),g_aoo1(4))
        octpole(5)  = octpole(5)  + ga_ddot(g_dens(i),g_aoo1(5))
        octpole(6)  = octpole(6)  + ga_ddot(g_dens(i),g_aoo1(6))
        octpole(7)  = octpole(7)  + ga_ddot(g_dens(i),g_aoo1(7))
        octpole(8)  = octpole(8)  + ga_ddot(g_dens(i),g_aoo1(8))
        octpole(9)  = octpole(9)  + ga_ddot(g_dens(i),g_aoo1(9))
        octpole(10) = octpole(10) + ga_ddot(g_dens(i),g_aoo1(10))
      enddo
      do i=1,ipol
        if (.not.ga_destroy(g_dens(i)))
     1  call errquit('tce_aoq1: failed to destroy g_dens',0,GA_ERR)
      enddo
      return
      end
