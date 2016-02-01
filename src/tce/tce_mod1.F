      subroutine tce_mod1(g_aod1,d_mod1)
!
! $Id$
!
!     Spin-spatial-symmetry blocked dipole matrix transformations.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer g_aod1(3)           ! AO dipole matrices
      integer l_aod1,k_aod1       ! AO dipole matrices
      integer l_mod1,k_mod1       ! MO dipole matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset   ! SF size and offset
      integer d_mod1(3)           ! File handle
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      integer axis                ! 1=x, 2=y, 3=z
      logical nodezero
      integer sym_abelian_axis    ! Move this to "sym.fh"
      external sym_abelian_axis
!
      nodezero = (ga_nodeid().eq.0)
      do axis = 1,3
        irrep_d = sym_abelian_axis(geom,axis)
        sf_size=0
        do g1b = 1,noa+nob+nva+nvb
          do g2b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              sf_size = sf_size + int_mb(k_range+g1b-1)
     1                          * int_mb(k_range+g2b-1)
            endif
            endif
            endif
          enddo
        enddo
        if (.not.ma_push_get(mt_dbl,sf_size,'MO dipole',
     1    l_mod1,k_mod1)) call errquit('tce_mod1: MA problem',3,MA_ERR)
        sf_offset=0
        do g2b = 1,noa+nob+nva+nvb
          do g1b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              range_g1 = int_mb(k_range+g1b-1)
              range_g2 = int_mb(k_range+g2b-1)
              offset_g1 = int_mb(k_offset+g1b-1)*nbf
              offset_g2 = int_mb(k_offset+g2b-1)*nbf
              if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',
     1          l_work,k_work)) call errquit('tce_mod1: MA problem',0,
     2          MA_ERR)
              if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO dipole',
     1          l_aod1,k_aod1)) call errquit('tce_mod1: MA problem',1,
     2          MA_ERR)
              call ga_get(g_aod1(axis),1,nbf,1,nbf,dbl_mb(k_aod1),nbf)
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_aod1),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_aod1))
     1          call errquit('tce_mod1: MA problem',2,MA_ERR)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_mod1+sf_offset),range_g1)
              if (nodezero.and.util_print('mod1',print_debug))
     1          call ma_print(dbl_mb(k_mod1+sf_offset),
     2          range_g1,range_g2,'Spin symmetry block of dipole')
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_mod1: MA problem',5,MA_ERR)
              sf_offset = sf_offset + range_g1 * range_g2
            endif
            endif
            endif
          enddo
        enddo
        call put_block(d_mod1(axis),dbl_mb(k_mod1),sf_size,0)
        if (.not.ma_pop_stack(l_mod1))
     1    call errquit('tce_mod1: MA problem',6,MA_ERR)
        if (nodezero.and.util_print('mod1',print_debug))
     1    call sf_print(d_mod1(axis),sf_size)
        if (.not.ga_destroy(g_aod1(axis)))
     1    call errquit('tce_mod1: GA problem',axis,GA_ERR)
      enddo
      return
      end


      subroutine tce_moq1(g_aoq1,d_moq1)
!
! $Id$
!
!     Spin-spatial-symmetry blocked quadrupole matrix transformations.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer g_aoq1(6)           ! AO quadrupole matrices
      integer l_aoq1,k_aoq1       ! AO quadrupole matrices
      integer l_moq1,k_moq1       ! MO quadrupole matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset   ! SF size and offset
      integer d_moq1(6)           ! File handle
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      integer axis                ! 
      logical nodezero
      integer sym_abelian_axis    ! Move this to "sym.fh"
      external sym_abelian_axis
!
      nodezero = (ga_nodeid().eq.0)
      do axis = 1,6
        irrep_d = sym_abelian_axis(geom,axis+3)
        sf_size=0
        do g1b = 1,noa+nob+nva+nvb
          do g2b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              sf_size = sf_size + int_mb(k_range+g1b-1)
     1                          * int_mb(k_range+g2b-1)
            endif
            endif
            endif
          enddo
        enddo
        if (.not.ma_push_get(mt_dbl,sf_size,'MO quadrupole',
     1    l_moq1,k_moq1)) call errquit('tce_moq1: MA problem',3,MA_ERR)
        sf_offset=0
        do g2b = 1,noa+nob+nva+nvb
          do g1b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              range_g1 = int_mb(k_range+g1b-1)
              range_g2 = int_mb(k_range+g2b-1)
              offset_g1 = int_mb(k_offset+g1b-1)*nbf
              offset_g2 = int_mb(k_offset+g2b-1)*nbf
              if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',
     1          l_work,k_work)) call errquit('tce_moq1: MA problem',0,
     2          MA_ERR)
              if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO quadrupole',
     1          l_aoq1,k_aoq1)) call errquit('tce_moq1: MA problem',1,
     2          MA_ERR)
              call ga_get(g_aoq1(axis),1,nbf,1,nbf,dbl_mb(k_aoq1),nbf)
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_aoq1),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_aoq1))
     1          call errquit('tce_moq1: MA problem',2,MA_ERR)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_moq1+sf_offset),range_g1)
              if (nodezero.and.util_print('moq1',print_debug))
     1          call ma_print(dbl_mb(k_moq1+sf_offset),
     2          range_g1,range_g2,'Spin symmetry block of quadrupole')
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_moq1: MA problem',5,MA_ERR)
              sf_offset = sf_offset + range_g1 * range_g2
            endif
            endif
            endif
          enddo
        enddo
        call put_block(d_moq1(axis),dbl_mb(k_moq1),sf_size,0)
        if (.not.ma_pop_stack(l_moq1))
     1    call errquit('tce_moq1: MA problem',6,MA_ERR)
        if (nodezero.and.util_print('moq1',print_debug))
     1    call sf_print(d_moq1(axis),sf_size)
        if (.not.ga_destroy(g_aoq1(axis)))
     1    call errquit('tce_moq1: GA problem',axis,GA_ERR)
      enddo
      return
      end



      subroutine tce_moo1(g_aoo1,d_moo1)
!
! $Id$
!
!     Spin-spatial-symmetry blocked octupole matrix transformations.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer g_aoo1(10)          ! AO octupole matrices
      integer l_aoo1,k_aoo1       ! AO octupole matrices
      integer l_moo1,k_moo1       ! MO octupole matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset   ! SF size and offset
      integer d_moo1(10)          ! File handle
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      integer axis                !
      logical nodezero
      integer sym_abelian_axis    ! Move this to "sym.fh"
      external sym_abelian_axis
!
      nodezero = (ga_nodeid().eq.0)
      do axis = 1,10
        irrep_d = sym_abelian_axis(geom,axis+9)
        sf_size=0
        do g1b = 1,noa+nob+nva+nvb
          do g2b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              sf_size = sf_size + int_mb(k_range+g1b-1)
     1                          * int_mb(k_range+g2b-1)
            endif
            endif
            endif
          enddo
        enddo
        if (.not.ma_push_get(mt_dbl,sf_size,'MO octupole',
     1    l_moo1,k_moo1)) call errquit('tce_moo1: MA problem',3,MA_ERR)
        sf_offset=0
        do g2b = 1,noa+nob+nva+nvb
          do g1b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              range_g1 = int_mb(k_range+g1b-1)
              range_g2 = int_mb(k_range+g2b-1)
              offset_g1 = int_mb(k_offset+g1b-1)*nbf
              offset_g2 = int_mb(k_offset+g2b-1)*nbf
              if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',
     1          l_work,k_work)) call errquit('tce_moo1: MA problem',0,
     2          MA_ERR)
              if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO octupole',
     1          l_aoo1,k_aoo1)) call errquit('tce_moo1: MA problem',1,
     2          MA_ERR)
              call ga_get(g_aoo1(axis),1,nbf,1,nbf,dbl_mb(k_aoo1),nbf)
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_aoo1),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_aoo1))
     1          call errquit('tce_moo1: MA problem',2,MA_ERR)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_moo1+sf_offset),range_g1)
              if (nodezero.and.util_print('moo1',print_debug))
     1          call ma_print(dbl_mb(k_moo1+sf_offset),
     2          range_g1,range_g2,'Spin symmetry block of octupole')
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_moo1: MA problem',5,MA_ERR)
              sf_offset = sf_offset + range_g1 * range_g2
            endif
            endif
            endif
          enddo
        enddo
        call put_block(d_moo1(axis),dbl_mb(k_moo1),sf_size,0)
        if (.not.ma_pop_stack(l_moo1))
     1    call errquit('tce_moo1: MA problem',6,MA_ERR)
        if (nodezero.and.util_print('moo1',print_debug))
     1    call sf_print(d_moo1(axis),sf_size)
        if (.not.ga_destroy(g_aoo1(axis)))
     1    call errquit('tce_moo1: GA problem',axis,GA_ERR)
      enddo
      return
      end



      subroutine tce_mom1(g_aom1,d_mom1,axes)
!
! $Id$
!
!     Spin-spatial-symmetry blocked multipole matrix transformations.
!
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer g_aom1(*)           ! AO dipole matrices
      integer l_aom1,k_aom1       ! AO dipole matrices
      integer l_mom1,k_mom1       ! MO dipole matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset   ! SF size and offset
      integer d_mom1(*)           ! File handle
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      integer axis,axes,ashift
      logical nodezero
      integer sym_abelian_axis    ! Move this to "sym.fh"
      external sym_abelian_axis
!
! Dipole:
! x=1 y=2 z=3
!
! Quadrupole:
! xx=4 yy=5 zz=6 xy=7 xz=8 yz=9
!
! Octupole:
! xxx=10 xxy=11 xxz=12 xyy=13 xyz=14 xzz=15 yyy=16 yyz=17 yzz=18 zzz=19
!
      if (axes.eq. 3) ashift = 0 ! Dipole
      if (axes.eq. 6) ashift = 3 ! Quadrupole
      if (axes.eq.10) ashift = 9 ! Octupole
      nodezero = (ga_nodeid().eq.0)
      do axis = 1,axes
        irrep_d = sym_abelian_axis(geom,axis+ashift)
        sf_size=0
        do g1b = 1,noa+nob+nva+nvb
          do g2b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              sf_size = sf_size + int_mb(k_range+g1b-1)
     1                          * int_mb(k_range+g2b-1)
            endif
            endif
            endif
          enddo
        enddo
        if (.not.ma_push_get(mt_dbl,sf_size,'MO multipole',
     1    l_mom1,k_mom1)) call errquit('tce_mom1: MA problem',3,MA_ERR)
        sf_offset=0
        do g2b = 1,noa+nob+nva+nvb
          do g1b = 1,noa+nob+nva+nvb
            if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1        +int_mb(k_spin+g2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. irrep_d) then
              range_g1 = int_mb(k_range+g1b-1)
              range_g2 = int_mb(k_range+g2b-1)
              offset_g1 = int_mb(k_offset+g1b-1)*nbf
              offset_g2 = int_mb(k_offset+g2b-1)*nbf
              if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',
     1          l_work,k_work)) call errquit('tce_mom1: MA problem',0,
     2          MA_ERR)
              if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO multipole',
     1          l_aom1,k_aom1)) call errquit('tce_mom1: MA problem',1,
     2          MA_ERR)
              call ga_get(g_aom1(axis),1,nbf,1,nbf,dbl_mb(k_aom1),nbf)
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_aom1),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_aom1))
     1          call errquit('tce_mom1: MA problem',2,MA_ERR)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_mom1+sf_offset),range_g1)
              if (nodezero.and.util_print('mom1',print_debug))
     1          call ma_print(dbl_mb(k_mom1+sf_offset),
     2          range_g1,range_g2,'Spin symmetry block of dipole')
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_mom1: MA problem',5,MA_ERR)
              sf_offset = sf_offset + range_g1 * range_g2
            endif
            endif
            endif
          enddo
        enddo
        call put_block(d_mom1(axis),dbl_mb(k_mom1),sf_size,0)
        if (.not.ma_pop_stack(l_mom1))
     1    call errquit('tce_mom1: MA problem',6,MA_ERR)
        if (nodezero.and.util_print('mom1',print_debug))
     1    call sf_print(d_mom1(axis),sf_size)
      enddo
      do axis = 1,axes
        if (.not.ga_destroy(g_aom1(axis)))
     1    call errquit('tce_mom1: GA problem',axis,GA_ERR)
      enddo
      return
      end
