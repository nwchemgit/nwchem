      subroutine tce_mod1(g_aod1,d_mod1)
c
c $Id: tce_mod1.F,v 1.5 2008-01-14 15:20:35 jhammond Exp $
c
c     Spin-spatial-symmetry blocked dipole matrix transformations.
c
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
c
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
      enddo
      if (.not.ga_destroy(g_aod1(1)))
     1  call errquit('tce_mod1: GA problem',2,GA_ERR)
      if (.not.ga_destroy(g_aod1(2)))
     1  call errquit('tce_mod1: GA problem',3,GA_ERR)
      if (.not.ga_destroy(g_aod1(3)))
     1  call errquit('tce_mod1: GA problem',4,GA_ERR)
      return
      end
