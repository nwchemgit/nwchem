      subroutine tce_mod1(g_aod1,d_mod1,dipole)
c
c $Id: tce_mod1.F,v 1.1 2003-07-25 17:56:57 sohirata Exp $
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
#include "tce.fh"
      double precision dipole(3)  ! Total dipole moments
      integer g_aod1(3)           ! AO dipole matrices
      integer l_aod1,k_aod1       ! AO dipole matrices
      integer l_mod1,k_mod1       ! MO dipole matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset   ! SF size and offset
      integer d_mod1(3)           ! File handle
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      integer g2
      integer axis                ! 1=x, 2=y, 3=z
      logical nodezero
      integer sym_abelian_axis    ! Move this to "sym.fh"
      external sym_abelian_axis
      double precision factor
c
c     =====================================
c     Determine the size of SF and allocate
c     =====================================
c
      nodezero = (ga_nodeid().eq.0)
      sf_size=0
      do g1b = 1,noa+nob+nva+nvb
        do g2b = 1,noa+nob+nva+nvb
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1      +int_mb(k_spin+g2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1      .eq. 0) then
            sf_size = sf_size + int_mb(k_range+g1b-1)
     1                        * int_mb(k_range+g2b-1)
          endif
          endif
          endif
        enddo
      enddo
      if (.not.ma_push_get(mt_dbl,sf_size,'MO dipole',
     1  l_mod1,k_mod1)) call errquit('tce_mod1: MA problem',3)
c
c     ==============
c     Transformation
c     ==============
c
      if (restricted) then
        factor = 2.0d0
      else
        factor = 1.0d0
      endif
      do axis = 1,3
        irrep_d = sym_abelian_axis(geom,axis)
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
     1          l_work,k_work)) call errquit('tce_mod1: MA problem',0)
              if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO dipole',
     1          l_aod1,k_aod1)) call errquit('tce_mod1: MA problem',1)
              call ga_get(g_aod1(axis),1,nbf,1,nbf,dbl_mb(k_aod1),nbf)
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_aod1),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_aod1))
     1          call errquit('tce_mod1: MA problem',2)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_mod1+sf_offset),range_g1)
              if (nodezero.and.util_print('mod1',print_debug))
     1          call ma_print(dbl_mb(k_mod1+sf_offset),
     2          range_g1,range_g2,'Spin symmetry block of dipole')
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_mod1: MA problem',5)
              if ((g2b.eq.g1b).and.(g2b.le.noa+nob)) then
                do g2 = 1,range_g2
                  dipole(axis) = dipole(axis) + factor *
     1            dbl_mb(k_mod1+sf_offset+(g2-1)*range_g2+g2-1)
                enddo
              endif
              sf_offset = sf_offset + range_g1 * range_g2
            endif
            endif
            endif
          enddo
        enddo
        call put_block(d_mod1(axis),dbl_mb(k_mod1),sf_size,0)
        if (nodezero.and.util_print('mod1',print_debug))
     1    call sf_print(d_mod1(axis),sf_size)
      enddo
      if (.not.ma_pop_stack(l_mod1))
     1  call errquit('tce_mod1: MA problem',6)
      if (.not.ga_destroy(g_aod1(1)))
     1  call errquit('tce_mod1: GA problem',2)
      if (.not.ga_destroy(g_aod1(2)))
     1  call errquit('tce_mod1: GA problem',3)
      if (.not.ga_destroy(g_aod1(3)))
     1  call errquit('tce_mod1: GA problem',4)
      return
      end
