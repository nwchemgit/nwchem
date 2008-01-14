      subroutine tce_mo1e(g_ao1e,d_mo1e)
c
c $Id: tce_mo1e.F,v 1.11 2008-01-14 15:20:34 jhammond Exp $
c
c     Spin-spatial-symmetry blocked Fock matrix transformations.
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
      integer g_ao1e(2)           ! AO Fock matrices
      integer l_ao1e,k_ao1e       ! AO Fock matrices
      integer l_mo1e,k_mo1e       ! MO Fock matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset ! SF size and offset
      integer d_mo1e              ! File handle
      integer spin                ! Spin
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      logical nodezero
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
      if (.not.ma_push_get(mt_dbl,sf_size,'MO Fock',
     1  l_mo1e,k_mo1e)) call errquit('tce_mo1e: MA problem',3,MA_ERR)
c
c     ==============
c     Transformation
c     ==============
c
      sf_offset=0
      do g2b = 1,noa+nob+nva+nvb
        do g1b = 1,noa+nob+nva+nvb
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1      +int_mb(k_spin+g2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1      .eq. 0) then
            spin = int_mb(k_spin+g1b-1)
            range_g1 = int_mb(k_range+g1b-1)
            range_g2 = int_mb(k_range+g2b-1)
            offset_g1 = int_mb(k_offset+g1b-1)*nbf
            offset_g2 = int_mb(k_offset+g2b-1)*nbf
            if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',
     1        l_work,k_work)) call errquit('tce_mo1e: MA problem',0,
     2        MA_ERR)
            if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO Fock',
     1        l_ao1e,k_ao1e)) call errquit('tce_mo1e: MA problem',1,
     2        MA_ERR)
            call ga_get(g_ao1e(spin),1,nbf,1,nbf,dbl_mb(k_ao1e),nbf)
            call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1        dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_ao1e),
     2        nbf,0.0d0,dbl_mb(k_work),range_g1)
            if (.not.ma_pop_stack(l_ao1e))
     1        call errquit('tce_mo1e: MA problem',2,MA_ERR)
            call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1        dbl_mb(k_work),range_g1,
     2        dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3        0.0d0,dbl_mb(k_mo1e+sf_offset),range_g1)
            if (nodezero.and.util_print('mo1e',print_debug))
     1        call ma_print(dbl_mb(k_mo1e+sf_offset),
     2        range_g1,range_g2,'Spin symmetry block of Fock')
            if (.not.ma_pop_stack(l_work))
     1        call errquit('tce_mo1e: MA problem',5,MA_ERR)
            sf_offset = sf_offset + range_g1 * range_g2
          endif
          endif
          endif
        enddo
      enddo
c
c     ===========
c     Write to SF
c     ===========
c
      call put_block(d_mo1e,dbl_mb(k_mo1e),sf_size,0)
c
c     ===================
c     Close SF and return
c     ===================
c
      if (nodezero.and.util_print('mo1e',print_debug)) then
        call sf_print(d_mo1e,sf_size)
      endif
      if (.not.ma_pop_stack(l_mo1e))
     1  call errquit('tce_mo1e: MA problem',6,MA_ERR)
      if (.not.ga_destroy(g_ao1e(1)))
     1  call errquit('tce_mo1e: GA problem',2,GA_ERR)
      if (.not.ga_destroy(g_ao1e(2)))
     1  call errquit('tce_mo1e: GA problem',3,GA_ERR)
      return
      end



      subroutine tce_mo1e_epsilon(d_mo1e)
c
c     Spin-spatial-symmetry blocked Fock matrix formation from epsilons.
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
      integer l_mo1e,k_mo1e       ! MO Fock matrices
      integer sf_size,sf_offset   ! SF size and offset
      integer d_mo1e              ! File handle
      integer spin                ! Spin
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer g1,g2
      logical nodezero
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
      if (.not.ma_push_get(mt_dbl,sf_size,'MO Fock',
     1  l_mo1e,k_mo1e)) call errquit('tce_mo1e: MA problem',3,MA_ERR)
c
c     ==============
c     Transformation
c     ==============
c
      sf_offset=0
      do g2b = 1,noa+nob+nva+nvb
        do g1b = 1,noa+nob+nva+nvb
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1      +int_mb(k_spin+g2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1      .eq. 0) then
            spin = int_mb(k_spin+g1b-1)
            range_g1 = int_mb(k_range+g1b-1)
            range_g2 = int_mb(k_range+g2b-1)
            do g2 = 1, range_g2
              do g1 = 1, range_g1
                if ((g1b.eq.g2b).and.(g1.eq.g2)) then
                  dbl_mb(k_mo1e+sf_offset+(g2-1)*range_g1+g1-1)=
     1            dbl_mb(k_evl_sorted+int_mb(k_offset+g1b-1)+g1-1)
                else
                  dbl_mb(k_mo1e+sf_offset+(g2-1)*range_g1+g1-1)=0.0d0
                endif
              enddo
            enddo
            if (nodezero.and.util_print('mo1e',print_debug))
     1        call ma_print(dbl_mb(k_mo1e+sf_offset),
     2        range_g1,range_g2,'Spin symmetry block of Fock')
            sf_offset = sf_offset + range_g1 * range_g2
          endif
          endif
          endif
        enddo
      enddo
c
c     ===========
c     Write to SF
c     ===========
c
      call put_block(d_mo1e,dbl_mb(k_mo1e),sf_size,0)
c
c     ===================
c     Close SF and return
c     ===================
c
      if (nodezero.and.util_print('mo1e',print_debug)) then
        call sf_print(d_mo1e,sf_size)
      endif
      if (.not.ma_pop_stack(l_mo1e))
     1  call errquit('tce_mo1e: MA problem',6,MA_ERR)
      return
      end
