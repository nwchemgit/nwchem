      subroutine tce_mo1e(g_ao1e,d_mo1e)
c
c $Id: tce_mo1e.F,v 1.2 2002-10-23 01:38:53 sohirata Exp $
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
#include "tce.fh"
      integer g_ao1e(2)           ! AO Fock matrices
      integer l_ao1e,k_ao1e       ! AO Fock matrices
      integer l_mo1e,k_mo1e       ! MO Fock matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset ! SF size and offset
      integer dummy               ! Not referenced
      integer d_mo1e,request      ! Used for SF
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
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. 0) then
              sf_size = sf_size + int_mb(k_range+g1b-1)
     1                            * int_mb(k_range+g2b-1)
            endif
          endif
        enddo
      enddo
      if (.not.ma_push_get(mt_dbl,sf_size,'MO Fock',
     1  l_mo1e,k_mo1e)) call errquit('tce_mo1e: MA problem',3)
c
c     ==============
c     Transformation
c     ==============
c
      sf_offset=0
      do g1b = 1,noa+nob+nva+nvb
        do g2b = 1,noa+nob+nva+nvb
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. 0) then
              spin = int_mb(k_spin+g1b-1)
              range_g1 = int_mb(k_range+g1b-1)
              range_g2 = int_mb(k_range+g2b-1)
              offset_g1 = int_mb(k_offset+g1b-1)*nbf
              offset_g2 = int_mb(k_offset+g2b-1)*nbf
              if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',
     1          l_work,k_work)) call errquit('tce_mo1e: MA problem',0)
              if (.not.ma_push_get(mt_dbl,nbf*nbf,'AO Fock',
     1          l_ao1e,k_ao1e)) call errquit('tce_mo1e: MA problem',1)
              call ga_get(g_ao1e(spin),1,nbf,1,nbf,dbl_mb(k_ao1e),nbf)
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_ao1e),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_ao1e))
     1          call errquit('tce_mo1e: MA problem',2)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_mo1e+sf_offset),range_g1)
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_mo1e: MA problem',5)
              sf_offset = sf_offset + range_g1 * range_g2
            endif
          endif
        enddo
      enddo
c
c     ===========
c     Write to SF
c     ===========
c
      dummy = sf_write(d_mo1e,0.0d0,dfloat(bytes*sf_size),
     1  dbl_mb(k_mo1e),request)
      dummy = sf_wait(request)
c
c     ===================
c     Close SF and return
c     ===================
c
      if (nodezero.and.util_print('mo1e',print_debug)) then
        call sf_print(d_mo1e,sf_size)
      endif
      if (.not.ma_pop_stack(l_mo1e))
     1  call errquit('tce_mo1e: MA problem',6)
      if (.not.ga_destroy(g_ao1e(1)))
     1  call errquit('tce_mo1e: GA problem',2)
      if (.not.ga_destroy(g_ao1e(2)))
     1  call errquit('tce_mo1e: GA problem',3)
      return
      end
