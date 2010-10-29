      subroutine tce_mod1_new(rtdb,d_mod1,mpole)
c
c $Id$
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
c
c     formerly in tce_aod1_new.F
c
      integer rtdb               ! Run-time database
      double precision origin(3)
      data origin/0d0,0d0,0d0/
c
      integer mpole,mpoles        ! order of multipole, number of blocks
      integer g_aod1              ! AO dipole matrices
      integer l_aod1,k_aod1       ! AO dipole matrices
      integer l_mod1,k_mod1       ! MO dipole matrices
      integer l_work,k_work       ! Work space
      integer sf_size,sf_offset   ! SF size and offset
      integer d_mod1(*)           ! File handle
      integer g1b,g2b             ! Block indexes
      integer range_g1,range_g2   ! Block ranges
      integer offset_g1,offset_g2 ! Block offsets
      integer axis                ! 1=x, 2=y, 3=z
      integer alo(3),ahi(3),ld(2)
      logical nodezero,debugprint
      integer sym_abelian_axis    ! Move this to "sym.fh"
      external sym_abelian_axis
c
      nodezero = (ga_nodeid().eq.0)
      debugprint = .false.
      if (multipole.eq.1) then
        debugprint = debugprint.or.(util_print('aod1',print_debug))
        debugprint = debugprint.or.(util_print('mod1',print_debug))
      elseif (multipole.eq.2) then
        debugprint = debugprint.or.(util_print('aoq1',print_debug))
        debugprint = debugprint.or.(util_print('moq1',print_debug))
      elseif (multipole.eq.3) then
        debugprint = debugprint.or.(util_print('aoo1',print_debug))
        debugprint = debugprint.or.(util_print('moo1',print_debug))
      endif
c
      mpoles = (mpole+1)*(mpole+2)/2
c
c     formerly in tce_aod1_new.F
c
c     ===============
c     Allocate g_aod1
c     ===============
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = mpoles ! formerly 3
      if (.not.nga_create(MT_DBL,3,ahi,'AO dipole ALL',alo,g_aod1))
     &    call errquit('aor_driver: nga_create failed g_aod1',0,GA_ERR)
      call ga_zero(g_aod1)
c
c     ====================
c     Make dipole matrices
c     ====================
c
#ifdef DEBUG_PRINT
      write(LuOut,*) 'calling int_mpole_1ega'
#endif
      if (mpole.eq.1) then
        call int_mpole_1ega(ao_bas_han,ao_bas_han,g_aod1,'dipole',
     &                      origin,.false.)
      elseif (mpole.eq.2) then
        call int_mpole_1ega(ao_bas_han,ao_bas_han,g_aod1,'quadrupole',
     &                      origin,.false.)
      elseif (mpole.eq.3) then
        call int_mpole_1ega(ao_bas_han,ao_bas_han,g_aod1,'octupole',
     &                      origin,.false.)
      endif
      if (debugprint) call ga_print(g_aod1)
c
c     end of what was formerly in tce_aod1_new.F
c
      do axis = 1,mpoles
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
        if (.not.ma_push_get(mt_dbl,sf_size,'MO dipole',l_mod1,k_mod1))
     &      call errquit('tce_mod1_new: MA problem',3,MA_ERR)
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
              if (.not.ma_push_get(mt_dbl,range_g1*nbf,'Work',l_work,
     1          k_work)) call errquit('tce_mod1_new: MA problem',0,
     2          MA_ERR)
              if (.not.ma_push_get(mt_dbl,nbf*nbf*1,'AO dipole',l_aod1,
     1          k_aod1)) call errquit('tce_mod1_new: MA problem',1,
     2          MA_ERR)
c
c             OLD CALL
c              call ga_get(g_aod1(axis),1,nbf,1,nbf,dbl_mb(k_aod1),nbf)
c
c             NEW CALL
              alo(1) = 1
              ahi(1) = nbf
              alo(2) = 1
              ahi(2) = nbf
              alo(3) = axis 
              ahi(3) = axis ! 3 for dipole, 6 for quadrupole, 10 for octupole
              ld(1) = nbf 
              ld(2) = 1
              call nga_get(g_aod1,alo,ahi,dbl_mb(k_aod1),ld)
c
              call dgemm('T','N',range_g1,nbf,nbf,1.0d0,
     1          dbl_mb(k_movecs_sorted+offset_g1),nbf,dbl_mb(k_aod1),
     2          nbf,0.0d0,dbl_mb(k_work),range_g1)
              if (.not.ma_pop_stack(l_aod1))
     1          call errquit('tce_mod1_new: MA problem',2,MA_ERR)
              call dgemm('N','N',range_g1,range_g2,nbf,1.0d0,
     1          dbl_mb(k_work),range_g1,
     2          dbl_mb(k_movecs_sorted+offset_g2),nbf,
     3          0.0d0,dbl_mb(k_mod1+sf_offset),range_g1)
c              if (nodezero.and.util_print('mod1',print_debug))
c     1          call ma_print(dbl_mb(k_mod1+sf_offset),
c     2          range_g1,range_g2,'Spin symmetry block of dipole')
              if (.not.ma_pop_stack(l_work))
     1          call errquit('tce_mod1_new: MA problem',5,MA_ERR)
              sf_offset = sf_offset + range_g1 * range_g2
            endif
            endif
            endif
          enddo
        enddo
        call put_block(d_mod1(axis),dbl_mb(k_mod1),sf_size,0)
        if (.not.ma_pop_stack(l_mod1))
     1    call errquit('tce_mod1_new: MA problem',6,MA_ERR)
c
c       Debug printing
c
        if (debugprint) then
          if (nodezero) write(LuOut,*) 'MO multipole integrals'
          call ga_print(d_mod1)
        endif
c
      enddo
      if (.not.ga_destroy(g_aod1))
     1  call errquit('tce_mod1_new: GA problem',2,GA_ERR)
c
#ifdef DEBUG_PRINT
      write(LuOut,*) 'end of tce_mod1_new'
#endif
c
      return
      end
