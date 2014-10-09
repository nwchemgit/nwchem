      subroutine int_2c_ga(ibas, jbas, g, integ_type, Alpha)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
C
      Integer NxTask            ! From TCGMSG
      External NxTask
c
c     Compute the desired type of integrals (kinetic, potential,  & overlap
C     1-electron integrals and 2-center 2-electron ERIs)
c     scale them by ALPHA and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c
c     arguments
c
      integer ibas, jbas        ! [input] bra and ket basis sets
      integer g                 ! [output] GA handle to array 
      character*(*) integ_type  ! [input] Name of integrals to compute
      Double Precision Alpha    ! [input] Scaling factor for these integrals
C
C     Paramters
C
      Integer Task_Chunk
      Parameter (Task_Chunk = 10)  ! Large enough to minimize overhead
c
c     local variables
c
      integer nshell_i, nshell_j
      integer ishell, jshell, nproc, mem1, max1e, Next
      integer ijshell, ilo, ihi, jlo, jhi, jdim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
c
      if (inp_compare(.false., integ_type, 'kinetic')) then
         type = 1
      else if (inp_compare(.false., integ_type, 'potential')) then
         type = 2
      else if (inp_compare(.false., integ_type, 'overlap')) then
         type = 3
      else if (inp_compare(.false., integ_type, 'eri')) then
         type = 4
      else
         write(6,*) ' integ_type = ', integ_type
         call errquit('int_2c_ga: unknown integ_type', 0, GA_ERR)
      endif
c
c     Get info about the basis sets
c
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('int_2c_ga: bas_numcont failed for ibas', ibas, BASIS_ERR)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('int_2c_ga: bas_numcont failed for jbas', jbas, BASIS_ERR)
c
c     allocate necessary local temporary arrays on the stack
c
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c
c     k_* are the offsets corrsponding to the l_* handles
c
      if ( type .ge. 1 .AND. type .le. 3) then
         call int_mem_1e(max1e, mem1)
      ElseIf ( type .eq. 4) then
         call int_mem_2e2c(max1e, mem1)
      EndIf
      if (.not. MA_push_get(MT_DBL,max1e,'int_2c_ga:buf',l_buf,k_buf))
     $     call errquit('int_2c_ga: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'int_2c_ga:scr',l_scr,k_scr))
     $     call errquit('int_2c_ga: ma failed', mem1, MA_ERR)
c
c     Loop thru shells with static parallel work decomposition
c
      call ga_sync
      nproc = ga_nnodes()
      ijshell = 0
      Next = NxTask( nproc, Task_Chunk)
C
      do ishell = 1, nshell_i
         if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) call errquit
     $        ('int_2c_ga: bas_cn2bfr failed ibas', ibas, BASIS_ERR)
         do jshell = 1, nshell_j
            if ( ijshell .eq. Next) then
               if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $              call errquit('int_2c_ga: bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
               jdim = jhi - jlo + 1
c
c     Generate the integrals
c
               if (type .eq. 1) then
                  call int_1eke (ibas, ishell, jbas, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else if (type .eq. 2) then
                  call int_1epe (ibas, ishell, jbas, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else if (type .eq. 3) then
                  call int_1eov (ibas, ishell, jbas, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else if (type .eq. 4) then
                  call int_2e2c (ibas, ishell, jbas, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else
                  call errquit('int_2c_ga: invalid type?', type,
     &       INPUT_ERR)
               endif
c
c     Add the integrals into the global array
c
               call ga_acc(g, jlo, jhi, ilo, ihi, dbl_mb(k_buf), jdim,
     $              Alpha)
               Next = NxTask( NProc, Task_Chunk)
            endif
            ijshell = ijshell + 1
         enddo
      enddo
      Next = NxTask( -NProc, Task_Chunk)
c
c     chop stack at first item allocated
c
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_2c_ga: pop failed', 0, MA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_2c_ga: pop failed', 0, MA_ERR)
c
      call ga_sync()            ! So that no nasty races can result
c
      end
