      subroutine int_1e_ga(ibas, jbas, g, integ_type)
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
c
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c
c     arguments
c
      integer ibas, jbas        ! [input] bra and ket basis sets
      integer g                 ! [output] GA handle to array 
      character*(*) integ_type  ! [input] Name of integrals to compute
c
c     local variables
c
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
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
      else
         write(6,*) ' integ_type = ', integ_type
         call errquit('int_1e_ga: unknown integ_type', 0)
      endif
c
c     Get info about the basis sets
c
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', ibas)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for jbas', jbas)
c
c     allocate necessary local temporary arrays on the stack
c
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c
c     k_* are the offsets corrsponding to the l_* handles
c
      call int_mem_1e(max1e, mem1)
      if (.not. MA_push_get(MT_DBL,max1e,'int_1e_ga:buf',l_buf,k_buf))
     $     call errquit('int_1e_ga: ma failed', max1e)
      if (.not. MA_push_get(MT_DBL, mem1,'int_1e_ga:scr',l_scr,k_scr))
     $     call errquit('int_1e_ga: ma failed', mem1)
c
c     Loop thru shells with static parallel work decomposition
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do ishell = 1, nshell_i
         if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) call errquit
     $        ('int_1e_ga: bas_cn2bfr failed ibas', ibas)
         do jshell = 1, nshell_j
            if (mod(ijshell, nproc) .eq. iproc) then
               if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $              call errquit('int_1e_ga: bas_cn2bfr ?', jbas)
               jdim = jhi - jlo + 1
c
c     Generate the integrals
c
               if (type .eq. 1) then
                  call int_1eke (ibas, jbas, ishell, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else if (type .eq. 2) then
                  call int_1epe (ibas, jbas, ishell, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else if (type .eq. 3) then
                  call int_1eov (ibas, jbas, ishell, jshell,
     $                 mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
               else
                  call errquit('int_1e_ga: invalid type?', type)
               endif
c
c     Add the integrals into the global array
c
               call ga_acc(g, jlo, jhi, ilo, ihi, dbl_mb(k_buf), jdim,
     $              1.0d0)
            endif
            ijshell = ijshell + 1
         enddo
      enddo
c
c     chop stack at first item allocated
c
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_1e_ga: pop failed', 0)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_1e_ga: pop failed', 0)
c
      call ga_sync()            ! So that no nasty races can result
c
      end
