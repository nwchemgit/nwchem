      subroutine int_1e_ga(ibas, jbas, g, integ_type, oskel)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
c     
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer ibas, jbas        ! [input] bra and ket basis sets
      integer g                 ! [output] GA handle to array 
      character*(*) integ_type  ! [input] Name of integrals to compute
      logical oskel             ! [input] If true generate symmetry unique list
c     
c     local variables
c     
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      logical odoit
      double precision q2
c     
      if (oscfps) call pstat_on(ps_int_1e)
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $        ('int_1e_ga: use of symmetry requires ibas=jbas', ibas)
      endif
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
      odoit = .true.
      q2 = 1.0d0
      do ishell = 1, nshell_i
         if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) call errquit
     $        ('int_1e_ga: bas_cn2bfr failed ibas', ibas)
         idim = ihi - ilo + 1
         do jshell = 1, nshell_j
c     
            if (oskel) 
     $           odoit = sym_shell_pair(ibas, ishell, jshell, q2)
c
            if (odoit) then
               if (mod(ijshell, nproc) .eq. iproc) then

*                  if (oskel)
*     $                 write(6,*) ga_nodeid(), ' getting bfr ', jshell
*                  call util_flush(6)

                  if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', jbas)
c     
c     Generate the integrals
c     

*                  if (oskel)
*     $                 write(6,*) ga_nodeid(), ' generating integrals '
*                  call util_flush(6)

                  if (type .eq. 1) then
                     call int_1eke (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 2) then
                     call int_1epe (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 3) then
                     call int_1eov (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else
                     call errquit('int_1e_ga: invalid type?', type)
                  endif
c     
c     Add the integrals into the global array
c     

*                  if (oskel)
*     $                 write(6,*) ga_nodeid(), ' ga_acc ', jshell, jlo,
*     $                 jhi
*                  call util_flush(6)

                  call ga_acc(g, ilo, ihi, jlo, jhi, dbl_mb(k_buf),
     $                 idim, q2)
               endif
               ijshell = ijshell + 1
            endif
            call ga_sync()
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
      if (oscfps) call pstat_off(ps_int_1e)
c     
      end
