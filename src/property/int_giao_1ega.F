      subroutine int_giao_1ega(ibas,jbas,g,integ_type,xyzpt,nat,
     &                         oskel)
C$Id: int_giao_1ega.F,v 1.1 2004-07-21 17:28:27 bert Exp $
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "sym.fh"
c     
c     Compute the desired type of 1e GIAO integrals 
c     and ADD them into the given global array.
c     This version computes the full square of integrals and should work
c     OK even if ibas != jbas.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer ibas, jbas            ! [input] bra and ket basis sets
      integer g                     ! [output] GA handle to array, one for each field direction (if needed)
      character*(*) integ_type      ! [input] Name of integrals to compute
      logical oskel                 ! [input] If true generate symmetry unique list
      integer nat                   ! [input] number of atoms for which we get the integrals (if needed)
      double precision xyzpt(3,nat) ! [input] coordinates of requested atoms (if needed)
c     
c     local variables
c     
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer alo(3), ahi(3), ld(2)
      integer type
      logical odoit
      double precision q2
      integer i, nblocks
c
      logical odbug
      logical osome
c
      odbug=.false.
      osome=.false.
      osome=osome.or.odbug
      odbug=odbug.and.(ga_nodeid().eq.0)
      osome=osome.and.(ga_nodeid().eq.0)
      if(osome) then
         write(6,*) 'in -int_giao_1ega- ... integ_type = ',
     $              integ_type,ga_nodeid()
         call util_flush(6)
      endif
c
      call ga_sync()
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $      ('int_giao_1ega: use of symmetry requires ibas=jbas', ibas,
     &       BASIS_ERR)
      end if
c     
      if (inp_compare(.false., integ_type, 's10')) then
         type = 1
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'l10')) then
         type = 2
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'tv10')) then
         type = 3
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'h01')) then
         type = 4
         nblocks = 3 * nat
      else if (inp_compare(.false., integ_type, 'h11 para'))then 
         type = 5 
         nblocks = 9 * nat
      else if (inp_compare(.false., integ_type, 'h11 dia'))then 
         type = 6 
         nblocks = 9 * nat
      else if (inp_compare(.false., integ_type, 'h11 all'))then 
         type = 7 
         nblocks = 9 * nat
      else
         write(6,*) ' integ_type = ', integ_type,ga_nodeid()
         call errquit('int_giao_1ega: unknown integ_type', 0, INT_ERR)
      end if
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell_i)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      if (.not. bas_numcont(jbas, nshell_j)) call errquit
     $     ('rhf_fock_1e: bas_numcont failed for jbas', jbas,
     &       BASIS_ERR)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
c
      call int_mem_1e(max1e, mem1)  
      max1e = max1e*nblocks
c
      if(.not.MA_push_get(MT_DBL,max1e,'int_giao_1ega:buf',l_buf,k_buf))
     $     call errquit('int_giao_1ega: ma failed', max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL, mem1,'int_giao_1ega:scr',l_scr,k_scr))
     $     call errquit('int_giao_1ega: ma failed', mem1, MA_ERR)
c     
c     Loop thru shells with static parallel work decomposition
c     
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      q2 = 1.0d0
      do jshell = 1, nshell_j
         do ishell = 1, nshell_i
c
            if (mod(ijshell, nproc) .eq. iproc) then
               odoit = .true.
               if (oskel) 
     $              odoit = sym_shell_pair(ibas, ishell, jshell, q2)
c
               if (odoit) then
                  if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
                  if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
c     
c     Generate the integrals
c     
                  call dfill(max1e,0.0d0,dbl_mb(k_buf),1)
                  if (type .eq. 1) then      ! 3
                     call int_giaos10 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 2) then ! 3
                     call int_giaol10 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 3) then ! 3
                     call int_giaotv10 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 4) then ! 3*nat
                     call int_giaoh01 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat)
                  else if (type .eq. 5) then ! 9*nat
                     call int_giaoh11 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf), 
     $                    xyzpt, nat, .true.,.false.)
                  else if (type .eq. 6) then ! 9*nat
                     call int_giaoh11 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat, .false.,.true.)
                  else if (type .eq. 7) then ! 9*nat
                     call int_giaoh11 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat, .true.,.true.)
                  else
                     call errquit('int_giao_1ega: invalid type?', type,
     &       GA_ERR)
                  end if
c     
c     Add the integrals into the global array
c    
                  alo(1) = jlo
                  ahi(1) = jhi
                  alo(2) = ilo
                  ahi(2) = ihi
                  alo(3) = 1
                  ahi(3) = nblocks
                  ld(1) = jhi-jlo+1
                  ld(2) = ihi-ilo+1
                  call nga_acc(g,alo,ahi,dbl_mb(k_buf),ld,1.0d0)
               end if
            endif
            ijshell = ijshell + 1
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_giao_1ega: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_giao_1ega: pop failed', 0, GA_ERR)
c     
      call ga_sync()            ! So that no nasty races can result
c     
      end
