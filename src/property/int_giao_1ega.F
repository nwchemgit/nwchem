      subroutine int_giao_1ega(ibas,jbas,g,integ_type,xyzpt,nat,
     &                         oskel)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "sym.fh"
#include "geom.fh"
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
      double precision xyzpt(3,*)   ! [input] coordinates of requested atoms (if needed)
c     
c     local variables
c     
      integer mone, two
      parameter (mone=-1, two=2)
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi, ilen, jlen
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer alo(3), ahi(3), ld(2)
      integer type
      logical odoit
      double precision q2
      integer nblocks
c
      logical odbug
      logical osome
c
      integer irtdb
      integer int_get_rtdb
      external int_get_rtdb
      integer iextbq
      integer nbq,nextbq,ncosbq,ntmp
      integer l_qbq,l_cbq
      integer k_qbq,k_cbq
      integer k_qextbq,k_cextbq
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
      elseif (inp_compare(.false., integ_type, 'srxRb')) then
         type = 2
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'l10')) then
         type = 3
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'tv10')) then
         type = 4
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'h01')) then
         type = 5
         nblocks = 3 * nat
      else if (inp_compare(.false., integ_type, 'h11 para'))then 
         type = 6 
         nblocks = 9 * nat
      else if (inp_compare(.false., integ_type, 'h11 dia'))then 
         type = 7 
         nblocks = 9 * nat
      else if (inp_compare(.false., integ_type, 'h11 all'))then 
         type = 8 
         nblocks = 9 * nat
      else if (inp_compare(.false., integ_type, 'dso'))then 
         type = 9 
         nblocks = 9 * nat ! nat is number of pairs in this case : nat * (nat-1)
      else if (inp_compare(.false., integ_type, 'pso'))then 
         type = 10        
         nblocks = 3 * nat
      else if (inp_compare(.false., integ_type, 'fc'))then 
         type = 11
         nblocks = nat
      else if (inp_compare(.false., integ_type, 'sd+fc'))then 
         type = 12
         nblocks = 6 * nat
      else if (inp_compare(.false., integ_type, 'velocity'))then 
         type = 13
         nblocks = 3      
      else if (inp_compare(.false., integ_type, 'angmom'))then 
         type = 14
         nblocks = 3      
      else if (inp_compare(.false., integ_type, 'bq10')) then
         type = 15
         nblocks = 3
c ... jochen: ops = pso-like integrals, for one nucleus
c        but all combinations of operators
      else if (inp_compare(.false., integ_type, 'ops')) then
         type = 16
         nblocks = 9 
      else
         write(6,*) ' integ_type = ', integ_type,ga_nodeid()
         call errquit('int_giao_1ega: unknown integ_type', 0, INT_ERR)
      end if

      if(osome) then
         write(6,*) 'in -int_giao_1ega- ... integ_type = ', type
      endif
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
      if (type.lt.9) then
         call int_mem_1e(max1e, mem1)  
         max1e = max1e*nblocks
      elseif (type.eq.9) then
         call int_init_dso(max1e,mem1,ibas,nat)
      elseif (type.eq.10) then
         call int_init_pso(max1e,mem1,ibas,nat)
      elseif (type.eq.16) then
         call int_init_ops(max1e,mem1,ibas)
      elseif (type.eq.11.or.type.eq.12) then
         call int_init_1eelec(max1e,mem1,ibas,2,nat)
      elseif (type.eq.15) then
         call int_mem_1e(max1e, mem1)  
         max1e = max1e*nblocks
      else
         call int_init_dip(max1e,mem1,ibas)
         max1e = max1e*nblocks
      endif
      mem1 = max(mem1,max1e)

      if(osome) then
         write(6,*) 'in -int_giao_1ega- max1e, mem1 = ', max1e, mem1
      endif
c
      if(.not.MA_push_get(MT_DBL,max1e,'int_giao_1ega:buf',l_buf,k_buf))
     $     call errquit('int_giao_1ega: ma failed', max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL, mem1,'int_giao_1ega:scr',l_scr,k_scr))
     $     call errquit('int_giao_1ega: ma failed', mem1, MA_ERR)
c
c     Get bq charges (external bq and cosmo)
c
      if (type.eq.15) then
        nbq=0
        nextbq = 0
        ncosbq = 0
        if(geom_extbq_on()) then
           nextbq = geom_extbq_ncenter()
           k_cextbq = geom_extbq_coord()
           k_qextbq = geom_extbq_charge()
           nbq = nextbq  ! external bq centers
        end if
c
c       Get rtdb handle
c
        irtdb = int_get_rtdb()   ! get rtdb handle
        if (rtdb_get(irtdb,'cosmo:nefc',mt_int,1,ncosbq))
     &       nbq = ncosbq  ! cosmo bq centers
c
        if (nextbq.gt.0.and.ncosbq.gt.0) 
     &       nbq = nextbq + ncosbq  ! all bq centers
c
c       Allocate memory for all the bq charges
c
        if(.not.MA_push_get(MT_DBL,nbq,'qbq',l_qbq,k_qbq))
     &    call errquit('int_giao_1ega: ma failed', l_qbq, MA_ERR)
        call dfill(nbq,0.d0,dbl_mb(k_qbq),1)
        if(.not.MA_push_get(MT_DBL,3*nbq,'cbq',l_cbq,k_cbq))
     $    call errquit('int_giao_1ega: ma failed', l_cbq, MA_ERR)
        call dfill(3*nbq,0.d0,dbl_mb(k_cbq),1)
c
c       Assign external bq charges and coordinates
c
        if (nextbq.gt.0) then
         call dcopy(nextbq,dbl_mb(k_qextbq),1,dbl_mb(k_qbq),1)
         call dcopy(3*nextbq,dbl_mb(k_cextbq),1,dbl_mb(k_cbq),1)
        end if ! nextbq 
c         
c       Get cosmo charges and coordinates
c
        if (ncosbq.gt.0) then
         call bq_fromrtdb(irtdb,'cosmo:efcz','cosmo:efcc',
     &       ntmp,dbl_mb(k_qbq+nextbq),dbl_mb(k_cbq+3*nextbq))
        end if ! ncosbq gt 0
      end if ! type .eq. 15
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
                  ilen = ihi-ilo+1
                  jlen = jhi-jlo+1
c     
c     Generate the integrals
c     
                  if (type .eq. 1) then      ! 3
                     call int_giaos10 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 2) then ! 3
                     call int_giaos100(jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 3) then ! 3
                     call int_giaol10 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 4) then ! 3
                     call int_giaotv10 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 5) then ! 3*nat
                     call int_giaoh01 (jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat)
                  else if (type .eq. 6) then ! 9*nat
                     call int_giaoh11 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf), 
     $                    xyzpt, nat, .true.,.false.)
                  else if (type .eq. 7) then ! 9*nat
                     call int_giaoh11 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat, .false.,.true.)
                  else if (type .eq. 8) then ! 9*nat
                     call int_giaoh11 (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat, .true.,.true.)
                  else if (type .eq. 9) then ! 9*nat*(nat-1)
                     call int_dso (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat)
                  else if (type .eq. 10) then ! 3*nat
                     call int_pso (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat)
                  else if (type .eq. 11) then ! nat
                     call int_1eelec(jbas, jshell, ibas,ishell,
     &                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     &                    mone,xyzpt,nat)
                  else if (type .eq. 12) then ! 6*nat
                     call int_1eelec(jbas, jshell, ibas,ishell,
     &                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     &                    two,xyzpt,nat)
                  else if (type .eq. 13) then ! 3
                     call int_veloc(jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 14) then ! 3
                     call int_angmom(jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
                  else if (type .eq. 15) then ! 3
                     call int_giaobq10(jbas, jshell, ibas, ishell, 
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     &                    dbl_mb(k_qbq),dbl_mb(k_cbq),nbq)
                  else if (type .eq. 16) then ! 9
                     call int_ops (jbas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf),
     $                    xyzpt, nat)
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
                  ld(1) = jlen
                  ld(2) = ilen
                  call nga_acc(g,alo,ahi,dbl_mb(k_buf),ld,1.0d0)
               end if
            endif
            ijshell = ijshell + 1
         end do
      end do

      if(osome) then
         write(6,*) 'in -int_giao_1ega- end of loops'
      endif
c     
c     chop stack at first item allocated
c     
      if (type.eq.15) then  ! for bq's
        if (.not. MA_pop_stack(l_cbq)) call errquit
     $     ('int_giao_1ega: pop failed', l_cbq, MA_ERR)
        if (.not. MA_pop_stack(l_qbq)) call errquit
     $     ('int_giao_1ega: pop failed', l_qbq, MA_ERR)
      end if ! type 15
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_giao_1ega: pop failed', l_scr, MA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_giao_1ega: pop failed', l_buf, MA_ERR)
     
      call ga_sync()            ! So that no nasty races can result
c     
      end
