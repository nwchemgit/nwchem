      subroutine int_mpole_1ega(ibas,jbas,g,integ_type,xyzpt,oskel)
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
c     
c     Compute the desired type of 1e multipole integrals 
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
      double precision xyzpt(3)     ! [input] coordinates of requested atoms (if needed)
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
      integer nblocks, mtpint, idim, offset, jdim
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
         write(6,*) 'in -int_mpole_1ega- ... integ_type = ',
     $              integ_type,ga_nodeid()
         call util_flush(6)
      endif
c
      call ga_sync()
c     
      if (oskel) then
         if (ibas.ne.jbas) call errquit
     $      ('int_mpole_1ega: use of symmetry requires ibas=jbas', ibas,
     &       BASIS_ERR)
      end if
c
      if (inp_compare(.false., integ_type, 'dipole'))then 
         type = 1
         nblocks = 3
      else if (inp_compare(.false., integ_type, 'quadrupole'))then 
         type = 2
         nblocks = 6
      else if (inp_compare(.false., integ_type, 'octupole'))then
         type = 3
         nblocks = 10
      else
         write(6,*) ' integ_type = ', integ_type,ga_nodeid()
         call errquit('int_mpole_1ega: unknown integ_type', 0, INT_ERR)
         nblocks = 0 ! meaningless, but makes silences gfortran error message
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
c     There might be something wrong here because it's not getting enough buf/scr for larger multipoles with smaller basis sets
c
      call int_init_1emtpl(max1e,mem1,ibas,type)
      mem1 = max(mem1,max1e)
c
      if(.not.MA_push_get(MT_DBL,max1e,'int_mpl_1ega:buf',l_buf,k_buf))
     $     call errquit('int_mpole_1ega: ma failed', max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL, mem1,'int_mpl_1ega:scr',l_scr,k_scr))
     $     call errquit('int_mpole_1ega: ma failed', mem1, MA_ERR)
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
     &                              BASIS_ERR)
                  if (.not. bas_cn2bfr(jbas, jshell, jlo, jhi))
     $                 call errquit('int_1e_ga: bas_cn2bfr ?', jbas,
     &                              BASIS_ERR)
c
                  ilen = ihi-ilo+1
                  jlen = jhi-jlo+1
c     
c     Generate the integrals
c     
c                  if (type .eq. 1) then      ! 3
c                     call int_mpolel(jbas,jshell,ibas,ishell,1,xyzpt,
c     &                    mem1,dbl_mb(k_scr),max1e,dbl_mb(k_buf),mtpint)
c                  else
c                     call errquit('int_mpl_1ega: invalid type?', type,
c     &                            GA_ERR)
c                  end if
c
c                  No type logic is necessary and it wastes cycles.
c                  The type has already been checked and an error is called if type>3.
c
                   call int_mpolel(jbas,jshell,ibas,ishell,type,xyzpt,
     &                  mem1,dbl_mb(k_scr),max1e,dbl_mb(k_buf),mtpint)
c     
c     Add the integrals into the global array
c    
                  offset = 0
                  do jdim = jlo, jhi
                     do idim = ilo, ihi
                        alo(1) = idim
                        ahi(1) = idim
                        alo(2) = jdim
                        ahi(2) = jdim
                        alo(3) = 1
                        ahi(3) = (type+1)*(type+2)/2
                        ld(1) = 1   
                        ld(2) = 1   
                        call nga_acc(g,alo,ahi,dbl_mb(k_buf+offset),
     &                               ld,1.0d0)
                        offset = offset + nblocks
                     enddo
                  enddo
               end if
            endif
            ijshell = ijshell + 1
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_mpole_1ega: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_mpole_1ega: pop failed', 0, GA_ERR)
c     
      call ga_sync()            ! So that no nasty races can result
c     
      end
