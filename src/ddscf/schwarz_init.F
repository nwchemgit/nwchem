      subroutine schwarz_init(geom, basis)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "cschwarz.fh"
#include "schwarz.fh"
#include "global.fh"
#include "geom.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer geom, basis       ! [input] handles
c
c     Initialize the schwarz inequality screening information
c     about this basis set.  Need to support these operations:
c
c     subroutine schwarz_init(basis)
c     subroutine schwarz_tidy()
c     double precision schwarz_max()
c     double precision schwarz_atom(iat, jat)
c     double precision schwarz_shell(ish, jsh)
c
      integer natoms, nshell, ish, jsh, 
     $     iat, jat, natri, nshtri, mem2, max2e,
     $     l_g, k_g, l_scr, k_scr, ishlo, ishhi, jshlo, jshhi,
     $     jshtop, idim, jdim, ish_bflo,  ish_bfhi, jsh_bflo, jsh_bfhi,
     $     numint, integ, len_csh, len_cat, ijat, me, nproc
      double precision at_max, sh_max
      logical status, odebug
      integer nxtask, next, task_size
      external nxtask
#include "itri.fh"
c
      odebug = .false.
c     
c     Get natoms, nshell, mem2, max2e
c     
      status = geom_ncent(geom, natoms)
      status = status .and. bas_numcont(basis, nshell)
      if (.not. status) call errquit('schwarz_init: info failed', 0,
     &       INPUT_ERR)
      call int_mem_2e4c(max2e, mem2)
c
c     Allocate space for the compressed arrays ... we shall have
c     1 entry every 4 bits and assume 8 bits per byte.
c     We also need to zero the arrays so values may be ORed in.
c
      nb_per_i = MA_sizeof(MT_INT, 1, MT_BYTE)
      natri = natoms*(natoms+1)/2
      len_cat = (natri - 1) / (2*nb_per_i) + 1
      if (.not. ma_alloc_get(MT_INT, len_cat, 'schwarz: atoms',
     $     l_cat, k_cat))
     $     call errquit('schwarz_init: ma_alloc 1', len_cat, MA_ERR)
      call ifill(len_cat, 0, int_mb(k_cat), 1)
c
      task_size = max(1,natri/(10*ga_nnodes()))
c
      nshtri = nshell*(nshell+1)/2
      len_csh = (nshtri - 1) / (2*nb_per_i) + 1
      if (.not. ma_alloc_get(MT_INT, len_csh, 'schwarz: shells',
     $     l_csh, k_csh))
     $     call errquit('schwarz_init: ma_alloc 2', len_csh, MA_ERR)

      call ifill(len_csh, 0, int_mb(k_csh), 1)
c     
c     scratch space needed for 2-electron routines
c     
      status = .true.
      status = status .and. 
     $     ma_push_get(MT_DBL, max2e, 'fock_2e: buf', l_g, k_g)
      status = status .and.
     $     ma_push_get(MT_DBL, mem2, 'fock_2e: scr', l_scr, k_scr)
      if (.not. status)
     $     call errquit('schwarz_init: insufficient memory', mem2,
     &       MA_ERR)
c
c     make sure that the integrals are being evaluated to sufficient precision
c
      call int_acc_high()
c     
      at_max = 0.0d0
      sch_max = 0.0d0
      ijat  = 0
      nproc = ga_nnodes()
      me    = ga_nodeid()
      next  = nxtask(nproc,task_size)
c
      do iat = 1, natoms
         if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $        call errquit('schwarz_init: bas_ce2cnr failed', ish,
     &       BASIS_ERR)
         do jat = 1, iat
            if (ijat .eq. next) then
               if (.not. bas_ce2cnr(basis, jat, jshlo, jshhi))
     $              call errquit('schwarz_init: bas_ce2cnr failed', jsh,
     &       BASIS_ERR)
               at_max = 0.0d0
               do ish = ishlo, ishhi
                  if (.not. bas_cn2bfr(basis, ish, ish_bflo, ish_bfhi))
     $                 call errquit('schwarz_init: bas_cn2bfr', ish,
     &       BASIS_ERR)
                  idim = ish_bfhi - ish_bflo + 1
                  jshtop = jshhi
                  if (iat .eq. jat) jshtop = ish
                  do jsh = jshlo, jshtop
                     if (.not. bas_cn2bfr(basis,jsh,jsh_bflo,jsh_bfhi))
     $                    call errquit('schwarz_init: bas_cn2bfr', jsh,
     &       BASIS_ERR)
                     jdim = jsh_bfhi - jsh_bflo + 1
c
                     call int_2e4c(basis, ish, jsh, basis, ish, jsh,
     $                    mem2, dbl_mb(k_scr), max2e, dbl_mb(k_g))
                     sh_max = 0.0d0
                     numint = idim*jdim*idim*jdim
* this for old max method which is not rotationally invariant
* so can break symmetry
*                     do integ = 0, numint-1
*                        sh_max = max(sh_max, abs(dbl_mb(k_g+integ)))
*                     enddo
*                     sh_max = sqrt(sh_max)
* this for new rms method that is invariant
                     do integ = 0, numint-1
                        sh_max = sh_max + dbl_mb(k_g+integ)**2
                     enddo
                     sh_max = sqrt(sqrt(sh_max))
c
                     at_max = max(at_max, sh_max)
c     
c     insert sh_max into the compressed shell array
c     
                     call comp4_insert(int_mb(k_csh), itri(ish,jsh),
     $                    sh_max, nb_per_i)
c     
                     if (odebug .and. sh_max.gt.0.0d0) then
                        write(6,1) ' shell ', ish, jsh, sh_max,
     $                       schwarz_shell(ish,jsh),
     $                       schwarz_shell(ish,jsh)/sh_max
                        call util_flush(6)
 1                      format(1x,a7,2i5,3d16.8)
                     endif
                  enddo
               enddo
               sch_max = max(sch_max, at_max)
c     
c     insert at_max into the compressed atom array
c     
               call comp4_insert(int_mb(k_cat), itri(iat,jat), at_max,
     $              nb_per_i)
c     
               if (odebug .and. at_max.gt.0.0d0) then
                  write(6,1) ' atoms ', iat, jat, at_max,
     $                 schwarz_atom(iat,jat),
     $                 schwarz_atom(iat,jat)/at_max
                  call util_flush(6)
               endif
               next = nxtask(nproc,task_size)
            endif
            ijat = ijat + 1
         enddo
      enddo
      next = nxtask(-nproc,task_size)
c     
c     A global OR is required to piece together the parts made in parallel
c
      call ga_igop(Msg_SchAt, int_mb(k_cat), len_cat, 'or')
      call ga_igop(Msg_SchSh, int_mb(k_csh), len_csh, 'or')
      call ga_dgop(Msg_SchMx, sch_max, 1, 'max')
c
      if (odebug .and. ga_nodeid().eq.0)
     $     write(6,*) ' sch_max ', sch_max
c
      if (.not. ma_pop_stack(l_scr)) call errquit('schwarz_init: pop',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_g)) call errquit('schwarz_init: pop', 1,
     &       MA_ERR)
c
      if (ga_nodeid().eq.0 .and. util_print('schwarz',print_high)) 
     $     then
         write(6,2) natri, nshtri
 2       format(/
     $        ' Schwarz integral screening summary '/
     $        ' ---------------------------------- '//
     $        '        No. of atom pairs =', i8/
     $        '       No. of shell pairs =', i8/)
c     
         write(6,3)
 3       format(' Atomic screening '/)
         call comp4_bins(int_mb(k_cat), natri, nb_per_i)
         write(6,4)
 4       format(/' Shell screening '/)
         call comp4_bins(int_mb(k_csh), nshtri, nb_per_i)
         write(6,*)
         call util_flush(6)
      endif
c
c     zero counts for screening information
c
      call ifill(nscreen, 0, iscreen, 1)
c
c     Restore normal accuracy for integral evaluation
c
      call int_acc_std()
c
      end
