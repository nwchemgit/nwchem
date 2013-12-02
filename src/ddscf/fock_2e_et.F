      subroutine fock_2e_est_task(geom, oskel, tol2e, task_size)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "geom.fh"
#include "global.fh"
#include "util.fh"
      integer geom
      logical oskel
      double precision tol2e   
      integer task_size      ! [output] returns recommended task size
c
c     Determine the appropriate task chunk size for the gradients.
c     [The fock build now use blocking]
c     For load imbalance to be about 1% we need about 100 tasks per
c     processor.  Many more than this is unecessary and may bottleneck
c     the shared counter.  Cut back to 20 tasks for experiment.
c
c     The default values of -1 is used to check for user input
c     
      integer nat, iat, jat
      integer npair, npair_sym, ntask
      double precision q2
c
      if (.not. geom_ncent(geom, nat)) call errquit
     $     ('fock_2e_est_task: geom bad?',0, GEOM_ERR)
c
      npair = 0
      npair_sym = 0
      do iat = 1, nat
         do jat = 1, iat
            if (schwarz_atom(iat, jat) .ge. tol2e) then
               npair = npair + 1
               if (oskel) then
                  if (sym_atom_pair(geom, iat, jat, q2)) 
     $                 npair_sym = npair_sym + 1
               endif
            endif
         enddo
      enddo
c
      if (.not. oskel) npair_sym = npair
      ntask = npair*npair_sym / (2*ga_nnodes())
c
      task_size = ntask
      if (ga_nnodes() .gt. 1) task_size = max(task_size/20,1)
      if (ga_nodeid().eq.0 .and.
     $     util_print('task_size',print_high)) then
         write(6,*) ' ntask, task_size ', ntask, task_size
         call util_flush(6)
      end if
c
      end
      subroutine fock_2e_block_atoms(basis, oskel, tol2e,
     $     blocks, nblock, maxblock)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "cfock.fh"
      integer basis
      logical oskel
      double precision tol2e   
      integer blocks(2,*)       ! [output] Blocking of atoms
      integer nblock            ! [output] No. of blocks
      integer maxblock          ! [output] Max bf in a block
c
      integer geom
      integer nat, nbf, task_size
      integer l_cnr, k_cnr, l_bfr, k_bfr
c
c     Get mapping from centers to shells/functions
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('fock_2e_block_atoms: bad basis handle?', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, nat)) call errquit
     $     ('fock_2e_block_atoms: geom bad?',0, GEOM_ERR)
      if (.not. ma_push_get(MT_INT, 2*nat, 'f2eba', l_cnr, k_cnr))
     $     call errquit('fock_2e_block: ma ', 2*nat, MA_ERR)
      if (.not. ma_push_get(MT_INT, 2*nat, 'f2eba', l_bfr, k_bfr))
     $     call errquit('fock_2e_block: ma ', 2*nat, MA_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('fock_2e_block_atoms: bad basis handle?', 0, BASIS_ERR)
c
      call bas_all_ce_info(basis, int_mb(k_cnr), int_mb(k_bfr))
c
c     Figure out mean density blocksize from the tasksize in 
c     atom quartets ... this is the task size in atom quartets
c     necessary to get about 100 tasks per process, taking into
c     account sparsity and symmetry.  Constrain this so that
c     a block is not larger than 128 so that we do not use 
c     too much memory.  Have a profound problem if there are
c     more than 255 functions on a single atom.
c
c     Also have the limit of 10 atoms per block.
c
      task_bf = min(task_bf,128,nbf)
      if (task_bf.eq.-1) then
         call fock_2e_est_task(geom, oskel, tol2e, task_size)
         task_bf = max(1,nint(nbf*sqrt(sqrt(dble(task_size)))/nat))
         task_bf = min(task_bf, 128)
      endif
c
      call fock_2e_do_block(basis, geom, nat,
     $     int_mb(k_bfr), task_bf, blocks, nblock, maxblock)
c
      if (.not. ma_pop_stack(l_bfr)) call errquit
     $     ('fock_2e_block_atoms: bad ma?',0, MA_ERR)
      if (.not. ma_pop_stack(l_cnr)) call errquit
     $     ('fock_2e_block_atoms: bad ma?',0, MA_ERR)
c
      end
      subroutine fock_2e_do_block(basis, geom, nat, bfr,
     $     task_bf, blocks, nblock, maxblock)
      implicit none
c
c     This subroutine groups atoms into blocks such that
c     1. no block contains more than 10 atoms
c     2. no block contains more than task_bf basis functions
c     3. no block starts or ends with atoms without basis functions
c        (e.g. point charges)
c     4. all blocks contain at least 1 basis function
c     The resulting blocks are the basis for the work distribution
c     in the distributed data Fock build.
c
c     The subroutine returns:
c     - blocks   -- the array of blocks such that blocks(1,i) is the
c                   atom of lowest rank and blocks(2,i) is the atom of
c                   highest rank
c     - nblock   -- the number of blocks
c     - maxblock -- the maximum number of basis functions in any block
c
c     Ad 3. This condition on the blocks has to be imposed because of
c           an inconsistency of the data in the array bfr. This array
c           has for any atom i bfr(1,i) the basis function of lowest
c           rank, and bfr(2,i) the basis function of highest rank on
c           that atom. However, if atom i has no basis functions then
c           bfr(1,i) = 0, and bfr(2,i) = -1, i.e. these are arbitrary
c           values. The result is that the number of basis functions
c           in a block cannot be established correctly by looking
c           at the first and last atom of the block if condition
c           3 is not met.
c
#include "global.fh"
#include "util.fh"
      integer basis
      integer geom
      integer nat
      integer bfr(2,nat)
      integer task_bf
      integer blocks(2,*), nblock, maxblock, minblock
c
      integer iat, n, iatlo, i, na, nnext
      logical firstprint
      data firstprint /.true./  ! So that only print once per job
c
      nblock = 0
      n = 0                     ! No. of bf in block
      na = 0                    ! No. of at in block
      iatlo = 1
      do iat = 1, nat
         na = na + 1
         n = n + bfr(2,iat) - bfr(1,iat) + 1
         nnext = n
         if (iat.lt.nat) nnext = n + bfr(2,iat+1) - bfr(1,iat+1) + 1
         if (nnext.gt.task_bf .or. na.ge.10) then
            nblock = nblock + 1
            blocks(1,nblock) = iatlo
            blocks(2,nblock) = iat
            iatlo = iat + 1
            n = 0
            na= 0
         endif
      enddo
      if (n .gt. 0) then
         nblock = nblock + 1
         blocks(1,nblock) = iatlo
         blocks(2,nblock) = nat
      endif
c
      n = 0
      do i = 1, nblock
c
c       Remove any leading atoms without basis functions
c
        do while (bfr(2,blocks(1,i))-bfr(1,blocks(1,i))+1 .le. 0 .and.
     $            blocks(1,i).lt.blocks(2,i))
          blocks(1,i) = blocks(1,i) + 1
        enddo
c
c       Remove any trailing atoms without basis functions
c
        do while (bfr(2,blocks(2,i))-bfr(1,blocks(2,i))+1 .le. 0 .and.
     $            blocks(2,i).gt.blocks(1,i))
          blocks(2,i) = blocks(2,i) - 1
        enddo
c
c       Keep a block only if it has at least one basis function
c
        if (bfr(2,blocks(2,i))-bfr(1,blocks(2,i))+1 .gt. 0) then
          n = n + 1
          blocks(1,n) = blocks(1,i)
          blocks(2,n) = blocks(2,i)
        endif
c
      enddo
      nblock = n
c
      maxblock = 0
      minblock = task_bf
      do i = 1, nblock
         maxblock = max(maxblock,
     $        bfr(2,blocks(2,i))-bfr(1,blocks(1,i))+1)
         minblock = min(minblock,
     $        bfr(2,blocks(2,i))-bfr(1,blocks(1,i))+1)
      enddo
c
      if (ga_nodeid().eq.0 .and. util_print('tasks',print_high)
     $     .and. firstprint) then
         firstprint = .false.
         write(6,1) minblock, maxblock, nblock
         write(6,2) (blocks(1,i), blocks(2,i),i=1, nblock)
         write(6,*)
         call util_flush(6)
 1       format(/' Fock matrix atom blocking: min nbf per block: ',i4/
     $           '                            max nbf per block: ',i4/
     $           '                                no. of blocks: ',i4/)
 2       format(5(:'(',i4,':',i4,')', 3x))
      end if
c
      end
