      subroutine uccsdt_acefterms()
c
c$Id: acefterms.F,v 1.3 2002-09-18 19:04:51 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
      integer qO_handles(8,0:7), qV_handles(8,0:7)
      integer lists(0:7,14)
      integer symiu, list_nr
      integer u_sh, u_block_size, u_sh_hi
c
c     Allocate qO_handles and qV_handles
c
      call qhandles_allocate(qO_handles,qV_handles)
c
c     Determine the blocking for SOs u in 
c     natural order ... and corresponding blocked range
c     also get v block info as we need to transform this whole block at once
c
      call soints_shblock_info(basis, 1, sh_n_uq, vrange)
c
c     do blocks of SOs u (integer urange(2,0:7)
c
      do u_sh = 1, sh_n_uq, u_block_size
         u_sh_hi = u_sh+u_block_size
         call soints_shblock_info(basis, u_sh, u_sh_hi, urange)
         call cctrans(u_sh,u_sh_hi,lists)
         call uccsdt_fterm(urange,vrange)
         call uccsdt_aterm(urange,vrange)
         do symiu = 0, nir-1
            do list_nr = 5, 8
               if (.not. ga_destroy(lists(symiu,list_nr))) call 
     &            errquit('uccsdt: acefterms ga_destroy failed',0)
            end do
         end do
         call uccsdt_cterm(urange,vrange,qO_handles,qV_handles)
      end do
      call uccsdt_eterm(qO_handles,qV_handles)
c
      end
c
      subroutine qhandles_allocate(qO_handles,qV_handles)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"

      integer qO_handles(8,0:7), qV_handles(8,0:7)
      integer symiu
c
c     do 2 * 8 * nir allocations of fock matrices of appropriate sizes
c     hence for each symmetry block:
c        get # of occ and # of virt
c        qO is dimension occ*occ
c        qV is dimension virt*virt
c        store handle in array
c
      do symiu = 0, nir-1
         if (.not. ga_create(MT_DBL,no_sym(symiu,1)*no_sym(symiu,1),
     &      'qO 1',0,0,qO_handles(1,symiu))) call 
     &      errquit('uccsd: failed to allocate qO 1', 0)
         if (.not. ga_create(MT_DBL,no_sym(symiu,1)*no_sym(symiu,1),
     &      'qO 2',0,0,qO_handles(2,symiu))) call 
     &      errquit('uccsd: failed to allocate qO 2', 0)
         if (.not. ga_create(MT_DBL,no_sym(symiu,2)*no_sym(symiu,2),
     &      'qO 4',0,0,qO_handles(4,symiu))) call 
     &      errquit('uccsd: failed to allocate qO 4', 0)
         if (.not. ga_create(MT_DBL,no_sym(symiu,2)*no_sym(symiu,2),
     &      'qO 5',0,0,qO_handles(5,symiu))) call 
     &      errquit('uccsd: failed to allocate qO 5', 0)
         if (.not. ga_create(MT_DBL,no_sym(symiu,2)*no_sym(symiu,2),
     &      'qO 8',0,0,qO_handles(8,symiu))) call 
     &      errquit('uccsd: failed to allocate qO 8', 0)
         if (.not. ga_create(MT_DBL,nv_sym(symiu,1)*nv_sym(symiu,1),
     &      'qV 1',0,0,qV_handles(1,symiu))) call 
     &      errquit('uccsd: failed to allocate qV 1', 0)
         if (.not. ga_create(MT_DBL,nv_sym(symiu,1)*nv_sym(symiu,1),
     &      'qV 2',0,0,qV_handles(2,symiu))) call 
     &      errquit('uccsd: failed to allocate qV 2', 0)
         if (.not. ga_create(MT_DBL,nv_sym(symiu,2)*nv_sym(symiu,2),
     &      'qV 4',0,0,qV_handles(4,symiu))) call 
     &      errquit('uccsd: failed to allocate qV 4', 0)
         if (.not. ga_create(MT_DBL,nv_sym(symiu,2)*nv_sym(symiu,2),
     &      'qV 5',0,0,qV_handles(5,symiu))) call 
     &      errquit('uccsd: failed to allocate qV 5', 0)
         if (.not. ga_create(MT_DBL,nv_sym(symiu,2)*nv_sym(symiu,2),
     &      'qV 8',0,0,qV_handles(8,symiu))) call 
     &      errquit('uccsd: failed to allocate qV 8', 0)
      enddo
c
      end
