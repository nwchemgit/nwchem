      subroutine uccsdt_acefterms(x_lists)
c
c$Id$
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
      integer qO_handles(0:7,8), qV_handles(0:7,8)
      integer lists(0:7,14), x_lists(3,0:7)
      integer symiu, list_nr
      integer u_sh, u_block_size, u_sh_hi
c
c     Allocate qO_handles and qV_handles
c
      call qhandles_allocate(qO_handles,qV_handles)
c
c     Get the correct scaling for the X-terms
c
      call uccsdt_aterm_init(x_lists)
c
c     Determine the blocking for SOs u in 
c     natural order ... and corresponding blocked range
c     also get v block info as we need to transform this whole block at once
c
      call soints_shblock_info(basis, 1, sh_n_uq, vrange)
      call init_iv_block_info(vrange)
c
c     do blocks of SOs u (integer urange(2,0:7)
c
      do u_sh = 1, sh_n_uq, u_block_size
         u_sh_hi = u_sh+u_block_size
         call soints_shblock_info(basis, u_sh, u_sh_hi, urange)
         call init_iu_block_info(urange)
         call cctrans(u_sh,u_sh_hi,lists)
         call uccsdt_fterm(urange,vrange,lists)
         call uccsdt_aterm_acc(urange,vrange,lists)
         call uccsdt_cterm(urange,vrange,qO_handles,qV_handles)
      end do
      call uccsdt_aterm_finish(x_lists)
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

      integer qO_handles(0:7,8), qV_handles(0:7,8)
      integer symiu, handles, ih, is
      integer nr_handle(5), spin_handle(5)
      data nr_handle /1,2,4,5,8/
      data spin_handle /1,1,2,2,2/
c
c     do 2 * 8 * nir allocations of fock matrices of appropriate sizes
c     hence for each symmetry block:
c        get # of occ and # of virt
c        qO is dimension virt*virt
c        qV is dimension occ*occ   
c        store handle in array
c
      do symiu = 0, nir-1
         do handles = 1, 5
           ih = nr_handle(handles)
           is = spin_handle(handles)
           if (.not. ga_create(MT_DBL,no_sym(symiu,is)*no_sym(symiu,is),
     &        'qO',0,0,qO_handles(symiu,ih))) call 
     &        errquit('uccsd: failed to allocate qO',ih)
           if (.not. ga_create(MT_DBL,nv_sym(symiu,is)*nv_sym(symiu,is),
     &        'qV',0,0,qV_handles(symiu,ih))) call 
     &        errquit('uccsd: failed to allocate qV',ih)
         enddo
      enddo
c
      end
c
      subroutine init_iv_block_info(vrange)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer vrange(2,0:7)
      integer spini,symiv,symv,symi,iv,ind
c
      do spini = 1, 2
         do symiv = 0, nir-1
            ind = 0
            do symv = 0, nir-1
               symi = ieor(symv,symiv)
               do iv = vrange(1,symv), vrange(2,symv)
                  oso_v_off(iv,symiv,spini) = ind
                  ind = ind + no_sym(symi,spini)
               enddo
            end do
            oso_v_len(symiv,spini) = ind
         end do
      end do
c
      end
c
      subroutine init_iu_block_info(urange)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7)
      integer spini,symiu,symu,symi,iu,ind
c
      do spini = 1, 2
         do symiu = 0, nir-1
            ind = 0
            do symu = 0, nir-1
               symi = ieor(symu,symiu)
               do iu = urange(1,symu), urange(2,symu)
                  oso_u_off(iu,symiu,spini) = ind
                  ind = ind + no_sym(symi,spini)
               enddo
            end do
            oso_u_len(symiu,spini) = ind
         end do
      end do
c
      end
