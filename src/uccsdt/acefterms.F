      subroutine uccsdt_acefterms(basis)
c
c$Id: acefterms.F,v 1.1 2002-08-28 22:18:18 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
      integer qO_handles(8,0:7), qV_handles(8,0:7)
      integer symiu, list_nr
      integer basis, u_sh, u_block_size
c
c     Determine the blocking for SOs u in 
c     natural order ... and corresponding blocked range
c     alos get v block info as we need to transform this whole block at once
      call soints_shblock_info(basis, 1, sh_n_uq, vrange)
c     do blocks of SOs u (integer urange(2,0:7)
      do u_sh = 1, sh_n_uq, u_block_size
         u_sh_hi = u_sh+u_block_size
         call soints_shblock_info(basis, u_sh, u_sh_hi, urange)
         call cctrans(u_sh,u_sh_hi)
         call uccsdt_fterm(urange,vrange)
         call uccsdt_aterm(urange,vrange)
         do symiu = 0, nir-1
            do list_nr = 5, 8
               if (.not. ga_destroy(list(symiu,list_nr))) call errquit()
            end do
            do list_nr = 1, 8
               qO_handles(list_nr,symiu) = -1
               qV_handles(list_nr,symiu) = -1
            end do
         end do
         call uccsdt_cterm(urange,vrange,qO_handles,qV_handles)
         call uccsdt_eterm(qO_handles,qV_handles)
      end do
c
      end

