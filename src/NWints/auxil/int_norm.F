      logical function int_normalize(basisin)
c $Id: int_norm.F,v 1.5 1999-07-14 00:28:45 d3e129 Exp $
      implicit none
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
      integer basisin             ! [input] basis handle
      integer basis, nucont, iucont, type, nprim, ngen, ept, cpt, igen
      logical norm_ok
c
#include "bas_exndcf_sfn.fh"
c
c     Normalize the contraction coefficients in this basis.
c     Eventually we need to store both the normalized and unnormalized
c     coefficients so that the user can print out something recognizable
c
      int_normalize = bas_check_handle(basisin, 'int_normalize')
      if (.not. int_normalize) return
c
      basis = basisin + BASIS_HANDLE_OFFSET

c
c check and set basis normalization
c
      norm_ok = bas_norm_id(basis).eq.BasNorm_STD
      if (norm_ok) return ! basis has been normalized properly
      norm_ok = bas_norm_id(basis).eq.BasNorm_UN
      if (.not.norm_ok) then
        if (.not. bas_norm_print(basisin)) call errquit
     &      ('int_normalize: bas_norm_print failed',911)
        write(luout,*)' basis has already been normalized somehow'
        call errquit('int_norm: basis normalization error ',911)
      endif
      bas_norm_id(basis) = BasNorm_STD
c
      nucont = infbs_head(HEAD_NCONT,basis)
      do iucont = 1, nucont
         type  = infbs_cont(CONT_TYPE, iucont,basis)
         nprim = infbs_cont(CONT_NPRIM,iucont,basis) 
         ngen  = infbs_cont(CONT_NGEN, iucont,basis)
         ept   = infbs_cont(CONT_IEXP, iucont,basis)
         cpt   = infbs_cont(CONT_ICFP, iucont,basis)
         if (type .ge. 0) then
            do igen = 0,ngen-1
               call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     $              dbl_mb(mb_exndcf(cpt+igen*nprim,basis)),
     &              type, nprim, 1)
            enddo
c....... sp
         else if (type .eq. -1) then
           call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt, basis)),0, nprim, 1) 
           call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt+nprim,basis)), 1, nprim, 1)
         else if (type .eq. -2) then
c........ spd
           call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt, basis)), 0, nprim, 1) 
            call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt+nprim,basis)),
     $           1, nprim, 1)
            call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &            dbl_mb(mb_exndcf(cpt+2*nprim,basis)),
     $           2, nprim, 1)
         else
            call errquit('int_normalize: invalid contraction type', 0)
         endif
      enddo
c
      end
