      logical function int_norm_2c(rtdb,basisin)
c $Id$
c::assumptions
*..   MA_init has been called
c
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "rtdb.fh"
      integer rtdb ! [in]
      integer basisin             ! [input] basis handle
      integer basis, nucont, iucont, type, nprim, ngen, ept, cpt, igen
      integer maxg2c, mem2c, lscr, h_scr, i_scr, maxang
      logical norm_ok
c
#include "bas_exndcf_sfn.fh"
c
c     Normalize the contraction coefficients in this basis.
c     Eventually we need to store both the normalized and unnormalized
c     coefficients so that the user can print out something recognizable
c
      int_norm_2c = bas_check_handle(basisin, 'int_norm_2c')
      if (.not. int_norm_2c) return

c
c check and set basis normalization
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
      norm_ok = bas_norm_id(basis).eq.BasNorm_2c
      if (norm_ok) return ! basis has been normalized properly
      norm_ok = bas_norm_id(basis).eq.BasNorm_UN
      if (.not.norm_ok) then
        if (.not. bas_norm_print(basisin)) call errquit
     &      ('int_norm_2c: bas_norm_print failed',911, BASIS_ERR)
        write(luout,*)' basis has already been normalized somehow'
        call errquit('int_norm_2c: basis normalization error ',911,
     &       BASIS_ERR)
      endif
      bas_norm_id(basis) = BasNorm_2c
c
c...  get memory requirements for one electron integrals
c
      if (.not. rtdb_get(rtdb,'int:mem2emax',Mt_Int,1,mem2c))
     &    mem2c = 14 000 000
      call emem_2e2c(basisin, basisin, maxg2c, mem2c)
      lscr = maxg2c + mem2c
      int_norm_2c =
     &    MA_push_get(mt_dbl,lscr,'int_norm_2c:scr',h_scr,i_scr)
      call dcopy(lscr,0.0d00,0,dbl_mb(i_scr),1)
c
c... define xyz for angular momentum      
      if (.not.bas_high_angular(basisin,maxang)) call errquit
     &    ('int_norm_2c: bas_high_angular failed?',911, BASIS_ERR)
      call defNxyz(maxang)
c
      basis = basisin + BASIS_HANDLE_OFFSET
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
              call nmcoef2c(dbl_mb(mb_exndcf(ept, basis)),
     &              dbl_mb(mb_exndcf(cpt+igen*nprim,basis)),
     &              type, nprim, dbl_mb(i_scr), lscr)
            enddo
c....... sp
         else if (type .eq. -1) then
            call nmcoef2c(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt, basis)), 
     &           0, nprim, dbl_mb(i_scr), lscr)
            call nmcoef2c(dbl_mb(mb_exndcf(ept, basis)),
     &            dbl_mb(mb_exndcf(cpt+nprim,basis)),
     &           1, nprim, dbl_mb(i_scr), lscr)
         else if (type .eq. -2) then
c........ spd
            call nmcoef2c(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt, basis)), 
     &           0, nprim, dbl_mb(i_scr), lscr)
            call nmcoef2c(dbl_mb(mb_exndcf(ept, basis)),
     &            dbl_mb(mb_exndcf(cpt+nprim,basis)),
     &           1, nprim, dbl_mb(i_scr), lscr)
            call nmcoef2c(dbl_mb(mb_exndcf(ept, basis)),
     &            dbl_mb(mb_exndcf(cpt+2*nprim,basis)),
     &           2, nprim, dbl_mb(i_scr), lscr)
         else
            call errquit('int_norm_2c: invalid contraction type', 0,
     &       BASIS_ERR)
         endif
      enddo
c
c... free memory used for 2c2e integrals
      int_norm_2c = MA_pop_stack(h_scr)
c
      end
