      logical function int_normalize(basisin)
c $Id: int_norm.F,v 1.1 1995-12-13 01:34:07 d3g681 Exp $
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
      integer basisin             ! [input] basis handle
      integer basis, nucont, iucont, type, nprim, ngen, ept, cpt, igen
      logical bas_check_handle
      external bas_check_handle
c
c     Normalize the contraction coefficients in this basis.
c     Eventually we need to store both the normalized and unnormalized
c     coefficients so that the user can print out something recognizable
c
      int_normalize = bas_check_handle(basisin, 'int_normalize')
      if (.not. int_normalize) return
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
               call nmcoef(exndcf(ept, basis),
     $              exndcf(cpt+igen*nprim,basis), type, nprim, 1)
            enddo
c....... sp
         else if (type .eq. -1) then
            call nmcoef(exndcf(ept, basis), exndcf(cpt, basis), 
     $           0, nprim, 1)
            call nmcoef(exndcf(ept, basis), exndcf(cpt+nprim,basis),
     $           1, nprim, 1)
         else if (type .eq. -2) then
c........ spd
            call nmcoef(exndcf(ept, basis), exndcf(cpt, basis), 
     $           0, nprim, 1)
            call nmcoef(exndcf(ept, basis), exndcf(cpt+nprim,basis),
     $           1, nprim, 1)
            call nmcoef(exndcf(ept, basis), exndcf(cpt+2*nprim,basis),
     $           2, nprim, 1)
         else
            call errquit('int_normalize: invalid contraction type', 0)
         endif
      enddo
c
      end
