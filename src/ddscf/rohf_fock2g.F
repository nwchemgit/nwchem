      subroutine rohf_get_grad(g_grad)
      implicit none
#include "cscf.fh"
#include "crohf.fh"
      integer g_grad
      call rohf_fock2grad(nbf, nclosed, nopen, crohf_g_fcv, 
     $     crohf_g_fpv, crohf_g_fcp, g_grad)
      end
      subroutine rohf_fock2grad( nbf, nclosed, nopen,
     $     g_fcv, g_fpv, g_fcp,
     $     g_grad)
C     $Id: rohf_fock2g.F,v 1.1 1995-12-13 01:25:43 d3g681 Exp $
      implicit none
#include "global.fh"
      integer nbf, nclosed, nopen ! [input] Occupation
      integer g_fcv, g_fpv, g_fcp ! [input] Fock matrices for each subspace
      integer g_grad            ! [output] Gradient vector
c     
      integer nvir, voff, ooff, oend, xoff, xend
c     
      nvir = nbf - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c     
      call ga_zero(g_grad)
      xend = nclosed*nvir
      call ga_copy_patch('n', g_fcv, voff, nbf, 1, nclosed,
     $     g_grad, 1, xend, 1, 1 )
      call ga_dscal_patch(g_grad, 1, xend, 1, 1, 4.d0 )
      if (nopen .gt. 0) then
         xoff = xend + 1
         xend = xend + nopen*nvir
         call ga_copy_patch('n', g_fpv, voff, nbf, ooff, oend,
     $        g_grad, xoff, xend, 1, 1 )
         call ga_dscal_patch(g_grad, xoff, xend, 1, 1, 2.d0 )
         xoff = xend + 1
         xend = xend + nopen*nclosed
         call ga_copy_patch('n', g_fcp, ooff, oend, 1, nclosed,
     $        g_grad, xoff, xend, 1, 1 )
         call ga_dscal_patch(g_grad, xoff, xend, 1, 1, 2.d0 )
      endif
c
      call ga_screen(g_grad, 1d-13)
c     
      end
