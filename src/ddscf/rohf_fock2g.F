      subroutine rohf_get_grad(g_grad)
      implicit none
#include "cscf.fh"
#include "crohf.fh"
      integer g_grad
      call rohf_fock2grad(nmo, nclosed, nopen, crohf_g_fcv, 
     $     crohf_g_fpv, crohf_g_fcp, g_grad)
      end
      subroutine rohf_fock2grad( nmo, nclosed, nopen,
     $     g_fcv, g_fpv, g_fcp,
     $     g_grad)
C     $Id: rohf_fock2g.F,v 1.3 1996-03-13 22:22:57 d3g681 Exp $
      implicit none
#include "global.fh"
      integer nmo, nclosed, nopen ! [input] Occupation and number of mos
      integer g_fcv, g_fpv, g_fcp ! [input] Fock matrices for each subspace
      integer g_grad            ! [output] Gradient vector
c     
      integer nvir, voff, ooff, oend, xoff, xend
      double precision tol
c     
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c     
      call ga_zero(g_grad)
      xend = nclosed*nvir
      call ga_copy_patch('n', g_fcv, voff, nmo, 1, nclosed,
     $     g_grad, 1, xend, 1, 1 )
      call ga_dscal_patch(g_grad, 1, xend, 1, 1, 4.d0 )
      if (nopen .gt. 0) then
         xoff = xend + 1
         xend = xend + nopen*nvir
         call ga_copy_patch('n', g_fpv, voff, nmo, ooff, oend,
     $        g_grad, xoff, xend, 1, 1 )
         call ga_dscal_patch(g_grad, xoff, xend, 1, 1, 2.d0 )
         xoff = xend + 1
         xend = xend + nopen*nclosed
         call ga_copy_patch('n', g_fcp, ooff, oend, 1, nclosed,
     $        g_grad, xoff, xend, 1, 1 )
         call ga_dscal_patch(g_grad, xoff, xend, 1, 1, 2.d0 )
      endif
c
c     Screen small numbers to ensure hard zeroes are not polluted
c     but don't mess up the quadratic convergence
c
      call ga_maxelt(g_grad, tol)
      tol = max(1d-12,min(1d-6,abs(tol)**2))
**      write(6,*) ' screening grad with ', tol
      call ga_screen(g_grad, tol)
c     
      end
