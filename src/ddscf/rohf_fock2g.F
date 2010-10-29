      subroutine rohf_get_grad(g_grad)
      implicit none
#include "cscf.fh"
#include "crohf.fh"
#include "mafdecls.fh"
      integer g_grad
c
      integer nvir
      double precision maxelt
c
      call rohf_fock2grad(nmo, nclosed, nopen, crohf_g_fcv, 
     $     crohf_g_fpv, crohf_g_fcp, g_grad)
c
      nvir = nmo - nclosed - nopen
      if(oadapt) call scf_sym_screen(nclosed, nopen, nvir,int_mb(k_irs),
     $  g_grad, 0, .true., maxelt)
c
      end
      subroutine rohf_fock2grad( nmo, nclosed, nopen,
     $     g_fcv, g_fpv, g_fcp,
     $     g_grad)
C     $Id$
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
      call ga_mat_to_vec(g_fcv, voff, nmo, 1, nclosed, g_grad, 1, 1,
     $    0.0d0, '=')
*      call ga_copy_patch('n', g_fcv, voff, nmo, 1, nclosed,
*     $     g_grad, 1, xend, 1, 1 )
      call ga_dscal_patch(g_grad, 1, xend, 1, 1, 4.d0 )
      if (nopen .gt. 0) then
         xoff = xend + 1
         xend = xend + nopen*nvir
         call ga_mat_to_vec(g_fpv, voff, nmo, ooff, oend,
     $      g_grad, xoff, 1, 0.0d0, '=')
*         call ga_copy_patch('n', g_fpv, voff, nmo, ooff, oend,
*     $        g_grad, xoff, xend, 1, 1 )
         call ga_dscal_patch(g_grad, xoff, xend, 1, 1, 2.d0 )
         xoff = xend + 1
         xend = xend + nopen*nclosed
         call ga_mat_to_vec(g_fcp, ooff, oend, 1, nclosed,
     $      g_grad, xoff, 1, 0.0d0, '=')
*         call ga_copy_patch('n', g_fcp, ooff, oend, 1, nclosed,
*     $        g_grad, xoff, xend, 1, 1 )
         call ga_dscal_patch(g_grad, xoff, xend, 1, 1, 2.d0 )
      endif
c
c     Screen small numbers to ensure hard zeroes are not polluted
c     but don't mess up the quadratic convergence
c
c     Cannot screen with a large number (greater than accuracy in integrals) 
c     since in a non-abelian group we could break symmetry ... either
c     screen with a small threshold or screen with a large one and then
c     enforce symmetry (which is expensive in high order groups).
c     Also, threshold must tend to zero for high convergence.
c
      call ga_maxelt(g_grad, tol)
      tol = min(tol*tol,1d-12)
      call ga_screen(g_grad, tol)
c     
      end
