C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_grad_b (l_b,n_p,n_o,n_m,n_cb,intp,intm,grad,n_a,
     &    ibug)
*                                                                      *
*   Accumulate gradients with respect to nuclear coordinates from      *
*   integrals with shifted angular momenta on centre B.                *
*                                                                      *
*   l_b (inp) - angular momentum on centre B (before shifting)         *
*   n_p (inp) - cartesian components for l_b+1 = (l_b+2)*(l_b+3)/2     *
*   n_o (inp) - cartesian components for l_b+0 = (l_b+1)*(l_b+2)/2     *
*   n_m (inp) - cartesian components for l_b-1 = (l_b+0)*(l_b+1)/2     *
*   n_cb (inp) - number of contracted functions on centre B            *
*   intp (inp) - integrals for l_b+1 with contraction coefficients     *
*                multiplied by exponents                               *
*   intm (inp) - integrals for l_b-1                                   *
*   grad (out) - gradient of integrals w.r.t. centre B                 *
*   n_a (inp) - total number of functions on centre A                  *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
      integer nall,i,j,k,ix,iy,iz,l_b,n_p,n_o,n_m,n_cb,n_a,ibug
      logical debug_gen,debug_addresses,debug_arrays
      double precision wi,wj,wk,two,
     &    intp(n_p,n_cb*n_a),
     &    intm(n_m,n_cb*n_a),
     &    grad(n_o,n_cb*n_a,3)
      parameter (two = 2.0d00)
*
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
      if (debug_gen) write (LuOut,*) 'Entering ecp_grad_b ...'
      if (debug_addresses) then 
         write (LuOut,*) n_p,n_o,n_m
         write (LuOut,*) n_a,n_cb
      end if
      if (debug_arrays) call ecp_matpr (intp,1,n_p*n_cb,1,n_a,
     &    1,n_p*n_cb,1,n_a,'l+1 integrals','E',78,4)
*
*     Add contributions from derivative of exponent
*
      nall = n_cb*n_a
      ix = 0
      do i = l_b,0,-1
        k = l_b-i
        do j = k,0,-1
          ix = ix+1
          iy = ix+k+1
          iz = iy+1
          call daxpy (nall,two,intp(ix,1),n_p,grad(ix,1,1),n_o)
          call daxpy (nall,two,intp(iy,1),n_p,grad(ix,1,2),n_o)
          call daxpy (nall,two,intp(iz,1),n_p,grad(ix,1,3),n_o)
        end do
      end do
      if (debug_arrays) then
        do i = 1,3
          call ecp_matpr (grad(1,1,i),1,n_o*n_cb,1,n_a,
     &        1,n_o*n_cb,1,n_a,'Gradient integrals','E',78,4)
        end do
      end if
*
*     Add contributions from derivative of monomial
*
      if (n_m .eq. 0) return
      if (debug_arrays) call ecp_matpr (intm,1,n_m*n_cb,1,n_a,
     &    1,n_m*n_cb,1,n_a,'l-1 integrals','E',78,4)
      ix = 0
      do i = l_b-1,0,-1
        k = l_b-i-1
        do j = k,0,-1
          ix = ix+1
          iy = ix+k+1
          iz = iy+1
          wi = -(i+1)
          wj = -(j+1)
          wk = -(k-j+1)
          call daxpy (nall,wi,intm(ix,1),n_m,grad(ix,1,1),n_o)
          call daxpy (nall,wj,intm(ix,1),n_m,grad(iy,1,2),n_o)
          call daxpy (nall,wk,intm(ix,1),n_m,grad(iz,1,3),n_o)
        end do
      end do
      if (debug_arrays) then
        do i = 1,3
          call ecp_matpr (grad(1,1,i),1,n_o*n_cb,1,n_a,
     &        1,n_o*n_cb,1,n_a,'Gradient integrals','E',78,4)
        end do
      end if
*
      return
      end
