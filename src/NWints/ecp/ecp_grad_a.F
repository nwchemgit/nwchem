C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_grad_a (l_a,n_p,n_o,n_m,n_ca,intp,intm,grad,n_b,
     &     ibug)
*                                                                      *
*   Accumulate gradients with respect to nuclear coordinates from      *
*   integrals with shifted angular momenta on centre A.                *
*                                                                      *
*   l_a (inp) - angular momentum on centre A (before shifting)         *
*   n_p (inp) - cartesian components for l_a+1 = (l_a+2)*(l_a+3)/2     *
*   n_o (inp) - cartesian components for l_a+0 = (l_a+1)*(l_a+2)/2     *
*   n_m (inp) - cartesian components for l_a-1 = (l_a+0)*(l_a+1)/2     *
*   n_ca (inp) - number of contracted functions on centre A            *
*   intp (inp) - integrals for l_a+1 with contraction coefficients     *
*                multiplied by exponents                               *
*   intm (inp) - integrals for l_a-1                                   *
*   grad (out) - gradient of integrals w.r.t. centre A                 *
*   n_b (inp) - total number of functions on centre B                  *
*   ibug - debug flag. 0 for no debug, 1 for address printing, 2 for   *
*          array printing, 3 for both.                                 *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "stdio.fh"
      integer i,j,k,m,ix,iy,iz,l_a,n_p,n_o,n_m,n_ca,n_b,ibug
      logical debug_gen,debug_addresses,debug_arrays
      double precision wi,wj,wk,two,
     &    intp(n_b,n_p,n_ca),
     &    intm(n_b,n_m,n_ca),
     &    grad(n_b,n_o,n_ca,3)
*
      parameter (two = 2.0d00)
      debug_gen = ibug .gt. 0
      debug_addresses = mod(ibug,2) .eq. 1
      debug_arrays = mod(ibug,10)/2 .eq. 1
      if (debug_gen) write (LuOut,*) 'Entering ecp_grad_a ...'
      if (debug_addresses) then 
         write (LuOut,*) n_p,n_o,n_m
         write (LuOut,*) n_b,n_ca
      end if
      if (debug_arrays) call ecp_matpr (intp,1,n_b,1,n_p*n_ca,
     &    1,n_b,1,n_p*n_ca,'l+1 integrals','E',78,4)
*
*     Add contributions from derivative of exponent
*
      ix = 0
      do i = l_a,0,-1
        k = l_a-i
        do j = k,0,-1
          ix = ix+1
          iy = ix+k+1
          iz = iy+1
          do m = 1,n_ca
            call daxpy (n_b,two,intp(1,ix,m),1,grad(1,ix,m,1),1)
            call daxpy (n_b,two,intp(1,iy,m),1,grad(1,ix,m,2),1)
            call daxpy (n_b,two,intp(1,iz,m),1,grad(1,ix,m,3),1)
          end do
        end do
      end do
      if (debug_arrays) then
        do i = 1,3
          call ecp_matpr (grad(1,1,1,i),1,n_b,1,n_o*n_ca,
     &        1,n_b,1,n_o*n_ca,'Gradient integrals','E',78,4)
        end do
      end if
*
*     Add contributions from derivative of monomial
*
      if (n_m .eq. 0) return
      if (debug_arrays) call ecp_matpr (intm,1,n_b,1,n_m*n_ca,
     &    1,n_b,1,n_m*n_ca,'l-1 integrals','E',78,4)
      ix = 0
      do i = l_a-1,0,-1
        k = l_a-i-1
        do j = k,0,-1
          ix = ix+1
          iy = ix+k+1
          iz = iy+1
          wi = -(i+1)
          wj = -(j+1)
          wk = -(k-j+1)
          if (debug_addresses) then
            write (LuOut,*) i,j,k-j
            write (LuOut,*) ix,iy,iz
            write (LuOut,*) wi,wj,wk
          end if
          do m = 1,n_ca
            call daxpy (n_b,wi,intm(1,ix,m),1,grad(1,ix,m,1),1)
            call daxpy (n_b,wj,intm(1,ix,m),1,grad(1,iy,m,2),1)
            call daxpy (n_b,wk,intm(1,ix,m),1,grad(1,iz,m,3),1)
          end do
        end do
      end do
      if (debug_arrays) then
        do i = 1,3
          call ecp_matpr (grad(1,1,1,i),1,n_b,1,n_o*n_ca,
     &        1,n_b,1,n_o*n_ca,'Gradient integrals','E',78,4)
        end do
      end if
*
      return
      end
