      subroutine xlmcoeff(lmax, d, normalize)
* $Id$
* 
* Note that the routine assumes that d is dimensioned in a
* very specific fashion that actually consumes twice as
* much memory as necessary ... however, simiplicity is
* probably the best thing at first.
* 
* d(i, m, l) = coefficient of the ith cartesian harmonic of
*              rank l in the expansion of the real solid 
*              spherical harmonic xlm (m = -l,...,l)
* 
* i.e.,
*        x(m,l) = sum(i) cartesian(i) * d(i, m, l)
*
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
c
*:functions
      double precision ddot
      external ddot
*:passed
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
      logical normalize
*:local
      integer i, j, k, l, m, lp1, ijk
      double precision twolp1

#if defined (REALLY_BAD_NORMALIZATION)
      double precision norm
#endif
c
c     compute the coefficients of cartesian polynomials for the
c     unnormalized real solid spherical harmonics up to order lmax
c
c::local
      integer h_xlm_e, k_xlm_e, h_xlm_qn, k_xlm_qn
*statement functions:
      integer itri
      integer ind
      itri(i,j)  = (i*(i-1))/2 + j
      ind(i,j,l) = itri(l-i+1,l-i-j+1) ! Index of x^i*y^j*z^(l-i-j)
c
      call dcopy((((lmax+1)*(lmax+2))/2)*(2*lmax+1)*(lmax+1),

     $     0.0d0, 0, d, 1)
c
      d(1, 0,0) = 1.0d0         ! Constant
c
      if (lmax .ge. 1) then
         d(1, 1,1) = 1.0d0      ! x
         d(2,-1,1) = 1.0d0      ! y
         d(3, 0,1) = 1.0d0      ! z
      end if
c
      do l = 1, lmax-1
         lp1 = l + 1
         twolp1 = 2*l + 1
c     
         do i = lp1, 0, -1
            do j = lp1-i, 0, -1
               k = lp1 - i - j
c
               ijk = ind(i,j,lp1)
c
*     q( lp1, lp1) = twolp1*(x*q( l,l) - y*q(-l,l))
               if (i.gt.0) d(ijk,lp1,lp1) =
     $              d(ijk,lp1,lp1) + twolp1*
     $              d(ind(i-1,j,l),l,l)
               if (j.gt.0) d(ijk,lp1,lp1) = 
     $              d(ijk,lp1,lp1) - twolp1*
     $              d(ind(i,j-1,l),-l,l)

*     q(-lp1, lp1) = twolp1*(x*q(-l,l) + y*q( l,l))
               if (i.gt.0) d(ijk,-lp1,lp1) =
     $              d(ijk,-lp1,lp1) + twolp1*
     $              d(ind(i-1,j,l),-l,l)
               if (j.gt.0) d(ijk,-lp1,lp1) = 
     $              d(ijk,-lp1,lp1) + twolp1*
     $              d(ind(i,j-1,l),l,l)

*     q(0,lp1) = (twolp1*z*q(0,l) - l*r2*q(0,l-1))/(l+1)
               if (k.gt.0) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) + twolp1*
     $              d(ind(i,j,l),0,l) / (l+1)
               if (i.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i-2,j,l-1),0,l-1) / (l+1)
               if (j.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i,j-2,l-1),0,l-1) / (l+1)
               if (k.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i,j,l-1),0,l-1) / (l+1)

               do m = 1,l
*     q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1))/(l-m+1)
*     q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1))/(l-m+1)

                  if (k.gt.0) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) +
     $                    twolp1*d(ind(i,j,l), m,l)/(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) +
     $                    twolp1*d(ind(i,j,l),-m,l)/(l-m+1)
                  end if
                  if (i.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i-2,j,l-1), m,l-1)/(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i-2,j,l-1),-m,l-1)/(l-m+1)
                  end if
                  if (j.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i,j-2,l-1), m,l-1)/(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i,j-2,l-1),-m,l-1)/(l-m+1)
                  end if
                  if (k.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i,j,l-1), m,l-1)/(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i,j,l-1),-m,l-1)/(l-m+1)
                  end if
               end do
c
            end do
         end do
      end do
c
#if defined(DEBUG)
      do lp1 = 0,lmax
        write(luout,*)' d matrix xmlcoeff'
        do m = -lp1,lp1
          do ijk = 1,(((lp1+1)*(lp1+2))/2)
*             if (d(ijk,m,lp1) .ne. 0.0d0)
*     $
            write(luout,1) lp1, m, ijk, d(ijk,m,lp1)
 1          format(1x,2i3,3x,i3,3x,f20.10)
          enddo
        enddo
      enddo
#endif
#if defined (REALLY_BAD_NORMALIZATION)
* normalize transformation coefficients (start with d functions)
*
      if (.not.normalize) return
      do lp1 = 2,lmax
        do m = -lp1,lp1
          norm = ddot((((lp1+1)*(lp1+2))/2),d(1,m,lp1),1,d(1,m,lp1),1)
          norm = 1.0d00/sqrt(norm)
          call dscal((((lp1+1)*(lp1+2))/2),norm,d(1,m,lp1),1)
        enddo
      enddo
*
#endif
      if (.not.normalize) then
        write(luout,*)' un-normalized transformation '
        return
      endif

c... allocate space for e
      if (.not.ma_push_get(mt_dbl,(2*lmax+1),
     &      ' tmp e buf xlm ',
     &      h_xlm_e, k_xlm_e)) call errquit
     &      ('xlmcoeff: error allocating e bug ',911, MEM_ERR)
* does not need to be zeroed      
c... allocate space for qn
      if (.not.ma_push_get(mt_dbl,((2*lmax+1)*(lmax+1)),
     &      ' tmp qn buf xlm ',
     &      h_xlm_qn, k_xlm_qn)) call errquit
     &      ('xlmcoeff: error allocating qn bug ',911, MEM_ERR)
* does not need to be zeroed      

      call nwints_xlm_norm(dbl_mb(k_xlm_qn),dbl_mb(k_xlm_e),d,lmax)
c... free tmp space
      if (.not.MA_pop_stack(h_xlm_qn)) call errquit
     &      ('xlmcoeff: error poping qn bug ',911, MEM_ERR)
      if (.not.MA_pop_stack(h_xlm_e)) call errquit
     &      ('xlmcoeff: error poping e bug ',911, MEM_ERR)
c
#if defined(DEBUG)
      write(luout,*)' normalized coefficients '
      do lp1 = 0,lmax
        do m = -lp1,lp1
          do ijk = 1,(((lp1+1)*(lp1+2))/2)
            write(luout,1) lp1, m, ijk, d(ijk,m,lp1)
          enddo
        enddo
      enddo
#endif
      end
