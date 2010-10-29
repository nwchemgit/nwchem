      program spherical
* $Id$
      implicit none
      integer memgot
      parameter (memgot=10000)
      double precision dd(memgot)
      integer lmax
      integer memused
c      
      lmax = 8
      memused = (((lmax+1)*(lmax+2))/2)*(2*lmax+1)*(lmax+1)
      write(6,*)' memused = ',memused
      if (memused.gt.memgot) then
        write(6,*)' not enough memory '
        write(6,*)' memgot  = ',memgot
        write(6,*)' memused = ',memused
        stop ' error '
      endif
      call xlmcoeff_test(lmax,dd,.true.)
      end
      subroutine xlmcoeff_test(lmax, d, normalize)
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
c
*:functions
      double precision ddot
      external ddot
*:passed
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
      logical normalize
*:local
      integer e_sz, qn_sz
      parameter (e_sz=2000,qn_sz=5000)
      double precision zzqn(qn_sz),zze(e_sz)
      integer i, j, k, l, m, lp1, ijk
      double precision twolp1

c
c     compute the coefficients of cartesian polynomials for the
c     unnormalized real solid spherical harmonics up to order lmax
c
c::local
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
      call xlm_print_d(d,lmax,'un-normalized')
*rak:      do lp1 = 0,lmax
*rak:        write(luout,*)' d matrix xmlcoeff un normalized'
*rak:        do m = -lp1,lp1
*rak:          do ijk = 1,(((lp1+1)*(lp1+2))/2)
*rak:*             if (d(ijk,m,lp1) .ne. 0.0d0)
*rak:*     $
*rak:            write(luout,1) lp1, m, ijk, d(ijk,m,lp1)
*rak: 1          format(1x,2i3,3x,i3,3x,f20.10)
*rak:          enddo
*rak:        enddo
*rak:      enddo
      if (.not.normalize) then
        write(luout,*)' un-normalized transformation '
        return
      endif

c... allocate space for e
      if (2*lmax+1.gt.e_sz) stop 'e error'
c... allocate space for qn
      if (((2*lmax+1)*(lmax+1)).gt.qn_sz) stop 'qn error'
      call nwints_xlm_norm(zzqn,zze,d,lmax)
      call xlm_print_d(d,lmax,'normalized')
      end
      subroutine xlm_print_d(d,lmax,msg)
      implicit none
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
      character*(*) msg
c
      integer l, m, ijk
c
      do l = 0,lmax
        do m = -l,l
          do ijk = 1,(((l+1)*(l+2))/2)
            write(6,10000)l,m,ijk,d(ijk,m,l),msg
          enddo
        enddo
      enddo
10000 format('<lval=',i2,'> <m=',i3,'>  <ijk=',i5,1pd20.10,1x,a)
      end
*rak:      call xlm_index_test(lmax)
*rak:      subroutine xlm_index_test(lmax)
*rak:      implicit none
*rak:      integer lmax
*rak:      integer lp1,m,ijk,count
*rak:      integer m1, test, lpm1
*rak:c
*rak:      write(6,*)' index check '
*rak:      count = 0
*rak:      do lp1 = 0,lmax
*rak:        do m = -lp1,lp1
*rak:          do ijk = 1,(((lp1+1)*(lp1+2))/2)
*rak:            count = count + 1
*rak:            m1 = m + lp1 + 1
*rak:            lpm1 = lp1-1
*rak:            test = max(0,abs(lpm1))
*rak:            test = test*(lpm1+1)*(lpm1+2)/2
*rak:            test = test*(2*lpm1+1)
*rak:*           test = test + m1 + ijk
*rak:            write(6,*)count,lp1,m,ijk,' m1=',m1,' test=',test
*rak:          enddo
*rak:        enddo
*rak:      enddo
*rak:      end
