      subroutine xlm_norm
* $Id$
      implicit double precision (a-h, o-z)
      include 'xlm.h'
c
c     initialize data for unnormalized solid spherical harmonics
c
c     b(m,l) = e(m) / (l + |m|)! , |m| <= l, 0 otherwise
c     c(m,l) = e(m) * (l - |m|)! , |m| <= l, 0 otherwise
c     e(m) = -1 if m > 0 and m odd
c             1 otherwise
c     s(m) = -1 if m < 0
c             0 if m = 0
c             1 if m > 0
c     phase(m) = -1 if m odd  
c              =  1 if m even
c     qn(m,l) = normalization constant required by Condon & Shortley
c               so that the unnormalized xml satisfy
c               qn( m,l) * x( m,l) = r^l * Re (Yml) m >= 0
c               qn(-m,l) * x(-m,l) = r^l * Im (Yml) m >= 0
c               |m| <= l, 0 otherwise
c     w, v are workspace
c
      do l = 0, maxl
         do m = -maxl, maxl
            qn(m,l) = 0.0d0
            b(m,l) = 0.0d0
            c(m,l) = 0.0d0
            w(m,l) = 0.0d0
            v(m,l) = 0.0d0
         end do
      end do
      do m = -maxl, 0
         e(m) = 1.0d0
      end do
      do m = 1, maxl, 2
         e(m) = -1.0d0
         e(m+1) = 1.0d0
         s(m) = 1.0d0
         s(m+1) = 1.0d0
      end do
      do m = -maxl2, maxl2
         s(m) = sign(1.0d0, dble(m))
      end do
      s(0) = 0.0d0
      phase(-maxl4) = 1.0d0
      do m = 1-maxl4, maxl4
         phase(m) = -phase(m-1)
      end do
c
      b(0,0) = 1.0d0
      c(0,0) = 1.0d0
      do l = 1, maxl
         b(0,l) = b(0,l-1)/l
         c(0,l) = c(0,l-1)*l
         do m = 1, l
            b( m, l) = -b(m-1, l) / (l+m)
            b(-m, l) = abs(b(m,l))
            c(m, l) = -c(m-1, l) / (l-m+1)
            c(-m,l) = abs(c(m,l))
         end do
      end do
c
      pi4 = 16.0d0*atan(1.0d0)
      do l = 0, maxl
         qn(0,l) = (2*l+1)/pi4
         do m = 1, l
            qn( m,l) = qn(m-1,l)/((l+m)*(l-m+1))
            qn(-m,l) = qn(m,l)
         end do
      end do
      do l = 0, maxl
         do m = -l, l
            qn(m,l) = e(m) * sqrt(qn(m,l))
         end do
      end do
c
*      write(6,*) ' b coefficients'
*      call xlmprint(maxl, b, maxl)
*      write(6,*) ' c coefficients'
*      call xlmprint(maxl, c, maxl)
*      write(6,*) ' qn coefficients'
*      call xlmprint(maxl, qn, maxl)
c
      end
