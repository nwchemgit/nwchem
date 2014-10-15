      subroutine xlm_norm(lmax, q, lq)
*
* $Id$
*
      implicit none
#include "xlm.fh"
      integer lmax, lq
      double precision q(-lq:lq,0:lmax)
c     
c     Return in q the normalization constants qn ... see xlm_init
c
      integer l, m
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = qn(m,l)
         enddo
      enddo
c
      end
      subroutine xlm_thresh(lmax, q, lq, thresh, dist, luse)
      implicit none
#include "xlm.fh"
      integer lmax, lq, luse
      double precision q(-lq:lq,0:lmax), thresh, dist
c
c     Return in luse the angular momentum required to obtain
c     the required accuracy in the potential at a distance dist.
c
c     For use in the FMM dist should be (ws+0.5)*box.  
c
c     If all components are negligble luse=-1 is returned.
c
      integer l, m
      double precision scale, test, sum
c
      scale = 1.0d0/dist**(lmax+1)
      sum = 0.0d0
      do l = lmax,0,-1
         test = 0.0d0
         do m = -l, l
            test = test + abs(q(m,l)*fieldmax(m,l))
         enddo
         sum = sum + test*scale
         if (sum .gt. thresh) goto 10
         scale = scale * dist
      enddo
 10   luse = l
c
      end
      subroutine xlm_m2l_thresh(lmax, q, lq, thresh, dist, box, luse)
      implicit none
#include "xlm.fh"
      integer lmax, lq, luse
      double precision q(-lq:lq,0:lmax), thresh, dist, box
c
c     Return in luse the angular momentum required to obtain
c     the required accuracy in the Taylor series of the potential 
c     at a distance dist in a box of size box.
c
c     For use in the FMM dist should be (ws+0.5)*box.  
c
c     If all components are negligble luse=-1 is returned.
c
      integer l, m, k
      double precision scale, test, sum, factorial
      external factorial
c
      scale = 1.0d0/dist**(lmax+1)
      sum = 0.0d0
      do l = lmax,0,-1
         test = 0.0d0
         do m = -l, l
            test = test + abs(q(m,l)*fieldmax(m,l))
         enddo
         sum = sum + test*scale
         if (sum .gt. thresh) goto 10
         scale = scale * dist
      enddo
 10   continue
c
c     l is the order of mpole required to get the potential
c     accurate ... but perhaps the taylor series needs even
c     higher order?
c
      if (l.gt.-1 .and. l.lt.lmax) then
         k = l+1
         test = test*(0.87d0*box/dist)**k * factorial(l+k)/
     $        (factorial(l)*factorial(k))
         do k = l+1, lmax
            if (test .lt. thresh) goto 20
            test = test*(l+k+1)*0.5d0*box/(dist*(k+1))
         end do
 20      k = k - 1
         if (l .lt. k) then
            write(6,*) ' Increasing l ', l, k
            l = k
         endif
      endif
c
      luse = l
c
      end
      subroutine xlm_slower(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      integer lmax, lq
      double precision x, y, z, q(-lq:lq, 0:lmax)
c
      integer l, m, lp1
      double precision r2, twolp1
c
c     compute the unnormalized real solid spherical harmonics
c     up to order lmax at the point x, y, z.
c     lq specifies the dimension of q.
c
      if (lmax .gt. lq) then
         write(6,*) ' invalid dimension for q ', lmax, lq
         stop 1
      end if
***      call dfill((2*lq+1)*(lmax+1), 0.0d0, q, 1)
      do l = 0,lmax             ! NOTE that |m| = l+1 is referenced
         do m = -lq,lq
            q(m,l) = 0.0d0
         end do
      end do
      q(0,0) = 1.0d0            ! MUST be after zero
      if (lmax.ge.1) then
         q(-1,1) = y
         q( 0,1) = z
         q( 1,1) = x
      end if
c
      r2 = x*x + y*y + z*z
      do l = 1, lmax-1
         lp1 = l + 1
         twolp1 = 2*l + 1
         q( lp1, lp1) = twolp1*(x*q( l,l) - y*q(-l,l))
         q(-lp1, lp1) = twolp1*(x*q(-l,l) + y*q( l,l))
         q(0,lp1) = (twolp1*z*q(0,l) - l*r2*q(0,l-1))*fastdiv(l+1)
c$$$         do m = 1,l
c$$$            q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1))/(l-m+1)
c$$$            q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1))/(l-m+1)
c$$$         end do
         do m = 1,l
            q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1)) *
     $           fastdiv(l-m+1)
            q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1)) *
     $           fastdiv(l-m+1)
         end do
      end do
c
      end
      subroutine xlm(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      integer lmax, lq
      double precision x, y, z, q(-lq:lq, 0:lmax)
c
      integer l, m, lp1
      double precision r2, twolp1, twolp1_z, lmr2
c
c     compute the unnormalized real solid spherical harmonics
c     up to order lmax at the point x, y, z.
c     lq specifies the dimension of q.
c
c     fast(er) version
c
      q(0,0) = 1.0d0
      if (lmax .eq. 0) return
      q(-1,1) = y
      q( 0,1) = z
      q( 1,1) = x
      if (lmax .eq. 1) return
      q(-2,2) = 6.0d0*x*y
      q(-1,2) = 3.0d0*y*z
      q( 0,2) = z*z - 0.5d0*(x*x+y*y)
      q( 1,2) = 3.0d0*x*z
      q( 2,2) = 3.0d0*(x*x-y*y)
      if (lmax .eq. 2) return
c
      r2 = x*x + y*y + z*z
      twolp1 = 5.0d0
      do l = 2, lmax-1
         lp1 = l + 1
         q( lp1, lp1) = twolp1*(x*q( l,l) - y*q(-l,l))
         q(-lp1, lp1) = twolp1*(x*q(-l,l) + y*q( l,l))
c
         twolp1_z = twolp1 * z
c
         lmr2 = dble(l)*r2
         q(0,lp1) = (twolp1_z*q(0,l) - lmr2*q(0,l-1))*fastdiv(l+1)
         do m = 1,l-1
            lmr2 = lmr2 + r2
            q( m,lp1) = (twolp1_z*q( m,l) - lmr2*q( m,l-1)) *
     $           fastdiv(l-m+1)
            q(-m,lp1) = (twolp1_z*q(-m,l) - lmr2*q(-m,l-1)) *
     $           fastdiv(l-m+1)
         end do
         q( l,lp1) = twolp1_z*q( l,l)
         q(-l,lp1) = twolp1_z*q(-l,l)
c
         twolp1 = twolp1 + 2.0d0
      end do
c
      end
      double precision function factorial(n)
      implicit none
      integer n, i
      double precision f
c
      f = 1.0d0
      do i = 1, n
         f = f * dble(i)
      enddo
c
      factorial = f
c
      end
      double precision function double_factorial(n)
      implicit none
      integer n, i
      double precision f
c
      f = 1.0d0
      do i = n, 1, -2
         f = f * dble(i)
      enddo
c
      double_factorial = f
c
      end
      subroutine xlm_coeff_inv(lmax,d,dinv)
      implicit none
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
      double precision 
     $     dinv(((lmax+1)*(lmax+2)*(lmax+3))/6,-lmax:lmax,0:lmax)
c
c     NOTE differing dimensions for d and dinv
c
c     Compute the coefficients of the solid harmonics for expanding
c     the cartesian basis ... the inverse of xlm_coeff.  I.e., return
c     in dinv
c
c     x^i y^j z^k = sum(l=i+j+k,0,-2) sum(m=-l..l) 
c     .                     r^(l-i-j-k) dinv(ijk,m,l) Xlm
c
c     Also returned (as an extra, added bonus for the wary shopper)
c     is d from xlm_coeff.
c
c     Xlm = sum(ijk) x^i y^j z^k d(ijk,m,l)
c     
c     dinv(ijk,lm) = <Xlm * x^iy^jz^k>/<Xlm^2>  
c
c     and use expansion of Xlm to obtain
c
c     dinv(ijk,lm) = sum(i'j'k') di'j'k',lm <x^(i+i')y^(j+j')z^(k+k')> 
c     .              / <Xlm^2> 
c
c     <x^i y^j z^k> = 4Pi (i-1)!!(j-1)!!(k-1)!!/(i+j+k+1)!!
c
c     <Xlm^2> = 4Pi [(l+|m|)!/(l-|m|)!] / [(2-delta(m,0))(2l+1)]
c
c     (<> = integral over unit sphere)
c
      integer i, j, k, l, lp, m, mp, itri, ind, ijk, loff,
     $     ip, jp, kp, ijkp
      double precision numerator, denominator, sum
      double precision double_factorial, factorial
      external double_factorial, factorial
c
      itri(i,j)  = (i*(i-1))/2 + j
      ind(i,j,l) = itri(l-i+1,l-i-j+1) ! Index of x^i*y^j*z^(l-i-j)
      loff(l) = ((l*(l+1)*(l+2))/6)
c
      call xlm_coeff(lmax, d)
c
      call dfill((((lmax+1)*(lmax+2)*(lmax+3))/6)*(2*lmax+1)*(lmax+1),
     $     0.0d0, dinv, 1)
c
      do l = 0, lmax
         do i = l, 0, -1
            do j = l-i, 0, -1
               k = l-i-j
               ijk = loff(l) + ind(i,j,l)
               do lp = l, 0, -2
                  do ip = lp, 0, -1
                     do jp = lp-ip, 0, -1
                        kp = lp-ip-jp
                        if ( (mod(i+ip,2).eq.0) .and. 
     $                       (mod(j+jp,2).eq.0) .and. 
     $                       (mod(k+kp,2).eq.0) ) then
                        numerator = 
     $                       double_factorial(i+ip-1)*
     $                       double_factorial(j+jp-1)*
     $                       double_factorial(k+kp-1)/
     $                       double_factorial(i+ip+j+jp+k+kp+1)
                        ijkp = ind(ip,jp,lp)
                        do mp = 0, lp
                           denominator = factorial(lp+mp)/
     $                          (factorial(lp-mp)*(lp+lp+1))
*                           write(6,*) l, i, j, k
*                           write(6,*) lp, ip, jp, kp
*                           write(6,*) ijk, ijkp
*                           write(6,*) numerator, denominator, 
*     $                          d(ijkp, mp,lp),
*     $                          d(ijkp, -mp,lp)
                           if (mp .ne. 0) 
     $                          denominator = denominator*0.5d0
                           dinv(ijk,mp,lp) = dinv(ijk,mp,lp) + 
     $                          d(ijkp, mp,lp)*numerator/denominator
                           if (mp .ne. 0) 
     $                          dinv(ijk,-mp,lp) = dinv(ijk,-mp,lp) + 
     $                          d(ijkp,-mp,lp)*numerator/denominator
                        enddo
                        endif
                     enddo
                  enddo
               enddo
*               write(6,22) i, j, k, ijk
* 22            format(' Expansion for ', 3i5, 2x, i5)
*               do lp = l, 0, -2
*                  do mp = -lp,lp
*                     if (dinv(ijk,mp,lp) .ne. 0.0d0)
*     $                    write(6,1) mp, lp, l-lp, dinv(ijk,mp,lp)
* 1                   format(1x,2i5,3x,i5,3x,g20.10)
*                  end do
*               enddo
            enddo
         enddo
      enddo
c
c     Check what we can for paranoia.  The norm of a solid harmonic
c     is unchanged by converting to cartesians and back again.
c
      do l = 0, lmax
         do m = -l, l
            sum = 0.0d0
            do i = l, 0, -1
               do j = l-i, 0, -1
                  k = l-i-j
                  ijk = ind(i,j,l)
                  sum = sum + d(ijk,m,l)*dinv(loff(l)+ijk,m,l)
               enddo
            enddo
            if (abs(sum-1.0d0) .gt. 1d-10) then
               write(0,*) ' errror ',abs(sum-1.0d0)
	    endif
         enddo
      enddo
c
      end
      
      subroutine xlm_coeff(lmax, d)
      implicit none
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
c
c     compute the coefficients of cartesian polynomials for the
c     unnormalized real solid spherical harmonics up to order lmax
c
      integer i, j, k, l, m, itri, ind, lp1, ijk
      double precision twolp1

      itri(i,j)  = (i*(i-1))/2 + j
      ind(i,j,l) = itri(l-i+1,l-i-j+1) ! Index of x^i*y^j*z^(l-i-j)
c
      call dfill((((lmax+1)*(lmax+2))/2)*(2*lmax+1)*(lmax+1),
     $     0.0d0, d, 1)
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
     $              d(ind(i,j,l),0,l) / dble(l+1)
               if (i.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i-2,j,l-1),0,l-1) / dble(l+1)
               if (j.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i,j-2,l-1),0,l-1) / dble(l+1)
               if (k.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i,j,l-1),0,l-1) / dble(l+1)

               do m = 1,l
*     q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1))/(l-m+1)
*     q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1))/(l-m+1)

                  if (k.gt.0) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) +
     $                    twolp1*d(ind(i,j,l), m,l)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) +
     $                    twolp1*d(ind(i,j,l),-m,l)/dble(l-m+1)
                  end if
                  if (i.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i-2,j,l-1), m,l-1)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i-2,j,l-1),-m,l-1)/dble(l-m+1)
                  end if
                  if (j.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i,j-2,l-1), m,l-1)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i,j-2,l-1),-m,l-1)/dble(l-m+1)
                  end if
                  if (k.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i,j,l-1), m,l-1)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i,j,l-1),-m,l-1)/dble(l-m+1)
                  end if
               end do
c
c$$$               do m = -lp1,lp1
c$$$                  if (d(ijk,m,lp1) .ne. 0.0d0)
c$$$     $                 write(6,1) lp1, m, i, j, k, d(ijk,m,lp1)
c$$$ 1                format(1x,2i5,3x,3i3,3x,g20.10)
c$$$               end do
            end do
         end do
      end do
c$$$      stop 31
c
      end
      subroutine xlm_init
      implicit none
#include "xlm.fh"
c
c     initialize data for unnormalized solid spherical harmonics
c
c     b(m,l) = e(m) / (l + |m|)! , |m| <= l, 0 otherwise
c     c(m,l) = e(m) * (l - |m|)! , |m| <= l, 0 otherwise
c     bcp0(m,l) = (b(m,l)*c(-m,l) + b(-m,l)*c(m,l))*phase(l+m) , m > 0
c     .           b(0,l)*c(0,l)*phase(l), m = 0
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
c     rqn = 1.0d0 / qn
c     ee(m2,m1) = used to accellerate multipole_to_local
c     w, v are workspace
c     fastdiv(l) = 1.0d0/l, l=1,max2l+1
c
      integer l, m, nsample
      double precision pi4, pi
      integer m1, m2, m1m2, m1mm2, itheta, iphi
      double precision pe1, pe1m, theta, phi, x, y, z, test
      logical initialized
      data initialized/.false./

      if (initialized) return
      initialized = .true.

      numphi = 0                ! For caching rotation info
      numtheta = 0
      call ifill(maxhash, -1, iphis, 1)
      call ifill(maxhash, -1, ithetas, 1)

      do l = 1, maxl2+1
         fastdiv(l) = 1.0d0/dble(l)
      end do
      do l = 0, maxl
         do m = -maxl, maxl
            qn(m,l) = 0.0d0
            w(m,l) = 0.0d0
            v(m,l) = 0.0d0
         end do
      end do
      do l = 0, maxl2
         do m = -maxl2, maxl2
            b(m,l) = 0.0d0
            c(m,l) = 0.0d0
         end do
      end do
      do m = -maxl2, 0
         e(m) = 1.0d0
      end do
      do m = 1, maxl2, 2
         e(m) = -1.0d0
         if (m.lt.maxl2) e(m+1) = 1.0d0
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
      do l = 1, maxl2
         b(0,l) = b(0,l-1)/dble(l)
         c(0,l) = c(0,l-1)*dble(l)
         do m = 1, l
            b( m, l) = -b(m-1, l) / dble(l+m)
            b(-m, l) = abs(b(m,l))
            c(m, l) = -c(m-1, l) / dble(l-m+1)
            c(-m,l) = abs(c(m,l))
         end do
      end do
c
      do l = 0, maxl
         bcp0(0,l) = b(0,l)*c(0,l)*phase(l)
         do m = 1, l
            bcp0(m,l) = (b(m,l)*c(-m,l) + b(-m,l)*c(m,l))*phase(l+m)
         end do
      end do
c
c     This defintion of the normalization constant is now correct
c     and agrees with numerical quadrature of the functions.
c     It also yields unitary rotation matrices.  
c
c     ANY changes must be checked for impact everywhere qn is used
c     since the correct normalization and unitary property 
c     are relied upon.
c
      pi4 = 4.0d0*3.1415926535897932d0
      do l = 0, maxl
         do m = -l, l
            qn( m,l) = sqrt(2.0d0*(2*l+1)*b(m,l)*c(m,l)/pi4)*e(m)
            rqn(m,l) = 1.0d0/qn(m,l)
         end do
         qn(0,l) = qn(0,l) / sqrt(2.0d0)
         rqn(0,l) = 1.0d0/qn(0,l)
      end do
c
      do l = 0, maxl
         do m = -l,l
            twobc(m,l) = 2.0d0*b(m,l)*c(m,l)*rqn(m,l)
            rtwobcqnbp(m,l) = b(m,l)*phase(l+m)/twobc(m,l)
         enddo
         rtwobcqnbp(0,l) = rtwobcqnbp(0,l) * 2.0d0
      enddo
c
      do l = 0,maxl
         brqnp(0,l) = b(0,l)*0.5d0*rqn(0,l)
         do m = 1, l
            brqnp( m,l) = b(m,l)*rqn( m,l)
            brqnp(-m,l) = b(m,l)*rqn(-m,l)*phase(m)
         enddo
      enddo
c
      do m1 = 0, maxl
         pe1 = phase(m1)*e(m1)
         pe1m= phase(-m1)
*         if (m1 .eq. 0) pe1 = pe1 * 0.5d0
*         if (m1 .eq. 0) pe1m = pe1m * 0.5d0
         do m2 = 0, maxl
            m1m2 = m1+m2
            m1mm2 = m1-m2
            ee( m2, m1) = (e(m2)*pe1*e(m1m2)+pe1m*e(-m1m2))*e(m1)
            ee(-m2, m1) = (pe1*e(m1mm2)+e(m2)*pe1m*e(-m1mm2))*e(m1)
         end do
      end do
c
c     Compute the maximum value taken on by the field terms over
c     the surface of a unit sphere - used for screening interactions
c
      do l = 0, maxl
         do m = -l, l
            fieldmax(m,l) = 0.0d0
         enddo
      enddo
      nsample = maxl*2          ! No. of sampling points in pi/2
      pi = pi4/4.0d0
      do itheta = 0, 2*nsample  ! Search in +ve octant
         do iphi = 0, nsample
            phi = dble(iphi)*0.5d0*pi/dble(nsample)
            theta = dble(itheta)*0.5d0*pi/dble(2*nsample)

            z = cos(theta)
            x = sin(theta)*cos(phi)
            y = sin(theta)*sin(phi)

            call xlm(maxl, x, y, z, w, maxl)
            do l = 0, maxl
               do m = -l, l
                  test = pi4*w(m,l)*qn(m,l)*qn(m,l)/dble(l+l+1)
                  fieldmax(m,l) = max(abs(test),fieldmax(m,l))
               enddo
            enddo
         enddo
      enddo
c
*      write(6,*) ' b coefficients'
*      call xlm_print(maxl, b, maxl)
*      write(6,*) ' c coefficients'
*      call xlm_print(maxl, c, maxl)
*      write(6,*) ' qn coefficients'
*      call xlm_print(maxl, qn, maxl)
c
      end
      double precision function xlm_combin(n,m)
      implicit none
      integer n, m
c
c     Return n!/((n-m)!m!)
c
      double precision result
      integer i
c
      result = 1.0d0
      do i = n-m+1, n
         result = result * dble(i)
      enddo
      do i = 1, m
         result = result / dble(i)
      enddo
c
      xlm_combin = result
c
      end
      subroutine xlm_apply_z_y_rotation(lmax,alpha,beta,q,lq)
      implicit none
#include "xlm.fh"
      integer lmax,lq
      double precision alpha, beta, q(-lq:lq,0:lmax)
c
c     Apply a rotation of alpha radians about Z-axis followed
c     by a rotation of beta radians about the new Y-axis.
c
      call xlm_apply_z_rotation(lmax,alpha,q,lq)
      call xlm_apply_y_rotation(lmax,beta,q,lq)
c
      end
      subroutine xlm_apply_z_rotation(lmax,alpha,q,lq)
      implicit none
#include "xlm.fh"
      integer lmax,lq
      double precision alpha, q(-lq:lq,0:lmax)
c
c     Apply a rotation about the z axis of alpha radians
c     to the coefficients of unnormalized solid harmonics
c
      double precision cosa(0:maxl), sina(0:maxl), xp, xm
      integer l, m
c
      do m = 0,lmax
         cosa(m) = cos(m*alpha)
         sina(m) = sin(m*alpha)
      enddo
c
      do l = 1, lmax
         do m = 1, l
            xp = cosa(m)*q(m,l) + sina(m)*q(-m,l)
            xm =-sina(m)*q(m,l) + cosa(m)*q(-m,l)
            q( m,l) = xp
            q(-m,l) = xm
         enddo
      enddo
c
      end
      subroutine xlm_apply_y_rotation(lmax,beta,q,lq)
      implicit none
#include "xlm.fh"
      integer lmax,lq
      double precision beta, q(-lq:lq,0:lmax)
c
c     Apply a rotation about the y axis of beta radians
c     to the coefficients of unnormalized solid harmonics
c
      double precision dp((maxl+1)*(maxl+1)), 
     $     dm(maxl*maxl), qq(-maxl:maxl)
      integer l, m, mp, ind
c
      do l = 0, lmax
         call xlm_y_rotation_matrix(l,beta,dp,dm)
         do m = -l, l
            qq(m) = q(m,l)*rqn(m,l) ! Precondition
            q(m,l) = 0.0d0
         enddo
         ind = 1
         do mp = 0, l
            do m = 0, l
               q(m,l) = q(m,l) + dp(ind)*qq(mp)
               ind = ind + 1
            enddo
         enddo
         if (l .gt. 0) then
            ind = 1
            do mp = -l,-1
               do m = -l,-1
                  q(m,l) = q(m,l) + dm(ind)*qq(mp)
                  ind = ind + 1
               enddo
            enddo
         endif
         do m = -l, l
            q(m,l) = q(m,l)*qn(m,l) ! Undo preconditioning
         enddo
      enddo
c
      end
      subroutine xlm_y_rotation_matrix(l,beta,dp,dm)
      implicit none
#include "errquit.fh"
#include "xlm.fh"
c
      integer l
      double precision beta, dp(0:l,0:l),dm(l,l)
c
c     Return in dp and dm the matrices corresponding to the
c     operator rotating about the y axis by angle beta.
c
c     For numerical stability, the rotation matrices are
c     PRECONDITIONED with qn(m,l), the normalization constant
c     for the funciton.
c
c     Note that the preconditioned matrices are unitary.
c
c     Thus,
c
c     qn(m,l)*Xlm (theta+angle) = sum(m') Xlm' (theta)*qn(m',l)*dp(m',m)
c     qn(-m,l)*Xl-m(theta+angle) = sum(m') Xl-m'(theta)*qn(-m',l)*dp(-m',-m)
c
c     Or, to apply the operator to a set of coefficients consider
c
c     f(theta) = sum(m) Cm Xlm(theta) 
c     .        = sum(m) (Cm/qn(m))*(Xlm(theta)*qn(m))
c     .        = sum(m,m') (Cm/qn(m))*(Xlm'(theta+angle)*qn(m')*Dm'm)
c     .        = sum(m') Xlm'(theta+angle)*qn(m')*[sum(m) Dm'm*(Cm/qn(m))]
c     
c     i.e., the coeff of Xlm'(theta+beta) is qn(m')*sum(m) Dm'm*(Cm/qn(m))
c     which are formed by first dividing the coeffs by qn(m), multiplying 
c     by Dm'm and then multiplying by qn(m').
c
c     This is much more numerically stable than not scaling.
c
c     This routine is meant to be correct and accurate, but there
c     are MUCH faster ways of generating these matrices.
c     
      integer k, info, ipiv(maxl+1), m
      double precision theta, pi, x, y, z, rx, ry, rz, phi
      double precision a(0:maxl,0:maxl)
c
      pi = 3.1415926535897932d0
c
      do k = 0, l
         theta = dble(k)*pi/dble(l+1)
         z = cos(theta+beta)
         x = sin(theta+beta)
         y = 0.0d0
         call xlm(l, x, y, z, w, maxl)
         do m = 0, l
            dp(k,m) = w(m,l)*qn(m,l) ! Condition with qn
         enddo
         z = cos(theta)
         x = sin(theta)
         y = 0.0d0
         call xlm(l, x, y, z, w, maxl)
         do m = 0, l
            a(k,m) = w(m,l)*qn(m,l) ! Condition with qn
         enddo
      enddo
c
      call dgesv(l+1,l+1,a,maxl+1,ipiv,dp,l+1,info)
      if (info .ne. 0) call errquit('xlm_y_rot_mat: Pinfo ', info,
     &       FMM_ERR)
c
c$$$      do m = 0, l
c$$$         dp(m,m) = dp(m,m) - 1.0d0
c$$$      enddo
c$$$      write(6,*) ' DP '
c$$$      call doutput(dp,1,l+1,1,l+1,l+1,l+1,1)
c$$$      call dfill((l+1)**2,0.0d0,dp,1)
c$$$      do m = 0, l
c$$$         dp(m,m) = dp(m,m) + 1.0d0
c$$$      enddo
c
      if (l .gt. 0) then
         phi = pi/dble(2*l+1)     ! -ve m components are zero at phi=0
         do k = 1, l
            theta = dble(k)*pi/dble(l+1)
            z = cos(theta)
            x = sin(theta)*cos(phi)
            y = sin(theta)*sin(phi)
            call xlm(l, x, y, z, w, maxl)
            do m = 1, l
               a(k,l-m+1) = w(-m,l)*qn(-m,l)
            enddo
c     
            rx = x*cos(-beta) - z*sin(-beta)
            ry = y
            rz = z*cos(-beta) + x*sin(-beta)
            call xlm(l, rx, ry, rz, w, maxl)
            do m = 1, l
               dm(k,l-m+1) = w(-m,l)*qn(-m,l)
            enddo
         enddo
c
         call dgesv(l,l,a(1,1),maxl+1,ipiv,dm,l,info)
         if (info .ne. 0) call errquit('xlm_y_rot_mat: Minfo ', info,
     &       FMM_ERR)
c
c$$$         do m = 1, l
c$$$            dm(m,m) = dm(m,m) - 1.0d0
c$$$         enddo
c$$$         write(6,*) ' DM '
c$$$         call doutput(dm,1,l,1,l,l,l,1)
c$$$         call dfill(l*l,0.0d0,dm,1)
c$$$         do m = 1, l
c$$$            dm(m,m) = dm(m,m) + 1.0d0
c$$$         enddo
      endif
c
      end
      subroutine xlm_translate2(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      double precision x, y, z
      integer lmax, lq
      double precision q(-lq:lq,0:lmax)
      integer lbig, l, m, mmp, mp, mmpa, lp, llp, mpa, lo, hi
      double precision sump, summ, bfac
c
c     given a set of unnormalized multipoles about the origin in q
c     translate the center to (x,y,z).
c
      if (lmax .gt. maxl) then
         write(6,*) ' parameter exceeded', lmax, maxl
         stop 1
      end if
c
      do l = lmax, 0, -1
         do m = -l, l
            if (abs(q(m,l)) .gt. 0.0d0) goto 10
         enddo
      enddo

 10   lbig = l
      if (lbig .lt. 0) then
         do l = 0, lmax
            do m = -l, l
               q(m,l) = 0.0d0
            enddo
         enddo
         return
      endif
c
      call xlm(lmax, -x, -y, -z, w, maxl)
c
c     transfer the charge, then everything else
c
      do l = 0, lmax
         do m = -l, l
*            q(m,l) = q(m,l) ! UH ?
*            w(m,l) = w(m,l)
            v(m,l) = w(m,l) * q(0,0) * b(abs(m),l) * b(0,0)
         end do
      end do
c
      do l = 1, lmax
         do m = 0, l
            sump = 0.0d0
            summ = 0.0d0
            do lp = 1, min(l,lbig)
               llp = l - lp
               hi = min(lp,max(l-lp+m,lp-l+m))
               lo = max(-lp,min(l-lp+m,lp-l+m,-lp))
               do mp = lo, hi
                  mmp = m - mp
                  mmpa = abs(mmp)
                     bfac = b(mmp,llp)*b(mp,lp)
                     mpa = abs(mp)
                     sump = sump + (w(mmpa,llp)*q(mpa,lp) - 
     $                    w(-mmpa,llp)*q(-mpa,lp)*s(mmp)*s(mp))*bfac
                     summ = summ + (w(mmpa,llp)*q(-mpa,lp)*s(mp) +
     $                    w(-mmpa,llp)*q(mpa,lp)*s(mmp))*bfac
               end do
            end do
            v(-m,l) = v(-m,l) + summ
            v( m,l) = v( m,l) + sump
         end do
      end do
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = v(m,l) / b(abs(m),l)
         end do
      end do
c
      end
      subroutine xlm_translate(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      double precision x, y, z
      integer lmax, lq
      double precision q(-lq:lq,0:lmax)
      integer lbig, l, m, mmp, mp, mmpa, lp, llp, mpa
      double precision sump, summ, bfac
c
c     given a set of unnormalized multipoles about the origin in q
c     translate the center to (x,y,z).
c
      if (lmax .gt. maxl) then
         write(6,*) ' parameter exceeded', lmax, maxl
         stop 1
      end if
c
      do l = lmax, 0, -1
         do m = -l, l
            if (abs(q(m,l)) .gt. 0.0d0) goto 10
         enddo
      enddo

 10   lbig = l
      if (lbig .lt. 0) then
         do l = 0, lmax
            do m = -l, l
               q(m,l) = 0.0d0
            enddo
         enddo
         return
      endif
c
      call xlm(lmax, -x, -y, -z, w, maxl)
c
c     transfer the charge, then everything else
c
      do l = 0, lmax
         do m = -l, l
*            q(m,l) = q(m,l) ! UH ?
*            w(m,l) = w(m,l)
            v(m,l) = w(m,l) * q(0,0) * b(abs(m),l) * b(0,0)
         end do
      end do
c
      do l = 1, lmax
         do m = 0, l
            sump = 0.0d0
            summ = 0.0d0
            do lp = 1, min(l,lbig)
               llp = l - lp
               do mp = -lp, lp
                  mmp = m - mp
                  mmpa = abs(mmp)
                  if (mmpa .le. llp) then
                     bfac = b(mmp,llp)*b(mp,lp)
                     mpa = abs(mp)
                     sump = sump + (w(mmpa,llp)*q(mpa,lp) - 
     $                    w(-mmpa,llp)*q(-mpa,lp)*s(mmp)*s(mp))*bfac
                     summ = summ + (w(mmpa,llp)*q(-mpa,lp)*s(mp) +
     $                    w(-mmpa,llp)*q(mpa,lp)*s(mmp))*bfac
                  end if
               end do
            end do
            v(-m,l) = v(-m,l) + summ
            v( m,l) = v( m,l) + sump
         end do
      end do
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = v(m,l) / b(abs(m),l)
         end do
      end do
c
      end
      subroutine xlm_local_translate(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      integer lmax, lq
      double precision x, y, z, q(-lq:lq,0:lmax)
c
c     Given a local spherical taylor series for the potential 
c     translate the center of the expansion to (x,y,z).  
c
c     If the routine xlm_translate does
c     .   Qml(r-a) <= sum(m'l')Tmlm'l'(-a)Qm'l'(r)
c     this routine does
c     .   Qml(r-a) <= sum(m'l')Qm'l'(r)Tm'l'ml(a)
c
      integer l, m, mp, mmp, mmpa, mpa, lp, llp, lo, hi
      double precision bfac
c
      if (lmax .gt. maxl) then
         write(6,*) ' parameter exceeded', lmax, maxl
         stop 1
      end if
c
      call xlm(lmax, x, y, z, w, maxl)
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = q(m,l) / b(abs(m),l)
            v(m,l) = 0.0d0
         end do
      end do
c
      do l = 0, lmax
         do m = 0, l
            do lp = 0, l
               llp = l - lp
*
               hi = min(lp,max(l-lp+m,lp-l+m))
               lo = max(-lp,min(l-lp+m,lp-l+m,-lp))
               do mp = lo, hi
*
*               do mp = -lp, lp
                  mmp = m - mp
                  mmpa = abs(mmp)
*                  if (mmpa .le. llp) then
                     bfac = b(mmp,llp)*b(mp,lp)
                     mpa = abs(mp)

*                     Tm,l,mpa,lp  +=  w(mmpa,llp)*bfac
*                     Tm,l,-mpa,lp += -w(-mmpa,llp)*s(mmp)*s(mp)*bfac
*                     T-m,l,-mpa,lp += w(mmpa,llp)*s(mp)*bfac
*                     T-m,l,mpa,lp  += w(-mmpa,llp)*s(mmp))*bfac

                     v( mpa,lp) = v( mpa,lp) + bfac*(
     $                    w(mmpa,llp)*q(m,l) + 
     $                    w(-mmpa,llp)*s(mmp)*q(-m,l))
                     v(-mpa,lp) = v(-mpa,lp) + bfac*(
     $                    w(mmpa,llp)*s(mp)*q(-m,l) - 
     $                    w(-mmpa,llp)*s(mmp)*s(mp)*q(m,l))

*                  end if
               end do
            end do
         end do
      end do
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = v(m,l)
         end do
      end do
c
      end
      double precision function xlm_coulomb(lmax, ax, ay, az, aq, la,
     $     bx, by, bz, bq, lb)
      implicit none
#include "xlm.fh"
c
      integer lmax, la, lb
      double precision ax, ay, az, bx, by, bz
      double precision aq(-la:la, 0:lmax), bq(-lb:lb, 0:lmax)
c
c     compute the coulomb interaction between two unnormalized
c     real solid multipole expansions at centers a and b
c
      double precision xx, yy, zz, rsq, r1, r12, energy, ei, ereal
      integer l1, l2, m1, m2, l2max, m1a, m2a, l12, m12, m12a
      double precision bcpr, bq_plus, bq_minus, aq_plus, aq_minus
      double precision qnorm, eimag
c
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(2*lmax, xx, yy, zz, u, maxl2)
c
      energy = 0.0d0
      ei = 0.0d0
      r1 = dsqrt(rsq)
      do l1 = 0, lmax
         l2max = lmax        ! Full summation
         do m1 = -l1, l1
            qnorm = abs(aq(m1,l1)) + abs(aq(-m1,l1))
            if (qnorm .ne. 0.0d0) then
               m1a = abs(m1)
               ereal = 0.0d0
               eimag = 0.0d0
               r12 = r1 * b(m1,l1) ! Included b here to keep no.s small
               do l2 = 0, l2max
                  do m2 = -l2, l2
                     m2a = abs(m2)
                     l12 = l1+l2
                     m12 = -(m1+m2)
                     m12a = abs(m12)
                     bcpr = (b(m2,l2)*r12)*(c(m12,l12)*phase(l2+m12))
                     bq_plus = (bq(m2a,l2)*bcpr)
                     bq_minus= (bq(-m2a,l2)*bcpr)
                     ereal = ereal + (bq_plus*u(m12a,l12) -
     $                    bq_minus*u(-m12a,l12)*s(m2)*s(m12))
                     eimag = eimag + (bq_minus*s(m2)*u(m12a,l12)+
     $                    bq_plus*u(-m12a,l12)*s(m12))
                  end do
                  r12 = r12 * rsq
               end do
               aq_plus = (aq(m1a,l1))
               aq_minus= (aq(-m1a,l1))
               energy = energy + 
     $              (ereal*aq_plus - eimag*aq_minus*s(m1))
               ei = ei + (ereal*aq_minus*s(m1)+eimag*aq_plus)
            end if
         end do
         r1 = r1 * rsq
      end do
*      if (abs(ei) .gt. 1.0d-10) then
*         write(6,*) ' xlmcoul: imaginary part of potential is large',
*     $        energy, ei
*      end if
c
      xlm_coulomb = energy
c
      end
      subroutine xlm_print(lmax, q, lq)
      implicit double precision (a-h, o-z)
      dimension q(-lq:lq, 0:lmax)
      parameter(mchunk = 6)
c
      do l = 0, lmax
         sum = 0.0d0
         do m = -l,l
            sum = sum + abs(q(m,l))
         end do
         if (sum .le. 1d-40) then
*            write(6,3) l
* 3          format(' l = ',i2,' is zero')
*            continue
         else
            write(6,1) l
 1          format(' l = ',i2)
            do mlow = -l, l, mchunk
               mhi = min(l, mlow+mchunk-1)
               write(6,2) (q(m,l),m=mlow,mhi)
 2             format('      ',1p,6d12.4)
            end do
         end if
      end do
c
      end
      double precision function xlm_potential(ax, ay, az, 
     $     bx, by, bz, bq, lb, lmax)
      implicit none
#include "xlm.fh"
c
      integer lb, lmax
      double precision ax, ay, az, bx, by, bz
      double precision bq(-lb:lb, 0:lmax)
c
c     compute the coulomb potential at A from the unnormalized 
c     real-solid multipole distribution at B.
c
      double precision xx, yy, zz, rsq
      double precision energy, r12
      integer l, m
c
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax, xx, yy, zz, u, maxl2)
c
      energy = 0.0d0
      r12 = dsqrt(rsq)
      do l = 0, lmax
         energy = energy + bcp0(0,l)*r12*bq(0,l)*u(0,l)
         do m = 1, l
            energy = energy + (bq(m,l)*u(m,l) +
     $           bq(-m,l)*u(-m,l))*bcp0(m,l)*r12
         end do
         r12 = r12 * rsq
      end do
c
      xlm_potential = energy
c
      end
      double precision function xlm_potential_test(ax, ay, az, 
     $     bx, by, bz, bq, lb, lmax)
      implicit none
#include "xlm.fh"
c
      integer lb, lmax
      double precision ax, ay, az, bx, by, bz
      double precision bq(-lb:lb, 0:lmax)
c
c     This version of xlm_potential just to make sure
c     the formula in the paper is correct!  It is.
c
      double precision xx, yy, zz, rsq
      double precision energy, r12, fourpi, sum
      integer l, m
c
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax, xx, yy, zz, u, maxl2)
      fourpi = 4.0d0*3.1415926535897932d0
c
      energy = 0.0d0
      r12 = dsqrt(rsq)
      do l = 0, lmax
         sum = 0.0d0
         do m = -l, l
            sum = sum + bq(m,l)*qn(m,l)*qn(m,l)*u(m,l)
         end do
         energy = energy + sum*fourpi*r12/dble(l+l+1)
         r12 = -r12 * rsq
      end do
c
      xlm_potential_test = energy
c
      end
      double precision function xlm_local_potential(x, y, z, 
     $     aq, la, lmax)
      implicit none
#include "xlm.fh"
      double precision x, y, z
      integer la, lmax
      double precision aq(-la:la,0:lmax)
c
c     Given in aq() an expansion of the local potential compute
c     it at the given point.  Note, that if the potential is 
c     expanded about some point A, then the coordinates relative
c     to A (r-a) should be passed in.
c
      integer l, m
      double precision sum
c
      call xlm(lmax, x, y, z, w, maxl)
      sum = 0.0d0
      do l = 0, lmax
         do m = -l,l
            sum = sum + w(m,l)*aq(m,l)
         end do
      end do
c
      xlm_local_potential = sum
c
      end
      subroutine xlm_unpack(lmax, p, q, lq)
      implicit none
      integer lmax, lq
      double precision p(*), q(-lq:lq, 0:lmax)
c
c     Unpack the multipoles in p (which are stored contiguously
c     with no space between) into the easier to use array q()
c
      integer l, m, ind
c
      ind = 1
      do l = 0, lmax
         do m = -l, l
            q(m,l) = p(ind)
            ind = ind + 1
         end do
      end do
c
      end
      subroutine xlm_pack(lmax, p, q, lq)
      implicit none
      integer lmax, lq
      double precision p(*), q(-lq:lq, 0:lmax)
c
c     Pack the multipoles into p (which are stored contiguously
c     with no space between) from the easier to use array q()
c
      integer l, m, ind
c
      ind = 1
      do l = 0, lmax
         do m = -l, l
            p(ind) = q(m,l)
            ind = ind + 1
         end do
      end do
c
      end
      subroutine xlm_accumulate_to_packed(lmax, p, q, lq)
      implicit none
      integer lmax, lq
      double precision p(*), q(-lq:lq, 0:lq)
c
c     Add into the packed multipoles p (which are stored contiguously
c     with no space between) from the easier to use array q()
c
      integer l, m, ind
c
      ind = 1
      do l = 0, lmax
         do m = -l, l
            p(ind) = p(ind) + q(m,l)
            ind = ind + 1
         end do
      end do
c
      end
      subroutine xlm_multipole_to_local(lmax, ax, ay, az, aq, la,
     $     bx, by, bz, bq, lb)
      implicit none
#include "xlm.fh"
c     
      integer lmax, la, lb
      double precision ax, ay, az, bx, by, bz
      double precision aq(-la:la, 0:la), bq(-lb:lb, 0:lb)
c     
c     Put into aq() the multipolar potential from the unnormalized
c     real solid multipole expansion at center B computed at center A.
c     
c     I.e., the potential at a point r near A due to the multipoles
c     at B may be computed as  
c     
c     .     sum(lm) Xlm(r-a)AQlm
c
c     
      double precision xx, yy, zz, rsq, r12, ereal, eimag
      integer l1, l2, m1, m2, l2max, l12, m12, m12a
      integer m1m2, m1mm2
      double precision urpp, urmm, uipp, uimm
**      integer flops
      double precision ww(0:maxl,0:maxl)
      double precision wv(0:maxl,0:maxl)
      double precision vw(0:maxl,0:maxl)
      double precision vv(0:maxl,0:maxl)
c
      do l1 = 0, la
         do m1 = -l1,l1
            aq(m1,l1) = 0.0d0
         enddo
      enddo
c
      do l2 = lmax,0,-1
         do m2 = -l2,l2
            if (abs(bq(m2,l2)) .gt. 0.0d0) goto 10
         enddo
      enddo
 10   l2max = l2
c
      if (l2max .lt. 0) return
c         
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax+l2max, xx, yy, zz, u, maxl2)
c
**      flops = 0
c
      do l2 = 0, l2max
         do m2 = 0, l2
            w(m2,l2) = bq( m2,l2)*b(m2,l2)*phase(l2+m2)*e(m2)
            v(m2,l2) = bq(-m2,l2)*b(m2,l2)*phase(l2+m2)*s(m2)*e(m2)
         enddo
         w(0,l2) = w(0,l2)*0.5d0
         v(0,l2) = v(0,l2)*0.5d0
      enddo
c
      do l12 = 0, lmax+l2max
         do m12 = -l12,l12
            m12a = abs(m12)
            ur(m12,l12) = u( m12a,l12)*c(-m12,l12)*e(m12)
            ui(m12,l12) = u(-m12a,l12)*c(-m12,l12)*s(-m12)*e(m12)
         end do
      end do
c
      r12 = sqrt(rsq)
      do l12 = 0, lmax+l2max
*         call dfill((maxl+1)**2,0.0d0,ww,1)
*         call dfill((maxl+1)**2,0.0d0,wv,1)
*         call dfill((maxl+1)**2,0.0d0,vw,1)
*         call dfill((maxl+1)**2,0.0d0,vv,1)
         do m1 = 0, min(l12,lmax)
            do m2 = 0, min(m1,l12-m1,l2max)
               m1m2 = m1+m2
               m1mm2 = m1-m2
               urpp = ur(m1m2, l12)*ee( m2,m1)
               urmm = ur(m1mm2,l12)*ee(-m2,m1)
               ww(m2,m1) = urpp + urmm
               ww(m1,m2) = ww(m2,m1)
               vv(m2,m1) = urpp - urmm
               vv(m1,m2) = vv(m2,m1)
c
               uipp = ui(m1m2, l12)*ee( m2,m1)
               uimm = ui(m1mm2,l12)*ee(-m2,m1)
               wv(m2,m1) = uimm - uipp
               vw(m1,m2) = -wv(m2,m1)
               vw(m2,m1) = uimm + uipp
               wv(m1,m2) = -vw(m2,m1)
            enddo
**            flops = flops + 8*(1+min(m1,l12-m1,l2max))
         enddo
c$$$         lll = min(l12+1,lmax+1)
c$$$         write(6,*) l12, ' ww '
c$$$         call doutput(ww, 1, lll, 1, lll, maxl+1, maxl+1, 1)
c$$$         write(6,*) l12, ' vv '
c$$$         call doutput(vv, 1, lll, 1, lll, maxl+1, maxl+1, 1)
c$$$         write(6,*) l12, ' wv '
c$$$         call doutput(wv, 1, lll, 1, lll, maxl+1, maxl+1, 1)
c$$$         write(6,*) l12, ' vw '
c$$$         call doutput(vw, 1, lll, 1, lll, maxl+1, maxl+1, 1)
         do l1 = max(0,l12-l2max), min(l12,lmax)
            l2 = l12 - l1
            ereal = ww(0,0)*w(0,l2)
            do m2 = 1, l2
               ereal = ereal+ww(m2,0)*w(m2,l2)+wv(m2,0)*v(m2,l2)
            end do
            aq( 0,l1) = aq( 0,l1) + ereal*r12
**            flops = flops + 3 + l2*4
            do m1 = 1, l1
               ereal = ww(0,m1)*w(0,l2)
               eimag = vw(0,m1)*w(0,l2)
               do m2 = 1, l2
                  ereal = ereal+ww(m2,m1)*w(m2,l2)+wv(m2,m1)*v(m2,l2)
                  eimag = eimag+vw(m2,m1)*w(m2,l2)+vv(m2,m1)*v(m2,l2)
               end do
               aq( m1,l1) = aq( m1,l1) + ereal*r12
               aq(-m1,l1) = aq(-m1,l1) - eimag*r12
            end do
**            flops = flops + l1*(6+l2*8)
         end do
         r12 = r12 * rsq
      end do
c
      do l1 = 0,lmax
         aq( 0,l1) = aq(0,l1)*b(0,l1)*0.5d0
         do m1 = 1, l1
            aq( m1,l1) = aq( m1,l1)*b(m1,l1)
            aq(-m1,l1) = aq(-m1,l1)*b(m1,l1)*s(m1)
         enddo
      enddo
c
c
**       write(6,*) ' oldflops ', flops
*       stop
c
      end
      subroutine xlm_multipole_to_local_z(lmax, az, aq, la, bq, lb)
      implicit none
#include "xlm.fh"
c     
      integer lmax, la, lb
      double precision az
      double precision aq(-la:la, 0:la), bq(-lb:lb, 0:lb)
c     
c     Put into aq() the multipolar potential from the unnormalized
c     real solid multipole expansion at center B computed at center A.
c
c     RESTRTICTED TO Z ONLY TRANSLATION
c     
c     I.e., the potential at a point r near A due to the multipoles
c     at B may be computed as  
c     
c     .     sum(lm) Xlm(r-a)AQlm
c     
      double precision r, r12, r12ww
      integer l1, l2, m1, l2max, l12, minl1l2
      double precision ww
c
      do l1 = 0, la
         do m1 = -l1,l1
            aq(m1,l1) = 0.0d0
         enddo
      enddo
c
c$$$      do l2 = lmax,0,-1
c$$$         do m2 = -l2,l2
c$$$            if (abs(bq(m2,l2)) .gt. 0.0d0) goto 10
c$$$         enddo
c$$$      enddo
c$$$ 10   l2max = l2
      l2max = lmax
c
c$$$      if (l2max .lt. 0) return
c
c$$$      do l2 = 0, l2max
c$$$         do m2 = -l2, l2
c$$$            w(m2,l2) = bq( m2,l2)*b(m2,l2)*phase(l2+m2)
c$$$         enddo
c$$$      enddo
c$$$      do l2 = 0, l2max
c$$$         do m2 = -l2, l2
c$$$            w(m2,l2) = bq(m2,l2)
c$$$         enddo
c$$$      enddo
c
      r = 1.0d0/az
      r12 = r
c
c     l12 = 0 
c
      ww = 2.0d0
      aq(0,0) = aq(0,0) + r12*ww*bq(0,0)
      r12 = r12 * r
c
c     l12 = 1
c
      if (lmax+l2max .ge. 1) then
         ww = -ww
         r12ww = r12*ww
         aq(0,0) = aq(0,0) + r12ww*bq(0,1) ! l1=0, l2=1
         aq(0,1) = aq(0,1) + r12ww*bq(0,0) ! l1=1, l2=0
         r12 = r12 * r
      endif
c
c     Can't this be done in L^2 ops?  
c
      do l12 = 2, lmax+l2max
         ww = -ww * dble(l12) ! (-1)^l * 2 * l!
         do l1 = max(0,l12-l2max), min(l12,lmax)
            l2 = l12 - l1
            minl1l2 = min(l1,l2)
            r12ww = r12*ww
            do m1 = -minl1l2, minl1l2
               aq( m1,l1) = aq( m1,l1) + r12ww*bq(m1,l2)
            end do
         end do
         r12 = r12 * r
      end do
c
c$$$      rsign = sign(1.0d0,r)
c$$$      do l1 = 0,lmax
c$$$         aq( 0,l1) = aq(0,l1)*b(0,l1)*0.5d0*rsign
c$$$         do m1 = 1, l1
c$$$            aq( m1,l1) = aq( m1,l1)*(b(m1,l1)*rsign)
c$$$            aq(-m1,l1) = aq(-m1,l1)*(b(m1,l1)*rsign)*phase(m1)
c$$$         enddo
c$$$      enddo
c
      end
      subroutine xlm_multipole_to_local_swap(lmax,
     $     ax, ay, az, aq, ap, la,
     $     bx, by, bz, bq, bp, lb)
      implicit none
#include "xlm.fh"
c     
      integer lmax, la, lb
      double precision ax, ay, az, bx, by, bz
      double precision aq(-la:la, 0:la), bq(-lb:lb, 0:lb)
      double precision ap(-la:la, 0:la), bp(-lb:lb, 0:lb)
c     
c     Put into ap() the multipolar potential from the unnormalized
c     real solid multipole expansion at center B computed at center A
c     and similarly pu the potential from aq() into bp()
c     
c     I.e., the potential at a point r near A due to the multipoles
c     at B may be computed as  
c     
c     .     sum(lm) Xlm(r-a)APlm
c
c     and vice versa
c     
      double precision xx, yy, zz, rsq, r12, r12ph
      double precision ereala, eimaga, erealb, eimagb 
      integer l1, l2, m1, m2, l12, m12, m12a
      integer m1m2, m1mm2, lmax2
      double precision urpp, urmm, uipp, uimm, ce, bpe
      integer flops
      double precision ww(0:maxl,0:maxl)
      double precision wv(0:maxl,0:maxl)
      double precision vw(0:maxl,0:maxl)
      double precision vv(0:maxl,0:maxl)
      double precision wa(0:maxl,0:maxl), va(0:maxl,0:maxl)
      double precision wb(0:maxl,0:maxl), vb(0:maxl,0:maxl)
c
      lmax2 = lmax+lmax
c
      do l1 = 0, la
         do m1 = -l1,l1
            ap(m1,l1) = 0.0d0
         enddo
      enddo
      do l1 = 0, lb
         do m1 = -l1,l1
            bp(m1,l1) = 0.0d0
         enddo
      enddo
c         
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax2, xx, yy, zz, u, maxl2)
c
      flops = 0
c
      do l2 = 0, lmax
         do m2 = 0, l2
            bpe = b(m2,l2)*phase(l2+m2)*e(m2)
            wa(m2,l2) = aq( m2,l2)*bpe
            va(m2,l2) = aq(-m2,l2)*(bpe*s(m2))
            wb(m2,l2) = bq( m2,l2)*bpe
            vb(m2,l2) = bq(-m2,l2)*(bpe*s(m2))
         enddo
         wa(0,l2) = wa(0,l2)*0.5d0
         wb(0,l2) = wb(0,l2)*0.5d0
      enddo
c
      do l12 = 0, lmax2
         do m12 = -l12,l12
            m12a = abs(m12)
            ce = c(-m12,l12)*e(m12)
            ur(m12,l12) = u( m12a,l12)*ce
            ui(m12,l12) = u(-m12a,l12)*ce*s(-m12)
         end do
      end do
c
      r12   = sqrt(rsq)
      r12ph = r12
      do l12 = 0, lmax2
         do m1 = 0, min(l12,lmax)
            do m2 = 0, min(m1,l12-m1,lmax)
               m1m2 = m1+m2
               m1mm2 = m1-m2
               urpp = ur(m1m2, l12)*ee( m2,m1)
               urmm = ur(m1mm2,l12)*ee(-m2,m1)
               ww(m2,m1) = urpp + urmm
               ww(m1,m2) = ww(m2,m1)
               vv(m2,m1) = urpp - urmm
               vv(m1,m2) = vv(m2,m1)
c
               uipp = ui(m1m2, l12)*ee( m2,m1)
               uimm = ui(m1mm2,l12)*ee(-m2,m1)
               wv(m2,m1) = uimm - uipp
               vw(m1,m2) = -wv(m2,m1)
               vw(m2,m1) = uimm + uipp
               wv(m1,m2) = -vw(m2,m1)
            enddo
*            flops = flops + 8*(1+min(m1,l12-m1,lmax))
         enddo
         do l1 = max(0,l12-lmax), min(l12,lmax)
            l2 = l12 - l1
            ereala = ww(0,0)*wa(0,l2)
            erealb = ww(0,0)*wb(0,l2)
            do m2 = 1, l2
               ereala = ereala + ww(m2,0)*wa(m2,l2)+wv(m2,0)*va(m2,l2)
               erealb = erealb + ww(m2,0)*wb(m2,l2)+wv(m2,0)*vb(m2,l2)
            end do
            ap( 0,l1) = ap( 0,l1) + erealb*r12
            bp( 0,l1) = bp( 0,l1) + ereala*r12ph
*            flops = flops + 3 + l2*4
            do m1 = 1, l1
               ereala = ww(0,m1)*wa(0,l2)
               eimaga = vw(0,m1)*wa(0,l2)
               do m2 = 1, l2
                  eimaga =eimaga+vw(m2,m1)*wa(m2,l2)+vv(m2,m1)*va(m2,l2)
                  ereala =ereala+ww(m2,m1)*wa(m2,l2)+wv(m2,m1)*va(m2,l2)
               enddo
               bp( m1,l1) = bp( m1,l1) + ereala*r12ph
               bp(-m1,l1) = bp(-m1,l1) - eimaga*r12ph
               erealb = ww(0,m1)*wb(0,l2)
               eimagb = vw(0,m1)*wb(0,l2)
               do m2 = 1, l2    ! SUN slower with merged loops ?
                  erealb =erealb+ww(m2,m1)*wb(m2,l2)+wv(m2,m1)*vb(m2,l2)
                  eimagb =eimagb+vw(m2,m1)*wb(m2,l2)+vv(m2,m1)*vb(m2,l2)
               end do
               ap( m1,l1) = ap( m1,l1) + erealb*r12
               ap(-m1,l1) = ap(-m1,l1) - eimagb*r12
            end do
*            flops = flops + l1*(6+l2*8)
         end do
         r12 = r12 * rsq
         r12ph = -r12ph * rsq
      end do
c
      do l1 = 0,lmax
         ap( 0,l1) = ap( 0,l1)*b(0,l1)*0.5d0
         bp( 0,l1) = bp( 0,l1)*b(0,l1)*0.5d0
         do m1 = 1, l1
            ap( m1,l1) = ap( m1,l1)*b(m1,l1)
            ap(-m1,l1) = ap(-m1,l1)*b(m1,l1)*s(m1)
            bp( m1,l1) = bp( m1,l1)*b(m1,l1)
            bp(-m1,l1) = bp(-m1,l1)*b(m1,l1)*s(m1)
         enddo
      enddo
c
c
*       write(6,*) ' flops ', flops
*       stop
c
      end
      subroutine xlm_new_multipole_to_local(lmax, 
     $     xp, yp, zp, p, lp, xq, yq, zq, q, lq)
      implicit none
#include "errquit.fh"
#include "xlm.fh"
c     
      double precision flops, flopsv, flopst, calls
      common /xlmflops/flops, flopsv, flopst, calls
c     
      integer lmax, lq, lp
      double precision xp, yp, zp, xq, yq, zq
      double precision q(-lq:lq,0:lmax), p(-lp:lp,0:lp)
c     
      double precision theta, phi, x, y, z, zz, phisign, qqp, qqm
      double precision rr, r12, ww, r12ww, pp, pm
      double precision pi, thetasign
      double precision qq(-maxl:maxl,0:maxl), tmp(-maxl:maxl)
      integer l, m, mp, ihash, iphi, itheta, indp, indm
      integer l1, l2, l12, m1, minl1l2
c     
c     From the multipole-potential q() at (xq,yq,zq) compute
c     the local Taylor series potential at (xp,yp,zp), in p() 
c     using a rotation-based algorithm.  
c     
c     Under the assumption that only a limited no. of angles
c     will be used, cache rotation matrix info.
c     
      x = xp - xq
      y = yp - yq
      z = zp - zq
      rr = 1.0d0/sqrt(x*x + y*y + z*z)
c     
c     Determine the required rotations
c     
      pi = 3.1415926535897932d0
      if (x.eq.0.0d0 .and. y.eq.0.0d0) then
         phi = 0.0d0
         theta = 0.0d0
         if (z .lt. 0.0d0) theta = pi
      else 
         zz = z*rr
         if (abs(zz) .gt. 1.0d0) zz = sign(1.0d0,zz)
         theta = acos(zz)
         phi = atan2(y,x)
      endif
c     
c     Rotate the source multipoles.  The actual potential
c     is roughly 2*qml*bml*cml*xml/r^(2l+1), so before we
c     rotate we must include these additional factors.
c     
c     Hash into the list of cached angles to get the rotation info.
c     The cached phis are all positive so need to track the sign.
c     The cached thetas are all 0-pi/2, so if phi is pi/2-pi then
c     must rotate with theta-pi and use negative translation.
c     
      ihash = int(abs(phi*1000.0d0))
      if (ihash .gt. maxhash) call errquit('xlm_n_m_l: phi hash',ihash,
     &       FMM_ERR)
      if (iphis(ihash) .gt. 0) then
         if (abs(phis(ihash)-abs(phi)) .gt. 1d-12) then
            write(6,*) ' ihash phi phis ', ihash, phi, phis(ihash)
            call errquit('xlm_n_m_l: duplicate phi hash',0,
     &       FMM_ERR)
         endif
      else
         call xlm_cache_phi(phi)
      endif
      iphi = iphis(ihash)
      phisign = sign(1.0d0,phi)
c     
      if (theta .le. 0.5d0*pi) then
         thetasign = 1.0d0
      else
         theta = theta - pi
         thetasign = -1.0d0
      endif
      ihash = int(abs(theta*1000.0d0))
      if (ihash .gt. maxhash)call errquit('xlm_n_m_l: theta hash',ihash,
     &       FMM_ERR)
      if (ithetas(ihash) .gt. 0) then
         if (abs(thetas(ihash)-abs(theta)) .gt. 1d-12) then
            write(6,*) ' ihash theta thetas ', ihash, theta, 
     $           thetas(ihash)
            call errquit('xlm_n_m_l: duplicate theta hash',0, FMM_ERR)
         endif
      else
         call xlm_cache_theta(theta)
      endif
      itheta = ithetas(ihash)
c     
c     Add extra weights, apply Z rotation, apply preconditioning
c     
      if (phisign .gt. 0.0d0) then
         do l = 0, lmax
            do m = 1, l
               qqp = cosmphi(m,iphi)*q( m,l) + sinmphi(m,iphi)*q(-m,l)
               qqm = cosmphi(m,iphi)*q(-m,l) - sinmphi(m,iphi)*q( m,l)
               qq( m,l) = qqp*twobc(m,l)
               qq(-m,l) = qqm*twobc(-m,l)
            enddo
            qq(0,l) = q(0,l)*rqn(0,l)
         enddo
      else
         do l = 0, lmax
            do m = 1, l
               qqp = cosmphi(m,iphi)*q( m,l) - sinmphi(m,iphi)*q(-m,l)
               qqm = cosmphi(m,iphi)*q(-m,l) + sinmphi(m,iphi)*q( m,l)
               qq( m,l) = qqp*twobc(m,l)
               qq(-m,l) = qqm*twobc(-m,l)
            enddo
            qq(0,l) = q(0,l)*rqn(0,l)
         enddo
      endif
*      flops = flops + 1 + 5*lmax + 4*lmax*lmax
c     
c     Apply Y rotation
c     
      indp = 2
      indm = 1
      do l = 1, lmax
         do m = -l, l
            tmp(m) = qq(m,l)
            qq(m,l) = 0.0d0
         enddo
c     
*         flopsv = flopsv + l*l*2 + (l+1)*(l+1)*2
c     
         if (l .gt. 3) then
            call xlm_m2l_rot1(l, qq(-l,l), tmp(-l),
     $           cdp(indp,itheta), cdm(indm,itheta), thetasign)
            indp = indp + (l+1)*(l+1)
            indm = indm + l*l
         else
            if (thetasign .gt. 0.0d0) then
               do mp = 0, l
                  do m = 0, l
                     qq(m,l) = qq(m,l) + cdp(indp+m,itheta)*tmp(mp)
                  enddo
                  indp = indp + l + 1
               enddo
               do mp = -l,-1
                  indm = indm + l
                  do m = -l,-1
                     qq(m,l) = qq(m,l) + cdm(indm+m,itheta)*tmp(mp)
                  enddo
               enddo
            else
               do m = 0, l
                  do mp = 0, l
                     qq(m,l) = qq(m,l) + cdp(indp+mp,itheta)*tmp(mp)
                  enddo
                  indp = indp + l + 1
               enddo
               do m = -l,-1
                  indm = indm + l
                  do mp = -l,-1
                     qq(m,l) = qq(m,l) + cdm(indm+mp,itheta)*tmp(mp)
                  enddo
               enddo
            endif
         endif
      enddo
c     
c     Remove precond and weights and add extra weights for translate
c
      do l = 0, lmax
         do m = -l, l
            qq(m,l) = qq(m,l)*rtwobcqnbp(m,l)
            p(m,l) = 0.0d0
         enddo
      enddo
*      flops = flops + (lmax+1)*(lmax+1)
c     
c     Translate 0,0,r
c
      rr = thetasign * rr
      r12 = rr
c
c     l12 = 0 
c
      ww = 2.0d0
      p(0,0) = p(0,0) + r12*ww*qq(0,0)
      r12 = r12 * rr
c
*      flopst = flopst + 3d0
c
c     l12 = 1
c
      if (lmax .gt. 0) then
         ww = -ww
         r12ww = r12*ww
         p(0,0) = p(0,0) + r12ww*qq(0,1) ! l1=0, l2=1
         p(0,1) = p(0,1) + r12ww*qq(0,0) ! l1=1, l2=0
         r12 = r12 * rr
*         flopst = flopst + 7d0
      endif
c
      do l12 = 2, lmax+lmax
         ww = -ww * dble(l12) ! (-1)^l * 2 * l12!
         do l1 = max(0,l12-lmax), min(l12,lmax)
            l2 = l12 - l1
            minl1l2 = min(l1,l2)
            r12ww = r12*ww
            do m1 = -minl1l2, minl1l2
               p( m1,l1) = p( m1,l1) + r12ww*qq(m1,l2)
            end do
*            flopst = flopst + 4*minl1l2 + 5d0
         end do
         r12 = r12 * rr
*         flopst = flopst + 1
      end do
c
      if (thetasign .gt. 0.0d0) then
         do l = 0, lmax
            do m = -l, l
               p(m,l) = p(m,l) * brqnp(m,l)
            enddo
         enddo
      else
         do l = 0, lmax
            do m = -l, l
               p(m,l) = -p(m,l) * brqnp(m,l)
            enddo
         enddo
      endif
*      flops = flops + 2*lmax + 1
c     
c     Rotate the potentials back
c     
c     Precondition
c
      indp = 2
      indm = 1
      do l = 1, lmax
c     
c     Apply Y rotation
c     
         do m = -l, l
            tmp(m) = p(m,l)
            p(m,l) = 0.0d0
         enddo
c     
*         flopsv = flopsv + l*l*2 + (l+1)*(l+1)*2
c$$$         if (l .gt. -1) then
            call xlm_m2l_rot1(l, p(-l,l), tmp(-l),
     $           cdp(indp,itheta), cdm(indm,itheta), -thetasign)
            indp = indp + (l+1)*(l+1)
            indm = indm + l*l
c$$$         else
c$$$            if (thetasign .lt. 0.0d0) then
c$$$               do mp = 0, l
c$$$                  do m = 0, l
c$$$                     p(m,l) = p(m,l) + cdp(indp+m,itheta)*tmp(mp)
c$$$                  enddo
c$$$                  indp = indp + l + 1
c$$$               enddo
c$$$               do mp = -l,-1
c$$$                  indm = indm + l
c$$$                  do m = -l,-1
c$$$                     p(m,l) = p(m,l) + cdm(indm+m,itheta)*tmp(mp)
c$$$                  enddo
c$$$               enddo
c$$$            else
c$$$               do m = 0, l
c$$$                  do mp = 0, l
c$$$                     p(m,l) = p(m,l) + cdp(indp+mp,itheta)*tmp(mp)
c$$$                  enddo
c$$$                  indp = indp + l + 1
c$$$               enddo
c$$$               do m = -l,-1
c$$$                  indm = indm + l
c$$$                  do mp = -l,-1
c$$$                     p(m,l) = p(m,l) + cdm(indm+mp,itheta)*tmp(mp)
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$         endif
      enddo
c
c     Undo preconditioning and apply Z rotation
c
      if (phisign .gt. 0.0d0) then
         do l = 0, lmax
            p(0,l) = p(0,l) * qn(0,l)
            do m = 1, l
               pp = p( m,l)*qn( m,l)
               pm = p(-m,l)*qn(-m,l)
               qqp = cosmphi(m,iphi)*pp - sinmphi(m,iphi)*pm
               qqm = cosmphi(m,iphi)*pm + sinmphi(m,iphi)*pp
               p( m,l) = qqp
               p(-m,l) = qqm
            enddo
         enddo
      else
         do l = 0, lmax
            p(0,l) = p(0,l) * qn(0,l)
            do m = 1, l
               pp = p( m,l)*qn( m,l)
               pm = p(-m,l)*qn(-m,l)
               qqp = cosmphi(m,iphi)*pp + sinmphi(m,iphi)*pm
               qqm = cosmphi(m,iphi)*pm - sinmphi(m,iphi)*pp
               p( m,l) = qqp
               p(-m,l) = qqm
            enddo
         enddo
      endif
*      flops = flops + 1 + 5*lmax + 4*lmax*l
c     
*      calls = calls + 1
c     
      end
      subroutine xlm_cache_phi(phi)
      implicit none
#include "errquit.fh"
#include "xlm.fh"
      double precision phi
c
c     Add to the list of phis with cached rotation info
c
      integer m, ihash
c
      if (numphi.lt.0 .or. numphi.ge.maxnumphi) call errquit
     $     ('xlm_cache_phi: numphi is invalid ', numphi, FMM_ERR)
      numphi = numphi + 1
c
      ihash = int(abs(phi*1000.0d0))
      if (iphis(ihash) .gt. 0) call errquit
     $     ('xlm_cache_phi: duplicate hash ', ihash, FMM_ERR)
      iphis(ihash) = numphi
      phis(ihash) = abs(phi)
c
      do m = 1, maxl
         cosmphi(m,numphi) = cos(dble(m)*abs(phi))
         sinmphi(m,numphi) = sin(dble(m)*abs(phi))
      end do
c
*      write(6,1) abs(phi), ihash, numphi
* 1    format(' xlm_cache_phi: added phi=',f9.6,
*     $     '   ihash=',i5,'   iphi=',i4)
c
      end
      subroutine xlm_cache_theta(theta)
      implicit none
#include "errquit.fh"
#include "xlm.fh"
      double precision theta
c
c     Add to the list of thetas with cached rotation info
c
      integer l, ihash, indp, indm
c
      if (numtheta.lt.0 .or. numtheta.ge.maxnumtheta) call errquit
     $     ('xlm_cache_theta: numtheta is invalid ', numtheta, FMM_ERR)
      numtheta = numtheta + 1
c
      ihash = int(abs(theta*1000.0d0))
      if (ithetas(ihash) .gt. 0) call errquit
     $     ('xlm_cache_theta: duplicate hash ', ihash, FMM_ERR)
      ithetas(ihash) = numtheta
      thetas(ihash) = abs(theta)
c
      indp = 1
      indm = 1
c$$$      write(6,*) ' theta ', theta
      do l = 0, maxl
         call xlm_y_rotation_matrix(l, abs(theta),
     $        cdp(indp,numtheta), cdm(indm,numtheta))
c$$$         if (l .lt. 4) then
c$$$            write(6,*)  l
c$$$            call output(cdp(indp,numtheta),1,l+1,1,l+1,l+1,l+1,1)
c$$$         endif
         indp = indp + (l+1)*(l+1)
         indm = indm + l*l
      end do
      if (indm-1 .ne. dmdim) call errquit
     $     ('xlm_theta_cache: bad indm ', indm, FMM_ERR)
      if (indp-1 .ne. dpdim) call errquit
     $     ('xlm_theta_cache: bad indp ', indp, FMM_ERR)
c
*      write(6,1) abs(theta), ihash, numtheta
* 1    format(' xlm_cache_theta: added theta=',f9.6,
*     $     '   ihash=',i5,'   itheta=',i4)
c
      end
c$$$      subroutine xlm_m2l_rot1(l,qq,tmp,cdp,cdm,thetasign)
c$$$      implicit none
c$$$      integer l, lmax, maxl
c$$$      double precision qq(-l:l), tmp(-l:l)
c$$$      double precision cdp(0:l,0:l), cdm(-l:-1,-l:-1)
c$$$      double precision thetasign
c$$$c     
c$$$      integer m, mp
c$$$c     
c$$$      if (thetasign .gt. 0.0d0) then
c$$$         do mp = 0, l
c$$$            do m = 0, l
c$$$               qq(m) = qq(m) + cdp(m,mp)*tmp(mp)
c$$$            enddo
c$$$         enddo
c$$$         do mp = -l,-1
c$$$            do m = -l,-1
c$$$               qq(m) = qq(m) + cdm(m,mp)*tmp(mp)
c$$$            enddo
c$$$         enddo
c$$$      else
c$$$         do m = 0, l
c$$$            do mp = 0, l
c$$$               qq(m) = qq(m) + cdp(mp,m)*tmp(mp)
c$$$            enddo
c$$$         enddo
c$$$         do m = -l,-1
c$$$            do mp = -l,-1
c$$$               qq(m) = qq(m) + cdm(mp,m)*tmp(mp)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$c
c$$$      end

      
C****************************************************************
C   Translated by Pacific-Sierra Research VAST-2          
C   Version 6.1C1 on  4/ 3/99 at 11:38:11
C****************************************************************
C
      subroutine xlm_m2l_rot1(l,qq,tmp,cdp,cdm,thetasign)
      implicit none
      integer l
      double precision qq(-l:l), tmp(-l:l)
      double precision cdp(0:l,0:l), cdm(-l:-1,-l:-1)
      double precision thetasign
c      
      integer m, mp
c      
      integer j1, j2, j3, j4, mp1, mp2, m1, m2, m3, m4
      doubleprecision d1, d2, d3, d4, d5, qq1, qq2, qq3, qq4, qq5, d11, 
     .d12, d13, d14, d15, qq6, qq7, qq8, qq9, qq10, qq61, qq62, qq63, qq
     .11, qq12, qq13, d27, d28, d29, d30, d31, d32, d33, d34, d35, d36, 
     .d37, d38, d39, d40, d41, d42
      if (thetasign .gt. 0.0d0) then
*      if (l + 1 .gt. 0) then
         j1 = iand(l + 1,3)
         do mp = 1, j1
            m4 = iand(max0(l + 1,0),3)
            do m = 1, m4
               qq(m-1) = qq(m-1) + cdp(m-1,mp-1)*tmp(mp-1)
            end do
            do m = m4 + 1, l + 1, 4
               d39 = qq(m-1) + cdp(m-1,mp-1)*tmp(mp-1)
               d40 = qq(m) + cdp(m,mp-1)*tmp(mp-1)
               d41 = qq(m+1) + cdp(m+1,mp-1)*tmp(mp-1)
               d42 = qq(m+2) + cdp(m+2,mp-1)*tmp(mp-1)
               qq(m-1) = d39
               qq(m) = d40
               qq(m+1) = d41
               qq(m+2) = d42
            end do
         end do
         do mp = j1 + 1, l + 1, 4
            m3 = iand(max0(l + 1,0),3)
            do m = 1, m3
               qq(m-1) = qq(m-1) + cdp(m-1,mp-1)*tmp(mp-1) + cdp(m-1,mp)
     1            *tmp(mp) + cdp(m-1,mp+1)*tmp(mp+1) + cdp(m-1,mp+2)*tmp
     2            (mp+2)
            end do
            do m = m3 + 1, l + 1, 4
               d35 = qq(m-1) + cdp(m-1,mp-1)*tmp(mp-1) + cdp(m-1,mp)*tmp
     1            (mp) + cdp(m-1,mp+1)*tmp(mp+1) + cdp(m-1,mp+2)*tmp(mp+
     2            2)
               d36 = qq(m) + cdp(m,mp-1)*tmp(mp-1) + cdp(m,mp)*tmp(mp)
     1             + cdp(m,mp+1)*tmp(mp+1) + cdp(m,mp+2)*tmp(mp+2)
               d37 = qq(m+1) + cdp(m+1,mp-1)*tmp(mp-1) + cdp(m+1,mp)*tmp
     1            (mp) + cdp(m+1,mp+1)*tmp(mp+1) + cdp(m+1,mp+2)*tmp(mp+
     2            2)
               d38 = qq(m+2) + cdp(m+2,mp-1)*tmp(mp-1) + cdp(m+2,mp)*tmp
     1            (mp) + cdp(m+2,mp+1)*tmp(mp+1) + cdp(m+2,mp+2)*tmp(mp+
     2            2)
               qq(m-1) = d35
               qq(m) = d36
               qq(m+1) = d37
               qq(m+2) = d38
            end do
         end do
*      endif
*      if (l .gt. 0) then
         j2 = iand(l,3)
         do mp = 1, j2
            m2 = iand(max0(l,0),3)
            do m = 1, m2
               qq(m-1-l) = qq(m-1-l) + cdm(m-1-l,mp-1-l)*tmp(mp-1-l)
            end do
            do m = m2 + 1, l, 4
               d31 = qq(m-1-l) + cdm(m-1-l,mp-1-l)*tmp(mp-1-l)
               d32 = qq(m-l) + cdm(m-l,mp-1-l)*tmp(mp-1-l)
               d33 = qq(m+1-l) + cdm(m+1-l,mp-1-l)*tmp(mp-1-l)
               d34 = qq(m+2-l) + cdm(m+2-l,mp-1-l)*tmp(mp-1-l)
               qq(m-1-l) = d31
               qq(m-l) = d32
               qq(m+1-l) = d33
               qq(m+2-l) = d34
            end do
         end do
         do mp = j2 + 1, l, 4
            m1 = iand(max0(l,0),3)
            do m = 1, m1
               qq(m-1-l) = qq(m-1-l) + cdm(m-1-l,mp-1-l)*tmp(mp-1-l) + 
     1            cdm(m-1-l,mp-l)*tmp(mp-l) + cdm(m-1-l,mp+1-l)*tmp(mp+1
     2            -l) + cdm(m-1-l,mp+2-l)*tmp(mp+2-l)
            end do
            do m = m1 + 1, l, 4
               d27 = qq(m-1-l) + cdm(m-1-l,mp-1-l)*tmp(mp-1-l) + cdm(m-1
     1            -l,mp-l)*tmp(mp-l) + cdm(m-1-l,mp+1-l)*tmp(mp+1-l) + 
     2            cdm(m-1-l,mp+2-l)*tmp(mp+2-l)
               d28 = qq(m-l) + cdm(m-l,mp-1-l)*tmp(mp-1-l) + cdm(m-l,mp-
     1            l)*tmp(mp-l) + cdm(m-l,mp+1-l)*tmp(mp+1-l) + cdm(m-l,
     2            mp+2-l)*tmp(mp+2-l)
               d29 = qq(m+1-l) + cdm(m+1-l,mp-1-l)*tmp(mp-1-l) + cdm(m+1
     1            -l,mp-l)*tmp(mp-l) + cdm(m+1-l,mp+1-l)*tmp(mp+1-l) + 
     2            cdm(m+1-l,mp+2-l)*tmp(mp+2-l)
               d30 = qq(m+2-l) + cdm(m+2-l,mp-1-l)*tmp(mp-1-l) + cdm(m+2
     1            -l,mp-l)*tmp(mp-l) + cdm(m+2-l,mp+1-l)*tmp(mp+1-l) + 
     2            cdm(m+2-l,mp+2-l)*tmp(mp+2-l)
               qq(m-1-l) = d27
               qq(m-l) = d28
               qq(m+1-l) = d29
               qq(m+2-l) = d30
            end do
         end do
*      endif
      else
*      if (l + 1 .gt. 0) then
         j3 = iand(l + 1,3)
         do m = 1, j3
            qq1 = qq(m-1)
            mp2 = iand(max0(l + 1,0),3)
            do mp = 1, mp2
               qq1 = qq1 + cdp(mp-1,m-1)*tmp(mp-1)
            end do
            qq11 = 0.D0
            qq12 = 0.D0
            qq13 = 0.D0
            do mp = mp2 + 1, l + 1, 4
               qq1 = qq1 + cdp(mp-1,m-1)*tmp(mp-1)
               qq11 = qq11 + cdp(mp,m-1)*tmp(mp)
               qq12 = qq12 + cdp(mp+1,m-1)*tmp(mp+1)
               qq13 = qq13 + cdp(mp+2,m-1)*tmp(mp+2)
            end do
            qq1 = qq1 + qq11 + qq12 + qq13
            qq(m-1) = qq1
         end do
cDEC$ NOVECTOR
         do m = j3 + 1, l + 1, 4
            qq2 = qq(m-1)
            qq3 = qq(m)
            qq4 = qq(m+1)
            qq5 = qq(m+2)
            do mp = 1, l + 1
               d1 = tmp(mp-1)
               d1 = tmp(mp-1)
               d1 = tmp(mp-1)
               d1 = tmp(mp-1)
               d2 = qq2 + cdp(mp-1,m-1)*d1
               d3 = qq3 + cdp(mp-1,m)*d1
               d4 = qq4 + cdp(mp-1,m+1)*d1
               d5 = qq5 + cdp(mp-1,m+2)*d1
               qq2 = d2
               qq3 = d3
               qq4 = d4
               qq5 = d5
            end do
            qq(m+2) = qq5
            qq(m+1) = qq4
            qq(m) = qq3
            qq(m-1) = qq2
         end do
*      endif
*      if (l .gt. 0) then
         j4 = iand(l,3)
         do m = 1, j4
            qq6 = qq(m-1-l)
            mp1 = iand(max0(l,0),3)
            do mp = 1, mp1
               qq6 = qq6 + cdm(mp-1-l,m-1-l)*tmp(mp-1-l)
            end do
            qq61 = 0.D0
            qq62 = 0.D0
            qq63 = 0.D0
            do mp = mp1 + 1, l, 4
               qq6 = qq6 + cdm(mp-1-l,m-1-l)*tmp(mp-1-l)
               qq61 = qq61 + cdm(mp-l,m-1-l)*tmp(mp-l)
               qq62 = qq62 + cdm(mp+1-l,m-1-l)*tmp(mp+1-l)
               qq63 = qq63 + cdm(mp+2-l,m-1-l)*tmp(mp+2-l)
            end do
            qq6 = qq6 + qq61 + qq62 + qq63
            qq(m-1-l) = qq6
         end do
cDEC$ NOVECTOR
         do m = j4 + 1, l, 4
            qq7 = qq(m-1-l)
            qq8 = qq(m-l)
            qq9 = qq(m+1-l)
            qq10 = qq(m+2-l)
            do mp = 1, l
               d11 = tmp(mp-1-l)
               d11 = tmp(mp-1-l)
               d11 = tmp(mp-1-l)
               d11 = tmp(mp-1-l)
               d12 = qq7 + cdm(mp-1-l,m-1-l)*d11
               d13 = qq8 + cdm(mp-1-l,m-l)*d11
               d14 = qq9 + cdm(mp-1-l,m+1-l)*d11
               d15 = qq10 + cdm(mp-1-l,m+2-l)*d11
               qq7 = d12
               qq8 = d13
               qq9 = d14
               qq10 = d15
            end do
            qq(m+2-l) = qq10
            qq(m+1-l) = qq9
            qq(m-l) = qq8
            qq(m-1-l) = qq7
         end do
*      endif
      endif
c      
      end
