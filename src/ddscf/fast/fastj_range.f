      double precision function fastj_gaussian_range(n, alpha, eps)
*
* $Id$
*
      implicit none
c
      integer n
      double precision alpha, eps
c
c     Return an approximation to the outer solution of 
c     .     r^n*exp(-ar^2) = eps
c     .     r = (n*ln(-ln(eps)) - n*ln(a) - 4*ln(eps)) /
c     .         4*sqrt(-alpha*ln(eps))
c
c     Accuracy improves with smaller eps.
c
      double precision logeps
c
      logeps = log(eps)
c      
      fastj_gaussian_range = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      end

      double precision function fnlm(n,l,m,alpha,a,b,c,x,y,z)
      implicit none
c
      integer n, l, m
      double precision alpha, a, b, c, x, y, z
c
c     r^(n-l) Xlm exp(-a*r*r) centered about (a,b,c)
c     
      double precision rsq, r, radial, q(-10:10,0:10)
c
      rsq = (x-a)**2 + (y-b)**2 + (z-c)**2
      r   = sqrt(rsq)
c
      radial = exp(-alpha*rsq)
      if (n .ne. l) radial = radial * r**(n-l)
c
      call xlm(l, x-a, y-b, z-c, q, 10)
c
      fnlm = radial * q(m,l)
c
*      if (r .lt. 0.5d0) then
*         write(6,*) n, l, m, alpha, a, b, c, x, y, z, r
*         call xlm_print(l, q, 10)
*         write(6,*) radial, q(m,l), radial
*      end if
c
**      fnlm = exp(-alpha*rsq) * (z-c)
c      
      end
      double precision function vnlm(n,l,m,alpha,a,b,c,x,y,z)
      implicit none
c
      integer n, l, m
      double precision alpha, a, b, c, x, y, z
c
c     Return the potential due to 
c         r^(n-l) Xlm exp(-a*r*r) centered about (a,b,c)
c
c     A(sqrt(a)*r) * 4*pi/(a**((n+2)/2)) * Xlm / r^l
c
c     A(x) = (0.5*x**(n-l) (n+l+1)F((n+l)/2,x**2) + (n-l)**I(n-l-1)/2)
c     .         / (2l+1)
c
c     
      double precision rsq, r, radial, q(-10:10,0:10)
      double precision anl_fit
      external anl
c
      rsq = (x-a)**2 + (y-b)**2 + (z-c)**2
      r   = sqrt(rsq)
c
      radial = anl_fit(
     $     n,l,sqrt(alpha)*r) * 4.0d0 * 3.1415926535897932d0 /
     $     (alpha**(dble(n-l+2)/2.0d0))
c
      call xlm(l, x-a, y-b, z-c, q, 10)
      vnlm = radial * q(m,l)

**      vnlm = radial * (z-c)
c
      end
      double precision function fastj_r_neglected(k, alpha, eps)
      implicit none
c
      integer k
      double precision alpha, eps
c
c     For a function f(r) = r^k*exp(-alpha*r^2) determine
c     the radial distance r such that the fraction of the 
c     function norm that is neglected if the 3D volume 
c     integration is terminated at a distance r is less
c     than or equal to eps.
c
c     The returned value is accurate to about 0.01.
c
      double precision r, test, fastj_neglected, step
c
      r = sqrt(-log(eps)/alpha)
      step = 0.125d0*r
 10   test = fastj_neglected(k,alpha,r)
**      write(6,*) ' r error ', r, test
      if (test .gt. eps) then
         r = r + step
      else 
         r = r - step
         if (r .lt. 0.0d0) r = 0.0d0
         step = step*0.5d0
         r = r + step
      endif
      if (step .gt. 0.01d0) goto 10
c
      fastj_r_neglected = r
c
      end
      double precision function fastj_neglected(k,alpha,r)
      implicit none
c
      integer k
      double precision alpha, r
c
c     For a function f(r) = r^k*exp(-alpha*r^2) determine
c     the fraction of the function norm that is neglected
c     if the 3D volume integration is terminated at a 
c     distance r.
c
c     neglected = int(t^2*f(t),t=r..infinity)/int(t^2*f(t),t=0..infinity)
c
      double precision fastj_ik
c
      fastj_neglected = fastj_ik(k+2,alpha,r) / 
     $     fastj_ik(k+2,alpha,0.0d0)
c
      end
      double precision function fastj_ik(k,alpha,r)
      implicit none
c
      integer k
      double precision alpha, r
c
c     I(k) = int(t^k exp(-alpha*t^2), t=0..infinity)
c
c     I(k) = [(k-1)*I(k-2) + r^(k-1)*exp(-alpha*r^2)]/(2*alpha)
c
      integer i, ilo
      double precision value
c
      double precision erfc
c
      ilo = mod(k,2)
c
      if (ilo .eq. 0) then
         value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $        erfc(sqrt(alpha)*r)
      else
         value = exp(-alpha*r*r)/(2.0d0*alpha)
      endif
c
      do i = ilo+2,k,2
         value = ((i-1)*value + r**(i-1)*exp(-alpha*r*r)) / 
     $        (2.0d0*alpha)
      enddo
c
      fastj_ik = value
c
      end
