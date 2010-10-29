C $Id$
      Subroutine drdy_pmcube (n3tm,n3,x,dx,fmat,fold,amass,dels,delsf,
     *    vec0,vec1,vec2)
C
C Given x, dx, fmat, and fold compute next step along the MEP using
C   the Page-McIver CUBE algorithm
C
C   x     - vector of mass-weighted coordinates
C   dx    - vector of mass-weighted gradient
C   fmat  - hessian matrix in mass-weighted coordinates
C   fold  - hessian matrix in mass-weighted coordinates from previous step
C   amass - vector of the sqrt root of the ratio of the atomic and reduced masses
C   dels  - step size in s
C   delsf - step size in s for forward difference derivative of f
C   vec0,vec1,vec2 - scratch arrays
C   ierr  - error flag
C
      Implicit none
      Integer n3tm,n3
      Double Precision x(n3tm),dx(n3tm),fmat(n3tm,n3tm),fold(n3tm,n3tm),
     *   amass(n3tm),dels,delsf,vec0(n3tm),vec1(n3tm),vec2(n3tm)
C
      Integer i,ierr,j
      Double Precision sum,sum1,sum2,vfv1,vfv,vgvv,xnorm
C
      ierr = 0
C Normalize gradient
      xnorm = 0.0d00
      do i = 1,n3
         xnorm = xnorm + dx(i)*dx(i)
      enddo
      xnorm = sqrt(xnorm)
      do i = 1,n3
         vec0(i) = -dx(i)/xnorm
      enddo
C Construct v0*f*v0 and v1 vector (vec1)
      vfv = 0.0d00
      do i = 1,n3
         sum = 0.0d00
         do j = 1,n3
            sum = sum + fmat(i,j)*vec0(j)
         enddo
         vec1(i) = -sum
         vfv = vfv + sum*vec0(i)
      enddo
      do i = 1,n3
         vec1(i) = (vec1(i) + vfv*vec0(i))/xnorm
      enddo
C use backward difference to approximate G0v0 by (fmat-fold)/delsf
      vgvv = 0.0d00
      vfv1 = 0.0d00
      do i = 1,n3
         sum1 = 0.0d00
         sum2 = 0.0d00
         do j = 1,n3
            sum1 = sum1 + (fmat(i,j)-fold(i,j))*vec0(j)/delsf
            sum2 = sum2 + fmat(i,j)*vec1(j)
         enddo
         vec2(i) = -sum1 - sum2
         vgvv = vgvv + sum1*vec0(i)
         vfv1 = vfv1 + sum2*vec0(i)
      enddo
      do i = 1,n3
         vec2(i) = (vec2(i) + 2.0d00*vfv*vec1(i) + 
     *    (2.0d00*vfv1+vgvv)*vec0(i))/xnorm
      enddo
C use Taylor series to take calculate next geometry
      do i = 1,n3
         x(i) = x(i) + dels*(vec0(i) + 0.5d00*dels*(vec1(i) + 
     *    dels*vec2(i)/3.0d00))
      enddo
      return
      end
