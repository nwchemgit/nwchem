C $Id$
      Subroutine drdy_pmclqa (n3tm,n3,ndbar,x,dx,fmat,fold,hess,amass,
     *    dels,delsf,vec0,vec1,u0,xlam,alpha,ierr)
C
C Given x, dx, fmat, and fold compute next step along the MEP using
C   the Page-McIver CLQA algorithm
C
C   x     - vector of mass-weighted coordinates
C   dx    - vector of mass-weighted gradient
C   fmat  - hessian matrix in mass-weighted coordinates
C   fold  - hessian matrix in mass-weighted coordinates from previous step
C   amass - vector of the sqrt root of the ratio of the atomic and reduced masses
C   dels  - step size in s
C   delsf - step size in s for forward difference derivative of f
C   vec0,vec1,u0,xlam,alpha - scratch arrays
C   ierr  - error flag
C
      Implicit none
      Integer n3tm,n3,ndbar,ierr
      Double Precision x(n3),dx(n3),fmat(n3tm,n3),fold(n3tm,n3),
     * hess(ndbar),amass(n3),dels,delsf,vec0(n3),vec1(n3),u0(n3tm,n3),
     * xlam(n3),alpha(n3)
C
      Integer i,j
      Double Precision sum,vgvv,xnorm
C
      Call drdy_pmlqa (n3tm,n3,ndbar,x,dx,fmat,hess,amass,dels,
     *    vec0,vec1,u0,xlam,alpha,ierr)
C Compute correction term
C   First calculate normalized force vector
      xnorm = 0.0d00
      do i = 1,n3
         xnorm = xnorm + dx(i)*dx(i)
      enddo
      xnorm = sqrt(xnorm)
      do i = 1,n3
         vec0(i) = dx(i)/xnorm
      enddo
C use backward difference to approximate G0v0 by (fmat-fold)/delsf
      vgvv = 0.0d00
      do i = 1,n3
         sum = 0.0d00
         do j = 1,n3
            sum = sum + (fmat(i,j)-fold(i,j))*vec0(j)/delsf
         enddo
         vec1(i) = -sum
         vgvv = vgvv + sum*vec0(i)
      enddo
      do i = 1,n3
         vec1(i) = (vec1(i) + vgvv*vec0(i))/xnorm
         x(i) = x(i) + vec1(i)*(dels**3)/6.0d00
      enddo
      return
      end
