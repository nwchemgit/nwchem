C $Id: pmlqa.F,v 1.3 2001-01-08 15:51:59 windus Exp $
      Subroutine drdy_pmlqa (n3tm,n3,ndbar,x,dx,fmat,hess,amass,dels,
     *    vec0,vec1,u0,xlam,alpha,ierr)
C
C Given x, dx, fmat, and fold compute next step along the MEP using
C   the Page-McIver CUBE algorithm
C
C   x     - vector of mass-weighted coordinates
C   dx    - vector of mass-weighted gradient
C   fmat  - hessian matrix in mass-weighted coordinates
C   amass - vector of the sqrt root of the ratio of the atomic and reduced masses
C   dels  - step size in s
C   vec0,vec1,u0,xlam,alpha - scratch arrays
C   ierr  - error flag
C
      Implicit Double Precision (A-H,O-Z)
      Dimension x(n3),dx(n3),fmat(n3tm,n3),hess(ndbar),amass(n3),
     *   vec0(n3),vec1(n3),u0(n3tm,n3),xlam(n3),alpha(n3)
C
      ierr = 0
C Normalize gradient
      xnorm = 0.0d00
      do i = 1,n3
         xnorm = xnorm + dx(i)*dx(i)
      enddo
      xnorm = sqrt(xnorm)
      do i = 1,n3
         vec0(i) = dx(i)/xnorm
      enddo
C      write (6,*) ' in pmlqa, fmat='
C      write (6,*) ((fmat(j,i),j=1,i),i=1,n3)
      Call drdy_diagf
     &      (n3tm,n3,fmat,hess,xlam,u0,alpha,vec1,ierr,0)
C Page-McIver LQA algorithm
C Construct u0(transpose)*vec0 and put into vec1, and scale xlam by xnorm
      do i = 1,n3
         xlam(i) = xlam(i)/xnorm
         vec1(i) = 0.0d00
         do j = 1,n3
            vec1(i) = vec1(i) + u0(j,i)*vec0(j)
         enddo
      enddo
C      write (6,669) (vec1(i),i=1,n3)
C Integrate expression that relates s to t
      delt = 0.2d00*dels
      t = 1.0d10
      ts = 1.0d00
      ic1 = 0
C Loop over step size for integration (delt) until converged or maximum
C   iterations
      do while (abs(1.0d0-t/ts).gt.1.d-6 .and. ic1.le.15) 
         ic1 = ic1 + 1
         ts = t
         delt = 0.5d00*delt
         t = 0.0d00
C Get integrand at t=0
         Call drdy_dsdt (n3,vec1,xlam,t,ft)
         s = 0.0d00
         ic2 = 0
C Integration loop, integrate until s>dels using trapezoidal rule
         do while (s.lt.dels .and. ic2.le.10000)
            ic2 = ic2 + 1
            slast = s
            ftlast = ft
            t = t + delt
            Call drdy_dsdt (n3,vec1,xlam,t,ft)
            s = s + 0.5d00*delt*(ft + ftlast)
         enddo
C Integration complete, linearly interpolate t to s=dels
         t = t - (dels-s)*delt/(slast-s)
C         write (6,670) ic1,ic2,delt,t,1.0-ts/t
      enddo
      if (ic1.le.15 .and. ic2.le.10000) then
C Construct alpha*g0'
         do i = 1,n3
            alpha(i) = (exp(-xlam(i)*t) - 1.0d00)*vec1(i)/xlam(i)
         enddo
C  Take step and unscale XLAM
         do i = 1,n3
            xlam(i) = xnorm*xlam(i)
            sum = x(i)
            do j = 1,n3
               sum = sum + u0(i,j)*alpha(j)
            enddo
            x(i) = sum
         enddo
      else
         write (6,6002)
         ierr = 1
      end if
      return
669   format (' go=',/,(1x,1p5e15.7))
670   format (' ic1,ic2,delt,t,err=',2i10,1p3e15.7)
6002  FORMAT (' IN SUBROUTINE LQA, PAGE-MCIVER LQA ALGORITHM', /,
     * ' PROBLEM WITH INTEGRATION OF S(T)')
      end
