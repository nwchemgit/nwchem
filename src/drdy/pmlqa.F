C $Id$
      Subroutine drdy_pmlqa (n3tm,n3,ndbar,x,dx,fmat,hess,amass,dels,
     *    vec0,vec1,u0,xlam,alpha,ierr)
C
C Given x, dx, fmat, and fold compute next step along the MEP using
C   the Page-McIver LQA algorithm
C
C   x     - vector of mass-weighted coordinates
C   dx    - vector of mass-weighted gradient
C   fmat  - hessian matrix in mass-weighted coordinates
C   amass - vector of the sqrt root of the ratio of the atomic and reduced masses
C   dels  - step size in s
C   vec0,vec1,u0,xlam,alpha - scratch arrays
C   ierr  - error flag
C
      Implicit none
#include "errquit.fh"
C
      Integer n3tm,n3,ndbar,ierr
      Double Precision x(n3),dx(n3),fmat(n3tm,n3),hess(ndbar),
     *   amass(n3),dels,vec0(n3),vec1(n3),u0(n3tm,n3),xlam(n3),
     *   alpha(n3)
C
      Integer i,ic1,ic2,j
      Double Precision arg,delt,ft,ftlast,s,slast,sum,t,ts,xnorm
C
      Integer icmax
      Data icmax/1000000/
      Integer fu6
      Data fu6/6/
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
C      write (fu6,*) ' in pmlqa, fmat='
C      write (fu6,*) ((fmat(j,i),j=1,i),i=1,n3)
      Call drdy_diagf
     &      (n3tm,n3,fmat,hess,xlam,u0,alpha,vec1,ierr)
      if (ierr.ne.0) then
         write(fu6,6000) ierr
         call errquit('drdy_pmlqa:drdy_diagf error',555, UNKNOWN_ERR)
      endif
C Page-McIver LQA algorithm
C Construct u0(transpose)*vec0 and put into vec1, and scale xlam by xnorm
      do i = 1,n3
         xlam(i) = xlam(i)/xnorm
         vec1(i) = 0.0d00
         do j = 1,n3
            vec1(i) = vec1(i) + u0(j,i)*vec0(j)
         enddo
      enddo
C      write (fu6,668) (xlam(i),i=1,n3)
C      write (fu6,669) (vec1(i),i=1,n3)
C Integrate expression that relates s to t
      delt = 0.2d00*dels
      t = 1.0d10
      ts = 0.0d00
      ic1 = 0
      ic2 = 0
C Loop over step size for integration (delt) until converged or maximum
C   iterations
      do while (abs(1.0d00-ts/t).gt.1.d-6.and.ic1.le.15.and.
     *          ic2.lt.icmax)
         ic1 = ic1 + 1
         ts = t
         delt = 0.5d00*delt
         t = 0.0d00
C Get integrand at t=0
         Call drdy_dsdt (n3,vec1,xlam,t,ft)
         s = 0.0d00
         ic2 = 0
C Integration loop, integrate until s>dels using trapezoidal rule
         do while (s.lt.dels .and. ic2.le.icmax)
            ic2 = ic2 + 1
            slast = s
            ftlast = ft
            t = t + delt
            Call drdy_dsdt (n3,vec1,xlam,t,ft)
            s = s + 0.5d00*delt*(ft + ftlast)
         enddo
C Integration complete, linearly interpolate t to s=dels
         t = t - (dels-s)*delt/(slast-s)
C         write (fu6,670) ic1,ic2,delt,t,1.0-ts/t
      enddo
      if (abs(1.0d00-t/ts).le.1.d-6) then
C Construct alpha*g0'
         do i = 1,n3
            arg = xlam(i)*t
            if (abs(arg).lt.1.d-8) then
               alpha(i) = -t*(1.0d00-arg*(0.5d00-arg/6.0d00))*vec1(i)
            else
               alpha(i) = (exp(-arg) - 1.0d00)*vec1(i)/xlam(i)
            endif
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
         write (fu6,6002) ic1,ic2
         ierr = 1
      end if
      return
668   format (' xlam=',/,(1x,1p5e15.7))
669   format (' go=',/,(1x,1p5e15.7))
670   format (' ic1,ic2,delt,t,err=',2i10,1p3e15.7)
6000  Format (' in subroutine diagf error in call to eispack routine',
     *   ' rsp, ierr=', I5)
6002  FORMAT (' In subroutine PMLQA, Page-McIver LQA algorithm,',
     * ' problem with integration of s(t)',/,
     * ' iterations, ic1,ic2=',2i10)
      end
