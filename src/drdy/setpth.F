      Subroutine drdy_setpth
C
C Setup for first step off saddle point
C
      Implicit None ! Double Precision (A-H,O-Z)
c
#include "drdyP.fh"
#include "stdio.fh"
*::local::
      Integer fu6
      Integer i, ierr, ilrg
      Double Precision u1max
*
      fu6 = luout
C
      Write (fu6,600)
C Calculation for step off saddle point
C Diagonalize hessian matrix stored in f
      Call drdy_diagf (n3tm,n3,fspsv,hess,vec0,u0,vec1,vec2,ierr,1)
C It is assumed that the smallest (most negative) eigenvalue is ordered first.
C    Determine phase of the eigenvector so that largest component is positive.
      u1max = 0.0d00
      do i = 1,n3
         if (abs(u0(i,1)) .ge. u1max) then
            u1max = abs(u0(i,1))
            ilrg = i
         end if
      enddo
C Switch sign if (largest element)*dir < 0
      if (u0(ilrg,1)*dir .lt. 0.0d00) then
         do i = 1,n3
            u0(i,1) = -u0(i,1)
         enddo
      endif
C Save initial direction vector
      do i = 1, n3
         vec0sv(i) = u0(i,1)
         vec1sv(i) = 0.0d00
      enddo
      write (fu6,601) (vec0sv(i),i=1,n3)
C           
      if (lgs(31).ge.3) then
C Calculate quadratic correction to step off saddle point
         call drdy_pmlqsp
         write (fu6,602) (vec1sv(i),i=1,n3)
      endif
      Return
600   Format(//,' Set up to take first step off saddle point')
601   Format(/,' Vector for initial step off saddle point=',/,
     *   (1x,1p6e13.5))
602   Format(/,' Second order correction for vector=',/,
     *   (1x,1p6e13.5))
      End
