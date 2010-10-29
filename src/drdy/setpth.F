C $Id$
      Subroutine drdy_setpth
C
C Setup for first step off saddle point
C
      Implicit None
#include "errquit.fh"
c
#include "drdyP.fh"
#include "stdio.fh"
*::local::
      Integer i,ierr,ilrg
      Double Precision u1max
C
      write(fu6,600)
C Calculation for step off saddle point
C Diagonalize hessian matrix stored in f
      Call drdy_diagf (n3tm,n3,fspsv,hess,vec0,u0,vec1,vec2,ierr)
      if (ierr.ne.0) then
         write(fu6,6000) ierr
         call errquit('drdy_setpth:drdy_diagf error',555, INPUT_ERR)
      endif
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
      write(fu6,601) (vec0sv(i),i=1,n3)
C           
      if (lgs(31).ge.3) then
C Calculate quadratic correction to step off saddle point
         call drdy_pmlqsp
         write(fu6,602) (vec1sv(i),i=1,n3)
      endif
      Return
600   Format(//,' Set up to take first step off saddle point')
601   Format(/,' Vector for initial step off saddle point=',/,
     *   (T35,0p9f9.4))
602   Format(/,' Second order correction for vector=',/,
     *   (T35,0p9f9.4))
6000  Format (' in subroutine diagf error in call to eispack routine',
     *   ' rsp, ierr=', I5)
      End
