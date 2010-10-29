C $Id$
      Subroutine drdy_diagf (n3tm,n3,fmat,hess,xlam,u0,
     &      scr1,scr2,ierr)
*
      Implicit None
*::passed::
      Integer n3tm,n3,ierr
      Double Precision fmat(n3tm,n3tm),hess(1),xlam(n3tm),u0(n3tm,n3tm),
     * scr1(n3tm),scr2(n3tm)
*::local::
      Integer i,ii,j,nbar
C      
      ii = 0
      do i = 1,n3
         do j = 1,i
            ii = ii + 1
            hess(ii) = fmat(j,i)
         enddo
      enddo
C Diagonalize F, eigenvalues in XLAM, eigenvectors in U0
      nbar = (n3*(n3+1))/2
      Call drdy_rsp (n3tm,n3,nbar,hess,xlam,1,u0,scr1,scr2,ierr)
      return
      end
