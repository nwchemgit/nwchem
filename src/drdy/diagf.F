C $Id: diagf.F,v 1.2 2001-01-08 15:51:59 windus Exp $
      Subroutine drdy_diagf (n3tm,n3,fmat,hess,xlam,u0,
     &      scr1,scr2,ierr,ipr)
*
      Implicit None ! Double Precision (A-H,O-Z)
*::passed::
      Integer ierr, ipr, n3tm, n3
      Double Precision fmat(n3tm,n3tm),hess(1),xlam(n3tm),u0(n3tm,n3tm)
      Double Precision scr1(n3tm),scr2(n3tm)
*::local::
      Integer i, ii, j, nbar
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
      if (ierr .ne. 0) then
         if (ipr.ne.0) write (6,6000) ierr
      else
         if(ipr.ne.0) write (6,668) (xlam(i),i=1,n3)
         if(ipr.ne.0) write (6,*) ' eigenvectors='
         do i = 1,n3
            if(ipr.ne.0) write (6,665) (u0(i,j),j=1,n3)
         enddo
      endif
      return
665   Format (1x,1p5e15.7)
668   Format (/,' Diagonalize hessian',/,' eigenvalues=',/,
     *   (1x,1p5e15.7))
6000  Format (' in subroutine diagf error in call to eispack routine',
     *   ' rsp, ierr=', I5)
      end
