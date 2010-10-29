C $Id$
      Subroutine drdy_moiten (x,amat,apack,eig,vec,scr1,scr2,lcol,
     *   natom,ndim)
C
C Construct moment of inertia tensor and diagonalize it
C
      Implicit none
#include "errquit.fh"
      Integer natom,ndim
      Logical lcol
      Double Precision x(ndim),amat(3,3),apack(ndim),
     *   eig(3),vec(3,3),scr1(3),scr2(3)
C
      Integer i,ierr,j,jmax,k,koff,ndima,nbar
      Double Precision eps,r1mag,vecmx,xc,yc,zc
      Data eps/1.e-6/
C
      Integer fu6
      Data fu6/6/
C
C      write(fu6,*) ' entered moiten with natom,ndim='
C      write(fu6,*) natom,ndim
C      write(fu6,*) ' x='
C      write(fu6,*) (x(i),i=1,ndim)
      do j = 1,3
         do k = 1,3
            amat(j,k) = 0.0
         enddo
      enddo
      do i = 1,natom
         koff = 3*(i-1)
         xc = x(koff+1)
         yc = x(koff+2)
         zc = x(koff+3)
         amat(1,1) = amat(1,1) + yc*yc + zc*zc
         amat(2,2) = amat(2,2) + xc*xc + zc*zc
         amat(3,3) = amat(3,3) + xc*xc + yc*yc
         amat(1,2) = amat(1,2) - xc*yc
         amat(1,3) = amat(1,3) - xc*zc
         amat(2,3) = amat(2,3) - yc*zc
      enddo
      amat(2,1) = amat(1,2)
      amat(3,1) = amat(1,3)
      amat(3,2) = amat(2,3)
      i = 0
      do j = 1,3
         do k = 1,j
            i = i + 1
            apack(i) = amat(j,k)
         enddo
      enddo
C      write(fu6,*) ' moment of inertia ='
C      write(fu6,*) ((amat(j,i),i=1,3),j=1,3)
C      write(fu6,*) ' packed ='
C      write(fu6,*) (apack(i),i=1,6)
C diagonalize the moment of inertia tensor to get the rr1 matrix
      ndima = 3
      nbar = 6
      call drdy_rsp(ndima,ndima,nbar,apack,eig,1,vec,scr1,scr2,ierr)
      if (ierr .ne. 0) then
         write(fu6,6000) ierr
         call errquit('drdy_moiten:problem with drdy_rsp',555,
     &       UNKNOWN_ERR)
      endif
C Change phase of eigenvectors so largest component is positive
      do i = 1,3
         vecmx = abs(vec(1,i))
         jmax = 1
         do j = 2,3
            if(abs(vec(j,i)).gt.vecmx) then
               jmax = j
               vecmx = abs(vec(j,i))
            endif
         enddo
         if (vec(jmax,i).lt.0.0) then
            do j = 1,3
               vec(j,i) = -vec(j,i)
            enddo
         endif
      enddo
C      write(fu6,*) ' eigenvalues of tensor='
C      write(fu6,*) (eig(i),i=1,3)
C      write(fu6,*) ' eigenvectors='
C      write(fu6,*) ((vec(j,i),i=1,3),j=1,3)
C test for collinear orientation (zero eigenvalue)
      lcol = abs(eig(1)).lt.eps
      if (lcol) then
C Construct normalized vector from atom 1 to 2 for orientation 1
         r1mag = 0.0
         do k = 1,3
            vec(k,1) = x(k+3)-x(k)
            r1mag = r1mag + vec(k,1)*vec(k,1)
         enddo
         r1mag = sqrt(r1mag)
         do k = 1,3
            vec(k,1) = vec(k,1)/r1mag
         enddo
C         write(fu6,*) ' collinear, vec(k,1)='
C         write(fu6,*) (vec(k,1),k=1,3)
      endif
      return
6000  Format (' in subroutine moiten, error in call to eispack routine',
     *   ' rsp, ierr=', I5)
      end
