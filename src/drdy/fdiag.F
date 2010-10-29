C $Id$
      Subroutine drdy_fdiag (ndim,nfreq,ierr,jprnt)
C
      Implicit none
      Integer ndim,nfreq,ierr,jprnt
C
#include "drdyP.fh"
C
      Integer i,j,jp,m,m1,m2,mp,n,nbar
      Double Precision temp
C
C  Diagonalize hessian and write out frequencies and eigenvectors
      n=1
      nbar = (ndim*(ndim+1))/2
      call drdy_RSP (n3tm,ndim,nbar,hess,scr1,n,eigvec,scr2,scr3,ierr)
      if (ierr.eq.0) then
C Sort eigenvalues in scr1 into descending order, and shift locations
C  of eigenvectors to match
         m = 2
         do while (m.le.ndim)
            temp = scr1(m)
            do i = 1,ndim
               scr2(i) = eigvec(i,m)
            enddo
            j = 1
            do while (j.lt.m .and. temp.lt.scr1(j))
               j = j + 1
            enddo
            jp = m
            do while (jp.gt.j)
               scr1(jp) = scr1(jp-1)
               do i = 1,ndim
                  eigvec(i,jp) = eigvec(i,jp-1)
               enddo
               jp = jp - 1
            enddo
            scr1(j) = temp
            do i = 1,ndim
               eigvec(i,j) = scr2(i)
            enddo
            m = m + 1
         enddo
C look for large negative eigenvalues to replace small positive or negative ones
         j = nfreq
         jp = ndim
         do while (jp.gt.nfreq+1 .and. j.ge.1
     *    .and. -scr1(jp).gt.abs(scr1(j)))
            temp = scr1(jp)
            scr1(jp) = scr1(j)
            scr1(j) = temp
            do i = 1,ndim
               temp = eigvec(i,jp)
               eigvec(i,jp) = eigvec(i,j)
               eigvec(i,j) = temp
            enddo
            j = j - 1
            jp = jp - 1
         enddo
C Calculate frequencies
         do m = 1,ndim
            scr2(m) = autocm*scr1(m)
            if (scr2(m).ne.0.0d00) scr2(m) = 
     *       scr2(m)/(sqrt(redm*abs(scr1(m))))
         enddo
         if (jprnt.ne.0) then
            write(fu6,600)
            m1 = 1
            m2 = min(6,ndim)
            do while (m1.le.ndim)
               Write(fu6,601) (m,m=m1,m2)
               Write(fu6,602) (scr1(m),m=m1,m2)
               Write(fu6,603) (scr2(m),m=m1,m2)
               Write(fu6,604)
               do m = 1,ndim
                  Write(fu6,605) m,(eigvec(m,mp), mp=m1,m2)
               enddo
               Write(fu6,*) ' '
               m1 = m2+1
               m2 = min(m1+5,ndim)
            enddo
         endif
      else
         write(fu6,606)
      endif
C
      return
600   Format(/,' Eigenvalues and eigenvectors of hessian')
601   Format(/,23x,'m=',i8,5i13)
602   Format(2x,'Eigenvalues (au)=',T26,6f13.6)
603   Format(2x,'Frequencies (cm^-1)=',T26,6f13.6)
604   Format(2x,'Eigenvectors=')
605   Format(20x,i5,6f13.6)
606   Format(/,' error in call to rsp, frequencies not calculated')
      end
