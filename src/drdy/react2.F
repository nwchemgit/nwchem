C $Id$
      Subroutine drdy_react2 (iop)
C
C  put correct parts of hessian matrix into fsv, calculate frequencies,
C   and call rphwrt
C
      Implicit None 
C
#include "drdyP.fh"
*::passed::
      Integer iop
*::local::
      Integer i,ierr,ifreq,ii,indx1,indx2,ip,ixyz,j,jndx1,jndx2,jxyz,
     * ndim,nfreq,nfreq1
C
      if (icode(iop).eq.1) then
         if (iop.le.2) write(fu6,600) iop
         if (iop.gt.2) write(fu6,601) iop-2
      else
         indx1 = 0
         do i = 1,nratom(iop)
            indx2 = 3*(iatsv(i,iop)-1)
            do ixyz = 1,3
               indx1 = indx1 + 1
               indx2 = indx2 + 1
               jndx1 = 0
               do j = 1,nratom(iop)
                  jndx2 = 3*(iatsv(j,iop)-1)
                  do jxyz = 1,3
                     jndx1 = jndx1 + 1
                     jndx2 = jndx2 + 1
                     fsv(jndx1,indx1) = f(jndx2,indx2)
                  enddo
               enddo
            enddo
         enddo
C
C  Diagonalize hessian and write out frequencies and eigenvectors
         ndim = 3*nratom(iop)
         ii = 0
         do i = 1,ndim
            do ip = 1,i
               ii = ii + 1
               hess(ii) = fsv(ip,i)
            enddo
         enddo
         if (iop.le.2) write(fu6,602) iop
         if (iop.gt.2) write(fu6,603) iop-2
         nfreq = nf(iop)
         call drdy_fdiag(ndim,nfreq,ierr,1)
         if (ierr.eq.0) then
C Store frequencies in freqr array
            nfreq = nf(iop)
            do ifreq = 1,nfreq
               freqr(ifreq,iop) = scr2(ifreq)
            enddo
            if (iop.eq.2 .or. iop.eq.4) then
               nfreq1 = nf(iop-1)
               do ifreq = 1,nfreq
                  freqr(ifreq+nfreq1,iop-1) = scr2(ifreq)
               enddo
            endif
         endif
      endif
C  Write out information to file 30
      call drdy_rphwrt(iop)
      return
600   Format(/,' Reactant species',i5,' is an atom')
601   Format(/,' Product species ',i5,' is an atom')
602   Format(/,' Reactant species',i5)
603   Format(/,' Product species ',i5)
      end
