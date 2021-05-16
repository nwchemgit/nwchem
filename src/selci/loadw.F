C> \ingroup selci
C> @{
      subroutine selci_loadw(iflwmt, w1, w2, nfmax, nfmax2, nsmax)
*
* $Id$
*
#include "implicit.fh"
      dimension w1(nfmax,nfmax,nsmax), w2(nfmax2,nfmax,nsmax-1)
c
c     read in the coupling coefficients off wmatrx
c
      read (iflwmt,*) 
     &     (((w1(iu,iv,i),iu=1,nfmax),iv=1,nfmax),i=nsmax,1,-1)
      if (nfmax2*nfmax*(nsmax-1).gt.0) then
         read (iflwmt,*)
     &        (((w2(iu,iv,i),iu=1,nfmax2),iv=1,nfmax),i=1,nsmax-1)
      endif
c
c      write(6,*) ' w1 '
c      do 10 i = 1,nsmax
c         write(6,*) i
c         call output(w1(1,1,i),1,nfmax,1,nfmax,nfmax,nfmax,1)
c 10   continue
c
c      write(6,*) ' w2 '
c      do 20 i = 1,nsmax-1
c         write(6,*) i
c         call output(w2(1,1,i),1,nfmax2,1,nfmax,nfmax2,nfmax,1)
c 20   continue
c
      end
C> @}
