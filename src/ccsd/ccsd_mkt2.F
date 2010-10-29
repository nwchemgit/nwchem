      subroutine ccsd_mkt2(ncor,nocc,nvir,nact,eorb,
     $                     t2,scra,g_nexch,g_t2,iprt)
C     $Id$
      implicit none
      integer ncor,nocc,nvir,nact,iprt,g_t2,g_nexch
      double precision eorb(*),t2(*),scra(*)
c
      integer iam
c
      integer jlo,jhi,ilo,ihi
      integer a,b,i,j,lhoov,ad1,ad2,lhoo
C
#include "ccsd_len.fh"
#include "global.fh"
#include "tcgmsg.fh"
C
      IAm = GA_NodeID()

      if (iprt.gt.5)print *,' making starting t2 '
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
c
      call ga_sync
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          do i=1,nocc
            ad1=(i-1)*nact+nocc+a
            call ga_get(g_nexch,1,nocc*nact,ad1,ad1,t2,nocc*nact)
            do j=1,i
            do b=1,nvir
              ad1=(b-1)*lhoo+i*(i-1)/2+j
              ad2=(j-1)*nact+nocc+b
              scra(ad1)=t2(ad2)/(eorb(ncor+i)+eorb(ncor+j)
     &                          -eorb(ncor+nocc+a)-eorb(ncor+nocc+b))
            enddo
            enddo
          enddo
          if (iprt.gt.50)then
           print *,'t2 start ',a
           write(*,435)(scra(i),i=1,lhoov)
 435       format(4e15.5)
          endif
          call ga_put(g_t2,1,lhoov,a,a,scra,lhoov)
        endif
      enddo
      call ga_sync
c
      return
      end

