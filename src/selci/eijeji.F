C> \ingroup selci
C> @{
C>
      subroutine selci_eijeji(e,odonly,i,j,ns,indbar,iocc,w1,w2,work)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension e(*),indbar(*),iocc(*),work(*),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1)
      logical odonly
c
      parameter (maxtmp = 2001)
      dimension tmp(maxtmp)
c
c     e(u,v) = <Iu|EijEji|Iv>, where u,v label the spin functions
c     and the orbital occupation of I is specified by ns, iocc
c     and indbar. Only for i.ne.j., and both i and j singly
c     occupied.
c
c     if (odonly) only form the actual diagonal ... THE REST IS JUNK
c     (this to optimize select)
c
c     e(u,v) = sum(q) <Iu|Eij|Sq><Sq|Eji|Iv>
c
c     work must be at least nf(ns)*nf(ns-2)
c     e is at least nf(ns)*nf(ns)
c
      if (iocc(i).ne.1 .or. iocc(j).ne.1)
     $     call errquit('eijeji: wrong occupancy',iocc(i)*10+iocc(j),
     &       INPUT_ERR)
c
      call selci_eij(work,i,j,ns,indbar,iocc,w1,w2)
c
      numf = nf(ns)
      numf2 = nf(ns-2)
c
      if (odonly) then
         if (numf .gt. maxtmp) call errquit
     $        ('eijeji: hard dim failed',maxtmp, CALC_ERR)
         do k = 1, numf
            tmp(k) = 0.0d0
         enddo
         mm = 0
         do m = 1, numf2
            do k = 1, numf
               wkm = work(k+mm)
               tmp(k) = tmp(k) + wkm*wkm
            enddo
            mm = mm + numf
         enddo
         kk = 1
         numfp1 = numf + 1
         do k = 1, numf
            e(kk) = tmp(k)
            kk = kk + numfp1
         enddo
      else
         call selci_axbt(work,numf,work,numf,e,numf,numf,numf2,numf)
      endif
c
      end
C>
C> @}
