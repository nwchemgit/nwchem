      subroutine selci_sigma(iflham,nci,nroot,hd,ci,si)
*
* $Id: sigma.F,v 1.4 1997-10-31 23:42:27 d3e129 Exp $
*
      implicit real*8 (a-h,o-z)
#include "cbuff.fh"
#include "global.fh"
      dimension hd(nci),ci(nci,nroot),si(nci,nroot)
c
c     compute sigma vectors for each root
c     
      rewind (iflham)
      me = ga_nodeid()
      if (me.eq.0) read (iflham)
      if (me.ne.0) then
        do 10 iroot = 1, nroot
           do 20 i = 1, nci
              si(i, iroot) = 0.0d0
 20        continue
 10     continue
      else
        do 11 iroot = 1, nroot
           do 21 i = 1,nci
              si(i,iroot) = hd(i)*ci(i,iroot)
 21        continue
 11     continue
      endif
c
c
 30   call selci_sread(iflham,rbuff,lenbuf)
      iptr = 0
      ipti = 0
c
 40   i = ibuff(ipti+1)
c      write(6,*) ' sigma: i=',i
      if (i.eq.iendr) goto 30
      if (i.eq.iendf) goto 1000
      i = i + imagic
      nj = ibuff(ipti+2)
      ipti = ipti + 2
c      write(6,*) ' sigma: nj=',nj,' iptr=',iptr,' ipti=',ipti
      if (nroot.eq.2) then
         do 62 jj = 1,nj
            j = ibuff(ipti+jj) + imagic
c            write(6,*) ' sigma: j=',j
            hij = rbuff(iptr+jj)
            si(i,1) = si(i,1) + ci(j,1)*hij
            si(j,1) = si(j,1) + ci(i,1)*hij
            si(i,2) = si(i,2) + ci(j,2)*hij
            si(j,2) = si(j,2) + ci(i,2)*hij
62       continue
      else if (nroot.eq.3) then
         do 63 jj = 1,nj
            j = ibuff(ipti+jj) + imagic
c            write(6,*) ' sigma: j=',j
            hij = rbuff(iptr+jj)
            si(i,1) = si(i,1) + ci(j,1)*hij
            si(j,1) = si(j,1) + ci(i,1)*hij
            si(i,2) = si(i,2) + ci(j,2)*hij
            si(j,2) = si(j,2) + ci(i,2)*hij
            si(i,3) = si(i,3) + ci(j,3)*hij
            si(j,3) = si(j,3) + ci(i,3)*hij
63       continue
      else
         do 50 iroot = 1,nroot
cvd$  nodepchk
            do 60 jj = 1,nj
               j = ibuff(ipti+jj) + imagic
c               write(6,*) ' sigma: j=',j
               hij = rbuff(iptr+jj)
               si(i,iroot) = si(i,iroot) + ci(j,iroot)*hij
               si(j,iroot) = si(j,iroot) + ci(i,iroot)*hij
 60         continue
 50      continue
      endif
      ipti = ipti + nj
      iptr = iptr + nj
      goto 40
c
1000  continue
c
      do 1010 iroot = 1, nroot
         call ga_dgop (99, si(1, iroot), nci, '+')
1010  continue
c
      end
