      subroutine ludcmp(a,n,np,indx,d)
#include "errquit.fh"
C$Id$
      call errquit('ludcmp: this routine needs to be replaced',0,
     &       UNKNOWN_ERR)
c$$$      implicit none
c$$$      integer n, np
c$$$      integer nmax
c$$$      double precision tiny
c$$$      parameter (nmax=100,tiny=1.0d-20)
c$$$      double precision a(np,np),vv(nmax)
c$$$      integer indx(n)
c$$$      double precision d
c$$$c
c$$$      double precision aamax, sum, dum
c$$$      integer i, j, k, imax
c$$$c***********************************************************************
c$$$c
c$$$c  lu decomposition of a matrix from numerical recipies
c$$$c
c$$$c***********************************************************************
c$$$      d=1.0d+00
c$$$      do 12 i=1,n
c$$$        aamax=0.0d+00
c$$$        do 11 j=1,n
c$$$          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
c$$$11      continue
c$$$        if (aamax.eq.0.0d+00) pause 'singular matrix.'
c$$$        vv(i)=1.0d+00/aamax
c$$$12    continue
c$$$      do 19 j=1,n
c$$$        if (j.gt.1) then
c$$$          do 14 i=1,j-1
c$$$            sum=a(i,j)
c$$$            if (i.gt.1)then
c$$$              do 13 k=1,i-1
c$$$                sum=sum-a(i,k)*a(k,j)
c$$$13            continue
c$$$              a(i,j)=sum
c$$$            endif
c$$$14        continue
c$$$        endif
c$$$        aamax=0.0d+00
c$$$        do 16 i=j,n
c$$$          sum=a(i,j)
c$$$          if (j.gt.1)then
c$$$            do 15 k=1,j-1
c$$$              sum=sum-a(i,k)*a(k,j)
c$$$15          continue
c$$$            a(i,j)=sum
c$$$          endif
c$$$          dum=vv(i)*abs(sum)
c$$$          if (dum.ge.aamax) then
c$$$            imax=i
c$$$            aamax=dum
c$$$          endif
c$$$16      continue
c$$$        if (j.ne.imax)then
c$$$          do 17 k=1,n
c$$$            dum=a(imax,k)
c$$$            a(imax,k)=a(j,k)
c$$$            a(j,k)=dum
c$$$17        continue
c$$$          d=-d
c$$$          vv(imax)=vv(j)
c$$$        endif
c$$$        indx(j)=imax
c$$$        if(j.ne.n)then
c$$$          if(a(j,j).eq.0.)a(j,j)=tiny
c$$$          dum=1./a(j,j)
c$$$          do 18 i=j+1,n
c$$$            a(i,j)=a(i,j)*dum
c$$$18        continue
c$$$        endif
c$$$19    continue
c$$$      if(a(n,n).eq.0.0d+00)a(n,n)=tiny
c$$$      return
      end
