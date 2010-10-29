      subroutine ccsd_mktask(basis,nsh,tklst)
C     $Id$
      implicit none
#include "errquit.fh"
      integer basis,nsh,tklst(nsh*(nsh+1)/2,2)
#include "bas.fh"
#include "rtdb.fh"
c
      integer itsk,ish,jsh,ibflo,ibfhi,jbflo,jbfhi
c
      itsk=0
      do ish=1,nsh
       if (.not. bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &      call errquit('mktask: bas_cn2bfr',ish, BASIS_ERR)
       do jsh=1,ish
        if (.not. bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &       call errquit('mktask: bas_cn2bfr',jsh, BASIS_ERR)
        itsk=itsk+1
        tklst(itsk,1)=(ibfhi-ibflo+1)*(jbfhi-jbflo+1)
        tklst(itsk,2)=(ish-1)*nsh+jsh
       enddo
      enddo
      call isortx(tklst(1,1),-1,itsk,tklst(1,2))
c
      return
      end
      subroutine isortx(ix,inc,n,indx)
      Implicit NONE
      Integer Inc, N
      Integer IX(*),indx(*)
c -- sort array ix in ascending/decending order and set up index
c    array indx, such that indx(i)=original position of i
c      do 100 i=1,n
c      indx(i)=i
c 100   continue
C
      Integer I, L, IXrra, Irra, J, IncX, Ir
C
      if (inc.eq.0)then
        write(6,105)inc,n
 105       format(' error isortx: inc, n ',2i10)
        stop 20
      else if(inc.lt.0)then
        incx=-inc
        do 110 i=0,n-1
        ix(i*incx+1)=-ix(i*incx+1)
 110       continue
      else
        incx=inc
      endif
c
      l=n/2+1
      ir=n
 10     continue
        if(l.gt.1) then
          l=l-1
          ixrra=ix((l-1)*incx+1)
          irra=indx(l)
        else
          ixrra=ix((ir-1)*incx+1)
          irra=indx(ir)
          ix((ir-1)*incx+1)=ix(1)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1) then
            ix(1)=ixrra
            indx(1)=irra
            goto 999
          end if
        end if
        i=l
        j=2*l
 20         if(j.le.ir) then
          if(j.lt.ir) then
            if(ix((j-1)*incx+1).lt.ix(j*incx+1)) j=j+1
          end if
          if(ixrra.lt.ix((j-1)*incx+1)) then
            ix((i-1)*incx+1)=ix((j-1)*incx+1)
            indx(i)=indx(j)
            i=j
            j=2*j
          else
            j=ir+1
          end if
        goto 20
        end if
        ix((i-1)*incx+1)=ixrra
        indx(i)=irra
      goto 10
c
 999   continue
      if(inc.lt.0)then
        do 210 i=0,n-1
        ix(i*incx+1)=-ix(i*incx+1)
 210       continue
      endif
      return
      end
