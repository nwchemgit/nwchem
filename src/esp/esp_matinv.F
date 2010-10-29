      logical function esp_matinv(a,n,is,iad1,iad2,d,mdm)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
c
c     inversion of square matrix a by means of the gauss-jordan
c     algorithm
c
      integer mdm,n
      real*8 a(mdm,mdm),d(mdm)
      integer is(2,mdm),iad1(mdm),iad2(mdm)
c
      integer i,j,k,l,m,ima
      real*8 b,e
c
      esp_matinv=.true.
      i=0
      k=0
c
      do 1 l=1,n
        is(1,l)=0
        is(2,l)=0
    1 continue
c
      do 2 ima=1,n
        b=zero
        do 3 l=1,n
          do 4 m=1,n
            if(is(1,l).ne.1.and.is(2,m).ne.1) then
              e=dabs(a(l,m))
              if(e.ge.b) then
                i=l
                k=m
              endif
              b=dmax1(b,e)
            endif
    4     continue
    3   continue
        is(1,i)=1
        is(2,k)=1
        iad1(k)=i
        iad2(i)=k
        b=a(i,k)
c
        if(dabs(b).lt.tiny) then
          b = 1
          a(i,k)=0
          esp_matinv=.false.
        else
          a(i,k)=one/b
        endif
        do 5 l=1,n
          if(l.ne.k) a(i,l)=-a(i,l)/b
    5   continue
c
        do 6 l=1,n
          do 7 m=1,n
            if(l.ne.i.and.m.ne.k) a(l,m)=a(l,m)+a(l,k)*a(i,m)
    7     continue
    6   continue
c
        do 8 l=1,n
          if(l.ne.i) a(l,k)=a(l,k)/b
    8   continue
    2 continue
c
      do 9 l=1,n
        do 10 j=1,n
          k=iad1(j)
          d(j)=a(k,l)
   10   continue
        do 11 j=1,n
          a(j,l)=d(j)
   11   continue
    9 continue
c
      do 12 l=1,n
        do 13 j=1,n
          k=iad2(j)
          d(j)=a(l,k)
   13   continue
        do 14 j=1,n
          a(l,j)=d(j)
   14   continue
   12 continue
c
      return
      end
