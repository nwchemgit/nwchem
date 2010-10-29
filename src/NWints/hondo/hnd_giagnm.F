      subroutine hnd_giagnm(gnm,ng,nmax,mmax,b00,b01,b10,c00,d00,f00)
c
c $Id$
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
c
      dimension gnm(ng,3,nmax,mmax)
      dimension c00(ng,3),d00(ng,3),f00(ng)
      dimension b00(ng,1),b01(ng,1),b10(ng,1)
c
c     ----- g(0,0) -----
c
      do 10 ig=1,ng
      gnm(ig,1,1,1)=f00(ig)
      gnm(ig,2,1,1)=1.0d0 
      gnm(ig,3,1,1)=1.0d0
   10 continue
      if(nmaxs.and.mmaxs) return
      if(nmaxs) go to 30
c
c     ----- g(1,0) = c00 * g(0,0) -----
c
      do 20 ig=1,ng
      gnm(ig,1,2,1)=c00(ig,1)*gnm(ig,1,1,1)
      gnm(ig,2,2,1)=c00(ig,2)*gnm(ig,2,1,1)
      gnm(ig,3,2,1)=c00(ig,3)*gnm(ig,3,1,1)
   20 continue
c
   30 continue
      if(mmaxs) go to 60
c
c     ----- g(0,1) = d00 * g(0,0) -----
c
      do 40 ig=1,ng
      gnm(ig,1,1,2)=d00(ig,1)*gnm(ig,1,1,1)
      gnm(ig,2,1,2)=d00(ig,2)*gnm(ig,2,1,1)
      gnm(ig,3,1,2)=d00(ig,3)*gnm(ig,3,1,1)
   40 continue
      if(nmaxs) go to 60
c
c     ----- g(1,1) = b00 * g(0,0) + d00 * g(1,0) -----
c
      do 50 ig=1,ng
      gnm(ig,1,2,2)=b00(ig,1)*gnm(ig,1,1,1)+d00(ig,1)*gnm(ig,1,2,1)
      gnm(ig,2,2,2)=b00(ig,1)*gnm(ig,2,1,1)+d00(ig,2)*gnm(ig,2,2,1)
      gnm(ig,3,2,2)=b00(ig,1)*gnm(ig,3,1,1)+d00(ig,3)*gnm(ig,3,2,1)
   50 continue
c
   60 continue
      max=max0(nmax-1,mmax-1)
      do 70 m=2,max
      do 70 ig=1,ng
   70 b00(ig,m)=b00(ig,m-1)+b00(ig,1)
c
      if(nmaxp) go to 120
c
c     ----- g(n+1,0) = n * b10 * g(n-1,0) + c00 * g(n,0) -----
c
      do 80 n=2,nmax-1
      do 80 ig=1,ng
      b10(ig,n)=b10(ig,n-1)+b10(ig,1)
   80 continue
      do 90 n=2,nmax-1
      do 90 ig=1,ng
      gnm(ig,1,n+1,1)=b10(ig,n-1)*gnm(ig,1,n-1,1)
     &                +c00(ig,1)*gnm(ig,1,n,1)
      gnm(ig,2,n+1,1)=b10(ig,n-1)*gnm(ig,2,n-1,1)
     &                +c00(ig,2)*gnm(ig,2,n,1)
      gnm(ig,3,n+1,1)=b10(ig,n-1)*gnm(ig,3,n-1,1)
     &                +c00(ig,3)*gnm(ig,3,n,1)
   90 continue
      if(mmaxs) go to 110
c
c     ----- g(n,1) = n * b00 * g(n-1,0) + d00 * g(n,0) -----
c
      do 100 n=2,nmax-1
      do 100 ig=1,ng
      gnm(ig,1,n+1,2)=b00(ig,n)*gnm(ig,1,n,1)+d00(ig,1)*gnm(ig,1,n+1,1)
      gnm(ig,2,n+1,2)=b00(ig,n)*gnm(ig,2,n,1)+d00(ig,2)*gnm(ig,2,n+1,1)
      gnm(ig,3,n+1,2)=b00(ig,n)*gnm(ig,3,n,1)+d00(ig,3)*gnm(ig,3,n+1,1)
  100 continue
c
  110 continue
c
  120 continue
      if(mmaxp) go to 170
c
c     ----- g(0,m+1) = m * b01 * g(0,m-1) + d00 * g(o,m) -----
c
      do 130 m=2,mmax-1
      do 130 ig=1,ng
      b01(ig,m)=b01(ig,m-1)+b01(ig,1)
  130 continue
      do 140 m=2,mmax-1
      do 140 ig=1,ng
      gnm(ig,1,1,m+1)=b01(ig,m-1)*gnm(ig,1,1,m-1)
     &                +d00(ig,1)*gnm(ig,1,1,m)
      gnm(ig,2,1,m+1)=b01(ig,m-1)*gnm(ig,2,1,m-1)
     &                +d00(ig,2)*gnm(ig,2,1,m)
      gnm(ig,3,1,m+1)=b01(ig,m-1)*gnm(ig,3,1,m-1)
     &                +d00(ig,3)*gnm(ig,3,1,m)
  140 continue
      if(nmaxs) go to 160
c
c     ----- g(1,m) = m * b00 * g(0,m-1) + c00 * g(0,m) -----
c
      do 150 m=2,mmax-1
      do 150 ig=1,ng
      gnm(ig,1,2,m+1)=b00(ig,m)*gnm(ig,1,1,m)+c00(ig,1)*gnm(ig,1,1,m+1)
      gnm(ig,2,2,m+1)=b00(ig,m)*gnm(ig,2,1,m)+c00(ig,2)*gnm(ig,2,1,m+1)
      gnm(ig,3,2,m+1)=b00(ig,m)*gnm(ig,3,1,m)+c00(ig,3)*gnm(ig,3,1,m+1)
  150 continue
c
  160 continue
  170 if(nmaxp.or.mmaxp) return
c
c     ----- g(n+1,m) = n * b10 * g(n-1,m  ) -----
c                    +     c00 * g(n  ,m  )
c                    + m * b00 * g(n  ,m-1)
c
      do 180 m=2,mmax-1
      do 180 n=2,nmax-1
      do 180 ig=1,ng
      gnm(ig,1,n+1,m+1)=b10(ig,n-1)*gnm(ig,1,n-1,m+1)+
     1                  c00(ig,1  )*gnm(ig,1,n  ,m+1)+
     2                  b00(ig,m  )*gnm(ig,1,n  ,m  )
      gnm(ig,2,n+1,m+1)=b10(ig,n-1)*gnm(ig,2,n-1,m+1)+
     1                  c00(ig,2  )*gnm(ig,2,n  ,m+1)+
     2                  b00(ig,m  )*gnm(ig,2,n  ,m  )
      gnm(ig,3,n+1,m+1)=b10(ig,n-1)*gnm(ig,3,n-1,m+1)+
     1                  c00(ig,3  )*gnm(ig,3,n  ,m+1)+
     2                  b00(ig,m  )*gnm(ig,3,n  ,m  )
  180 continue
c
      return
      end
