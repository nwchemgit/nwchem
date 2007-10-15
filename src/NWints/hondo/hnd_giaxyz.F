      subroutine hnd_giaxyz(gijkl,hijkl,gnkl,hnkl,fnkl,gnm,
     1                  ng,nmax,mmax,nimax,njmax,nkmax,nlmax)
c
c $Id: hnd_giaxyz.F,v 1.3 2007-10-15 23:37:24 bert Exp $
c
      implicit none 
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
c
      integer ng,nmax,mmax,nimax,njmax,nkmax,nlmax
      double precision gijkl(ng,3,3,nlmax*nkmax,njmax,nimax)
      double precision hijkl(ng,3,3,nlmax*nkmax*njmax,nimax)
      double precision  gnkl(ng,3,nlmax,nkmax,nmax)
      double precision  hnkl(ng,3,nlmax*nkmax,nmax)
      double precision  fnkl(ng,3,nlmax*nkmax*nmax)
      double precision   gnm(ng,3,nmax,mmax)
c
      integer nk,nl,n,ig,m,ni,lkj,nj,lk,maxi
c
c     ----- g(n,k,l) -----
c
      if(expndk) go to 40
c
      do 30 nk=1,nkmax
      do 10 nl=1,nlmax
      do 10  n=1,nmax
      do 10 ig=1,ng
      gnkl(ig,1,nl,nk,n)=gnm(ig,1,n,nl)
      gnkl(ig,2,nl,nk,n)=gnm(ig,2,n,nl)
      gnkl(ig,3,nl,nk,n)=gnm(ig,3,n,nl)
   10 continue
      if(nk.eq.nkmax) go to 30
      maxi=mmax-nk
      do 20  m=1,maxi
      do 20  n=1,nmax
      do 20 ig=1,ng
      gnm(ig,1,n,m)=dxkl*gnm(ig,1,n,m)+gnm(ig,1,n,m+1)
      gnm(ig,2,n,m)=dykl*gnm(ig,2,n,m)+gnm(ig,2,n,m+1)
      gnm(ig,3,n,m)=dzkl*gnm(ig,3,n,m)+gnm(ig,3,n,m+1)
   20 continue
   30 continue
c
      go to 100
   40 continue
c
      do 70 nl=1,nlmax
      do 50 nk=1,nkmax
      do 50  n=1,nmax
      do 50 ig=1,ng
      gnkl(ig,1,nl,nk,n)=gnm(ig,1,n,nk)
      gnkl(ig,2,nl,nk,n)=gnm(ig,2,n,nk)
      gnkl(ig,3,nl,nk,n)=gnm(ig,3,n,nk)
   50 continue
      if(nl.eq.nlmax) go to 70
      maxi=mmax-nl
      do 60  m=1,maxi
      do 60  n=1,nmax
      do 60 ig=1,ng
      gnm(ig,1,n,m)=dxkl*gnm(ig,1,n,m)+gnm(ig,1,n,m+1)
      gnm(ig,2,n,m)=dykl*gnm(ig,2,n,m)+gnm(ig,2,n,m+1)
      gnm(ig,3,n,m)=dzkl*gnm(ig,3,n,m)+gnm(ig,3,n,m+1)
   60 continue
   70 continue
c
  100 continue
c
c     ----- g(i,j,k,l) -----
c
      if(expndi) go to 140
c
      do 130 ni=1,nimax
      do 110 lkj=1,nlmax*nkmax*njmax
      do 110 ig=1,ng
      hijkl(ig,1,1,lkj,ni)=fnkl(ig,1,lkj)
      hijkl(ig,2,1,lkj,ni)=fnkl(ig,2,lkj)
      hijkl(ig,3,1,lkj,ni)=fnkl(ig,3,lkj)
  110 continue
      if(ni.eq.nimax) go to 130
      maxi=nmax-ni
      do 120  n=1,maxi
      do 120 nk=1,nkmax
      do 120 nl=1,nlmax
      do 120 ig=1,ng
      gnkl(ig,1,nl,nk,n)=dxij*gnkl(ig,1,nl,nk,n)+gnkl(ig,1,nl,nk,n+1)
      gnkl(ig,2,nl,nk,n)=dyij*gnkl(ig,2,nl,nk,n)+gnkl(ig,2,nl,nk,n+1)
      gnkl(ig,3,nl,nk,n)=dzij*gnkl(ig,3,nl,nk,n)+gnkl(ig,3,nl,nk,n+1)
  120 continue
  130 continue
c
      return
  140 continue
c
      do 170 nj=1,njmax
      do 150 ni=1,nimax
      do 150 lk=1,nlmax*nkmax
      do 150 ig=1,ng
      gijkl(ig,1,1,lk,nj,ni)=hnkl(ig,1,lk,ni)
      gijkl(ig,2,1,lk,nj,ni)=hnkl(ig,2,lk,ni)
      gijkl(ig,3,1,lk,nj,ni)=hnkl(ig,3,lk,ni)
  150 continue
      if(nj.eq.njmax) go to 170
      maxi=nmax-nj
      do 160  n=1,maxi
      do 160 nk=1,nkmax
      do 160 nl=1,nlmax
      do 160 ig=1,ng
      gnkl(ig,1,nl,nk,n)=dxij*gnkl(ig,1,nl,nk,n)+gnkl(ig,1,nl,nk,n+1)
      gnkl(ig,2,nl,nk,n)=dyij*gnkl(ig,2,nl,nk,n)+gnkl(ig,2,nl,nk,n+1)
      gnkl(ig,3,nl,nk,n)=dzij*gnkl(ig,3,nl,nk,n)+gnkl(ig,3,nl,nk,n+1)
  160 continue
  170 continue
c
      return
      end
