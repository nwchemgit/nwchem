      subroutine hnd_giaxyz(gijkl,hijkl,gnkl,hnkl,fnkl,gnm,
     1                  ng,nmax,mmax,nimax,njmax,nkmax,nlmax,
     2                  dij,dkl)
c
c $Id: hnd_giaxyz.F,v 1.2 2004-07-22 16:31:41 edo Exp $
c
      implicit none 
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
c
      integer ng,nmax,mmax,nimax,njmax,nkmax,nlmax
      double precision gijkl(ng,3,nlmax*nkmax,njmax,nimax)
      double precision hijkl(ng,3,nlmax*nkmax*njmax,nimax)
      double precision  gnkl(ng,nlmax,nkmax,nmax)
      double precision  hnkl(ng,nlmax*nkmax,nmax)
      double precision  fnkl(ng,nlmax*nkmax*nmax)
      double precision   gnm(ng,nmax,mmax)
      double precision   dij(ng), dkl(ng)
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
   10 gnkl(ig,nl,nk,n)=gnm(ig,n,nl)
      if(nk.eq.nkmax) go to 30
      maxi=mmax-nk
      do 20  m=1,maxi
      do 20  n=1,nmax
      do 20 ig=1,ng
   20 gnm(ig,n,m)=dkl(ig)*gnm(ig,n,m)+gnm(ig,n,m+1)
   30 continue
c
      go to 100
   40 continue
c
      do 70 nl=1,nlmax
      do 50 nk=1,nkmax
      do 50  n=1,nmax
      do 50 ig=1,ng
   50 gnkl(ig,nl,nk,n)=gnm(ig,n,nk)
      if(nl.eq.nlmax) go to 70
      maxi=mmax-nl
      do 60  m=1,maxi
      do 60  n=1,nmax
      do 60 ig=1,ng
   60 gnm(ig,n,m)=dkl(ig)*gnm(ig,n,m)+gnm(ig,n,m+1)
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
  110 hijkl(ig,1,lkj,ni)=fnkl(ig,lkj)
      if(ni.eq.nimax) go to 130
      maxi=nmax-ni
      do 120  n=1,maxi
      do 120 nk=1,nkmax
      do 120 nl=1,nlmax
      do 120 ig=1,ng
  120 gnkl(ig,nl,nk,n)=dij(ig)*gnkl(ig,nl,nk,n)+gnkl(ig,nl,nk,n+1)
  130 continue
c
      return
  140 continue
c
      do 170 nj=1,njmax
      do 150 ni=1,nimax
      do 150 lk=1,nlmax*nkmax
      do 150 ig=1,ng
  150 gijkl(ig,1,lk,nj,ni)=hnkl(ig,lk,ni)
      if(nj.eq.njmax) go to 170
      maxi=nmax-nj
      do 160  n=1,maxi
      do 160 nk=1,nkmax
      do 160 nl=1,nlmax
      do 160 ig=1,ng
  160 gnkl(ig,nl,nk,n)=dij(ig)*gnkl(ig,nl,nk,n)+gnkl(ig,nl,nk,n+1)
  170 continue
c
      return
      end
