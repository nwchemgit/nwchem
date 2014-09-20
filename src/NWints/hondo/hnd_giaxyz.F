      subroutine hnd_giaxyz(gijkl,hijkl,gnkl,hnkl,fnkl,gnm,
     1                  ng,nmax,mmax,nimax,njmax,nkmax,nlmax)
c
c $Id$
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
      if(.not.expndk) then
c
         do nk=1,nkmax
            do nl=1,nlmax
               do n=1,nmax
                  do ig=1,ng
                     gnkl(ig,1,nl,nk,n)=gnm(ig,1,n,nl)
                     gnkl(ig,2,nl,nk,n)=gnm(ig,2,n,nl)
                     gnkl(ig,3,nl,nk,n)=gnm(ig,3,n,nl)
                  enddo
               enddo
            enddo
            if(nk.lt.nkmax) then
               maxi=mmax-nk
CDEC$ NOINLINE
               call hnd_giaxyz3(maxi,nmax,ng,
     D              dxkl,dykl,dzkl,gnm)
            endif
         enddo
c     
      else
c
         do nl=1,nlmax
            do nk=1,nkmax
               do n=1,nmax
                  do ig=1,ng
                     gnkl(ig,1,nl,nk,n)=gnm(ig,1,n,nk)
                     gnkl(ig,2,nl,nk,n)=gnm(ig,2,n,nk)
                     gnkl(ig,3,nl,nk,n)=gnm(ig,3,n,nk)
                  enddo
               enddo
            enddo
            if(nl.lt.nlmax) then
               maxi=mmax-nl
CDEC$ NOINLINE
               call hnd_giaxyz3(maxi,nmax,ng,
     D              dxkl,dykl,dzkl,gnm)
            endif
         enddo
c     
      endif
c     
c     ----- g(i,j,k,l) -----
c
      if(.not.expndi) then
c     
         do ni=1,nimax
            do lkj=1,nlmax*nkmax*njmax
               do ig=1,ng
                  hijkl(ig,1,1,lkj,ni)=fnkl(ig,1,lkj)
                  hijkl(ig,2,1,lkj,ni)=fnkl(ig,2,lkj)
                  hijkl(ig,3,1,lkj,ni)=fnkl(ig,3,lkj)
               enddo
            enddo
            if(ni.lt.nimax) then
               maxi=nmax-ni
CDEC$ NOINLINE
               call hnd_giaxyz4(maxi,nkmax,nlmax,ng,
     D              dxij,dyij,dzij,gnkl)
            endif
         enddo
c     
      else
c     
         do nj=1,njmax
            do ni=1,nimax
               do lk=1,nlmax*nkmax
                  do  ig=1,ng
                     gijkl(ig,1,1,lk,nj,ni)=hnkl(ig,1,lk,ni)
                     gijkl(ig,2,1,lk,nj,ni)=hnkl(ig,2,lk,ni)
                     gijkl(ig,3,1,lk,nj,ni)=hnkl(ig,3,lk,ni)
                  enddo
               enddo
            enddo
            if(nj.lt.njmax) then
               maxi=nmax-nj
CDEC$ NOINLINE
               call hnd_giaxyz4(maxi,nkmax,nlmax,ng,
     D              dxij,dyij,dzij,gnkl)
            endif
         enddo
      endif
c     
      return
      end
      subroutine hnd_giaxyz4(maxi,nkmax,nlmax,ng,
     D     dxij,dyij,dzij,gnkl)
      implicit none
      integer maxi,nkmax,nlmax,ng
      double precision dxij,dyij,dzij
      double precision gnkl(ng,3,nlmax,nkmax,*)
c
      integer n,nk,nl,ig
      do n=1,maxi
         do nk=1,nkmax
            do nl=1,nlmax
               do ig=1,ng
                  gnkl(ig,1,nl,nk,n)=dxij*
     *                 gnkl(ig,1,nl,nk,n)+gnkl(ig,1,nl,nk,n+1)
                  gnkl(ig,2,nl,nk,n)=dyij*
     *                 gnkl(ig,2,nl,nk,n)+gnkl(ig,2,nl,nk,n+1)
                  gnkl(ig,3,nl,nk,n)=dzij*
     *                 gnkl(ig,3,nl,nk,n)+gnkl(ig,3,nl,nk,n+1)
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine hnd_giaxyz3(maxi,nmax,ng,
     D     dxkl,dykl,dzkl,gnm)
      implicit none
      integer maxi,nmax,ng
      double precision dxkl,dykl,dzkl
      double precision gnm(ng,3,nmax,*)
c
      integer m,n,ig
      do m=1,maxi
         do n=1,nmax
            do ig=1,ng
               gnm(ig,1,n,m)=dxkl*gnm(ig,1,n,m)+gnm(ig,1,n,m+1)
               gnm(ig,2,n,m)=dykl*gnm(ig,2,n,m)+gnm(ig,2,n,m+1)
               gnm(ig,3,n,m)=dzkl*gnm(ig,3,n,m)+gnm(ig,3,n,m+1)
            enddo
         enddo
      enddo
      return
      end
