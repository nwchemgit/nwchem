      subroutine hnd_hlfd_ij(ii,jj,iat,kat,didvij)
      implicit double precision (a-h,o-z)
#include "hnd_rys.fh"
#include "hnd_tol.fh"
      parameter (mxatom=500)
      parameter (mxprim=2048)
      parameter (mxshel=512)
      character*8 errmsg
      logical iandj
      logical norm
      logical out
      common/hnd_iofile/ir,iw,ip
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_basnum/num
      common/hnd_nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      dimension pij(225)
      dimension didvij(3,3,225)
      dimension w2(5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension     xv(6,5,5),    yv(6,5,5),    zv(6,5,5)
      dimension    dxv(6,5,5),   dyv(6,5,5),   dzv(6,5,5)
      dimension  dxvdi(5,5,5), dyvdi(5,5,5), dzvdi(5,5,5)
      dimension ddxvdi(5,5,5),ddyvdi(5,5,5),ddzvdi(5,5,5)
      dimension dnam(3)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','- hlfd -'/
      data dnam   /3he'x,3he'y,3he'z/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
c
c     ----- order of components in shells in hondo -----
c
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      out =.false.
      norm=normf.ne.1.or.normp.ne.1
      tol =rln10*itol
c
c     ----- calculate derivatives of -hlf- term -----
c
      nder=2
c
c     ----- ishell -----
c
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
      litder=lit+nder-1
      iat   =i
c
c     ----- jshell -----
c
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      jat   =j
c
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(iw,9997) maxrys,lit,ljt,nroots
         call hnd_hnderr(3,errmsg)
      endif
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            do kxyz=1,3
               do ixyz=1,3
                  didvij(ixyz,kxyz,ij)=zero
               enddo
            enddo
         enddo
      enddo
c
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do 6000 jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
      pij(ij)=dum2
  360 continue
c
c     ----- nuclear attraction derivatives -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
c                      
      znuc=-zan(kat)
      cx=c(1,kat)
      cy=c(2,kat)
      cz=c(3,kat)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call hnd_rt123
      if(nroots.eq.4) call hnd_root4
      if(nroots.eq.5) call hnd_root5
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      ww=w(iroot)*znuc
      w2(iroot)=ww*uu*two
      tt=one/(aa+uu)
      tx= sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
c
      do 410 j=1,ljt
      nj=j
      do 410 i=1,litder
      ni=i
      call hnd_dsxyz
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint
      call hnd_dvxyz
      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint
  410 continue
c
      call hnd_deri( dxvdi(1,1,iroot), dyvdi(1,1,iroot),
     1              dzvdi(1,1,iroot),
     1              xv(1,1,iroot),    yv(1,1,iroot),    zv(1,1,iroot),
     2          lit,ljt,ai)
      call hnd_deri(ddxvdi(1,1,iroot),ddyvdi(1,1,iroot),
     1             ddzvdi(1,1,iroot),
     1             dxv(1,1,iroot),   dyv(1,1,iroot),   dzv(1,1,iroot),
     2          lit,ljt,ai)
c
  420 continue
c
      ij=0
      do 450 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      do 440 j=minj,maxj
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      dumxx=zero
      dumxy=zero
      dumxz=zero
      dumyx=zero
      dumyy=zero
      dumyz=zero
      dumzx=zero
      dumzy=zero
      dumzz=zero
      do 430 iroot=1,nroots
      dumxx=dumxx+w2(iroot)*
     1     (ddxvdi(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxy=dumxy+w2(iroot)*
     1     (dxvdi (ix,jx,iroot)*dyv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxz=dumxz+w2(iroot)*
     1     (dxvdi (ix,jx,iroot)* yv   (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumyx=dumyx+w2(iroot)*
     1     (dxv   (ix,jx,iroot)*dyvdi (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyy=dumyy+w2(iroot)*
     1     ( xv   (ix,jx,iroot)*ddyvdi(iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyz=dumyz+w2(iroot)*
     1     ( xv   (ix,jx,iroot)*dyvdi (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumzx=dumzx+w2(iroot)*
     1     (dxv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdi (iz,jz,iroot))
      dumzy=dumzy+w2(iroot)*
     1     ( xv   (ix,jx,iroot)*dyv   (iy,jy,iroot)*dzvdi (iz,jz,iroot))
      dumzz=dumzz+w2(iroot)*
     1     ( xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*ddzvdi(iz,jz,iroot))
  430 continue
      ij=ij+1
      didvij(1,1,ij)=didvij(1,1,ij)+dumxx*(aa1*pi212*pij(ij))
      didvij(1,2,ij)=didvij(1,2,ij)+dumxy*(aa1*pi212*pij(ij))
      didvij(1,3,ij)=didvij(1,3,ij)+dumxz*(aa1*pi212*pij(ij))
      didvij(2,1,ij)=didvij(2,1,ij)+dumyx*(aa1*pi212*pij(ij))
      didvij(2,2,ij)=didvij(2,2,ij)+dumyy*(aa1*pi212*pij(ij))
      didvij(2,3,ij)=didvij(2,3,ij)+dumyz*(aa1*pi212*pij(ij))
      didvij(3,1,ij)=didvij(3,1,ij)+dumzx*(aa1*pi212*pij(ij))
      didvij(3,2,ij)=didvij(3,2,ij)+dumzy*(aa1*pi212*pij(ij))
      didvij(3,3,ij)=didvij(3,3,ij)+dumzz*(aa1*pi212*pij(ij))
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
      return
 9997 format(' in -hlfd- , the rys quadrature is not implememented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots= ',3i3)
      end
