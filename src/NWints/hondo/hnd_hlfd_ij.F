c
c $Id: hnd_hlfd_ij.F,v 1.2 1999-05-06 17:25:50 windus Exp $
c
c Taken from HONDO
c
      subroutine hnd_hlfd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1  xyzj,expj,coefj,j_nprim,j_ngen,Lj,kat,didvij,djdvij,zan,xyz,
     2  nder,nint)
c
c  This is a routine that computes additional cross 2nd derivative Helman-
c  Feynman/orbital terms.  It is assumed that the buffers for the integrals
c  have been zeroed before entering this routine.
c
      implicit none
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
      integer i_nprim   ! [input] num. prims on function i
      integer i_ngen    ! [input] num general conts on func. i
      integer Li        ! [input] angular momentum of func. i
      integer j_nprim   ! [input] num. prims on function j
      integer j_ngen    ! [input] num general conts on func. j
      integer Lj        ! [input] angular momentum of func. j
      integer nder      ! [input] should always be 2=2nd der
      integer nint      ! [input] number of base integrals
      integer kat       ! [input] lexical number of an atom
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim) ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen) ! [input] j coeffs
      double precision xyz(3,*) ! [input] all atom positions
      double precision didvij(nint,3,3) ! [output] 2nd der integrals
      double precision djdvij(nint,3,3) ! [output] 2nd der integrals
      double precision zan(*)   ! [input] nuclear charges
c
      character*8 errmsg
      logical out
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      double precision xint, yint, zint, tx, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision ijx, ijy, ijz
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      double precision dumxx, dumyy, dumzz, dumxy, dumxz, dumyz
      double precision dumyx, dumzx, dumzy
      double precision two, four, pi212, aax, aay, aaz, znuc
      double precision uu, u2, u4, ww, w2, w4, tt, xv, yv, zv
      double precision dxv, dyv, dzv, ddxv, ddyv, ddzv
      double precision dxvdi, dyvdi, dzvdi
      double precision dxvdj, dyvdj, dzvdj
      double precision ddxvdi, ddyvdi, ddzvdi
      double precision ddxvdj, ddyvdj, ddzvdj
      integer ni, nj, minprim, maxprim, iroot
      integer i1, i2, lit, mini, maxi, j1, j2, ljt, minj, maxj
      integer ig, jg, ij, i, j, ix, iy, iz, jx, jy, jz, maxrys
      integer litder, ljtder
      dimension pij(225)
      dimension w2(5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension minprim(5), maxprim(5)
      dimension     xv(6,6,5),    yv(6,6,5),    zv(6,6,5)
      dimension    dxv(6,6,5),   dyv(6,6,5),   dzv(6,6,5)
      dimension  dxvdi(5,5,5), dyvdi(5,5,5), dzvdi(5,5,5)
      dimension  dxvdj(5,5,5), dyvdj(5,5,5), dzvdj(5,5,5)
      dimension ddxvdi(5,5,5),ddyvdi(5,5,5),ddzvdi(5,5,5)
      dimension ddxvdj(5,5,5),ddyvdj(5,5,5),ddzvdj(5,5,5)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','- hlfd -'/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
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
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 3,10,20,35/
c
      tol =rln10*itol
c
c     ----- calculate derivatives of -hlf- term -----
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      i1=1
      i2=i_nprim
      lit=Li + 1
      mini=minprim(lit)
      maxi=maxprim(lit)
c
      litder=lit+nder-1
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      j1=1
      j2=j_nprim
      ljt=Lj + 1
      minj=minprim(ljt)
      maxj=maxprim(ljt)
c
      ljtder=ljt+nder-1
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_hlfd_ij2: maxrys is too large!',555)
      endif
c
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=expi(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=zero
      cpi=zero
      cdi=zero
      cfi=zero
      cgi=zero
      if(lit.eq.1) then
        csi=coefi(ig,i_ngen)
      elseif(lit.eq.2) then
        cpi=coefi(ig,i_ngen)
      elseif(lit.eq.3) then
        cdi=coefi(ig,i_ngen)
      elseif(lit.eq.4) then
        cfi=coefi(ig,i_ngen)
      elseif(lit.eq.5) then
        cgi=coefi(ig,i_ngen)
      else
        call errquit('hnd_hlfd_ij2: illegal -lit- ',lit)
      endif
c
c     ----- j primitive -----
c
      do 6000 jg=j1,j2
      aj=expj(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=zero
      cpj=zero
      cdj=zero
      cfj=zero
      cgj=zero
      if(ljt.eq.1) then
        csj=coefj(jg,j_ngen)
      elseif(ljt.eq.2) then
        cpj=coefj(jg,j_ngen)
      elseif(ljt.eq.3) then
        cdj=coefj(jg,j_ngen)
      elseif(ljt.eq.4) then
        cfj=coefj(jg,j_ngen)
      elseif(ljt.eq.5) then
        cgj=coefj(jg,j_ngen)
      else
        call errquit('hnd_hlfd_ij2: illegal -ljt- ',ljt)
      endif
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,220,220,220,
     1       150,220,220,220,220,220,220,220,220,220,
     2       180,220,220,120,220,220,220,220,220,220,
     3       220,220,220,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  150 dum1=cfi*fac
      go to 220
  180 dum1=cgi*fac
      go to 220
  220 continue
c
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,350,350,350,
     1       280,350,350,350,350,350,350,350,350,350,
     2       310,350,350,350,350,350,350,350,350,350,
     3       350,350,350,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  280 dum2=dum1*cfj
      go to 350
  310 dum2=dum1*cgj
      go to 350
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
      cx=xyz(1,kat)
      cy=xyz(2,kat)
      cz=xyz(3,kat)
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
      do 410 j=1,ljtder
      nj=j
      do 410 i=1,litder
      ni=i
      call hnd_sxyz
c     call original hnd_dsyxz instead?
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint
      call hnd_dvxyz     
      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint
  410 continue
c
c probably need to adjust this for the j derivative
c
      call hnd_deriaj( dxvdi(1,1,iroot), dyvdi(1,1,iroot),
     1              dzvdi(1,1,iroot), dxvdj(1,1,iroot),
     1              dyvdj(1,1,iroot), dzvdj(1,1,iroot),
     1              xv(1,1,iroot),    yv(1,1,iroot),    zv(1,1,iroot),
     2          lit,ljt,ai,aj)
      call hnd_deriaj(ddxvdi(1,1,iroot),ddyvdi(1,1,iroot),
     1             ddzvdi(1,1,iroot), ddxvdj(1,1,iroot),
     1             ddyvdj(1,1,iroot), ddzvdj(1,1,iroot),
     1             dxv(1,1,iroot),   dyv(1,1,iroot),   dzv(1,1,iroot),
     2          lit,ljt,ai,aj)
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
      didvij(ij,1,1)=didvij(ij,1,1)+dumxx*(aa1*pi212*pij(ij))
      didvij(ij,1,2)=didvij(ij,1,2)+dumxy*(aa1*pi212*pij(ij))
      didvij(ij,1,3)=didvij(ij,1,3)+dumxz*(aa1*pi212*pij(ij))
      didvij(ij,2,1)=didvij(ij,2,1)+dumyx*(aa1*pi212*pij(ij))
      didvij(ij,2,2)=didvij(ij,2,2)+dumyy*(aa1*pi212*pij(ij))
      didvij(ij,2,3)=didvij(ij,2,3)+dumyz*(aa1*pi212*pij(ij))
      didvij(ij,3,1)=didvij(ij,3,1)+dumzx*(aa1*pi212*pij(ij))
      didvij(ij,3,2)=didvij(ij,3,2)+dumzy*(aa1*pi212*pij(ij))
      didvij(ij,3,3)=didvij(ij,3,3)+dumzz*(aa1*pi212*pij(ij))
c
c      didvij(1,1,ij)=didvij(1,1,ij)+dumxx*(aa1*pi212*pij(ij))
c      didvij(1,2,ij)=didvij(1,2,ij)+dumxy*(aa1*pi212*pij(ij))
c      didvij(1,3,ij)=didvij(1,3,ij)+dumxz*(aa1*pi212*pij(ij))
c      didvij(2,1,ij)=didvij(2,1,ij)+dumyx*(aa1*pi212*pij(ij))
c      didvij(2,2,ij)=didvij(2,2,ij)+dumyy*(aa1*pi212*pij(ij))
c      didvij(2,3,ij)=didvij(2,3,ij)+dumyz*(aa1*pi212*pij(ij))
c      didvij(3,1,ij)=didvij(3,1,ij)+dumzx*(aa1*pi212*pij(ij))
c      didvij(3,2,ij)=didvij(3,2,ij)+dumzy*(aa1*pi212*pij(ij))
c      didvij(3,3,ij)=didvij(3,3,ij)+dumzz*(aa1*pi212*pij(ij))
c
c handle the djdvij part
c
      dumxx=zero
      dumxy=zero
      dumxz=zero
      dumyx=zero
      dumyy=zero
      dumyz=zero
      dumzx=zero
      dumzy=zero
      dumzz=zero
      do 435 iroot=1,nroots
      dumxx=dumxx+w2(iroot)*
     1     (ddxvdj(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxy=dumxy+w2(iroot)*
     1     (dxvdj (ix,jx,iroot)*dyv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxz=dumxz+w2(iroot)*
     1     (dxvdj (ix,jx,iroot)* yv   (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumyx=dumyx+w2(iroot)*
     1     (dxv   (ix,jx,iroot)*dyvdj (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyy=dumyy+w2(iroot)*
     1     ( xv   (ix,jx,iroot)*ddyvdj(iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyz=dumyz+w2(iroot)*
     1     ( xv   (ix,jx,iroot)*dyvdj (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumzx=dumzx+w2(iroot)*
     1     (dxv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdj (iz,jz,iroot))
      dumzy=dumzy+w2(iroot)*
     1     ( xv   (ix,jx,iroot)*dyv   (iy,jy,iroot)*dzvdj (iz,jz,iroot))
      dumzz=dumzz+w2(iroot)*
     1     ( xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*ddzvdj(iz,jz,iroot))
  435 continue
      ij=ij+1
      djdvij(ij,1,1)=djdvij(ij,1,1)+dumxx*(aa1*pi212*pij(ij))
      djdvij(ij,1,2)=djdvij(ij,1,2)+dumxy*(aa1*pi212*pij(ij))
      djdvij(ij,1,3)=djdvij(ij,1,3)+dumxz*(aa1*pi212*pij(ij))
      djdvij(ij,2,1)=djdvij(ij,2,1)+dumyx*(aa1*pi212*pij(ij))
      djdvij(ij,2,2)=djdvij(ij,2,2)+dumyy*(aa1*pi212*pij(ij))
      djdvij(ij,2,3)=djdvij(ij,2,3)+dumyz*(aa1*pi212*pij(ij))
      djdvij(ij,3,1)=djdvij(ij,3,1)+dumzx*(aa1*pi212*pij(ij))
      djdvij(ij,3,2)=djdvij(ij,3,2)+dumzy*(aa1*pi212*pij(ij))
      djdvij(ij,3,3)=djdvij(ij,3,3)+dumzz*(aa1*pi212*pij(ij))
c
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
