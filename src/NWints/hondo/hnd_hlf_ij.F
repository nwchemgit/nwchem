c
c $Id: hnd_hlf_ij.F,v 1.2 1999-05-06 17:25:50 windus Exp $
c
c Taken from HONDO
c
      subroutine hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj, kat, dvij, ddvij,
     2  zan, xyz, nder, nint)
c
c  This is the routine that actually computes the 1rst and 2nd
c  derivatives for the Helman-Feynman term.  It is assumed that the buffers
c  for the integrals have been zeroed before entering this routine.
c
      implicit none
#include "stdio.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nder     ! [input] 1=1rst der; 2=2nd der
      integer nint     ! [input] number of base integrals
      integer kat      ! [input] lexical number of an atom
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,*)  ! [input] all atom positions
c     double precision dvij(3,nint) ! [output] 1rst. der integrals
      double precision dvij(nint,3) ! [output] 1rst. der integrals
      double precision ddvij(nint,3,3) ! [output] 2nd der integrals
      double precision zan(*)   ! [input] nuclear charges
c
      character*8 errmsg
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      double precision xint, yint, zint, tx, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision ijx, ijy, ijz
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      double precision dumx, dumy, dumz, dumxx, dumyy, dumzz
      double precision dumxy, dumxz, dumyz
      double precision two, four, pi212, aax, aay, aaz, znuc
      double precision uu, u2, u4, ww, w2, w4, tt, xv, yv, zv
      double precision dxv, dyv, dzv, ddxv, ddyv, ddzv
      integer ni, nj, minprim, maxprim, iroot
      integer i1, i2, lit, mini, maxi, j1, j2, ljt, minj, maxj
      integer ig, jg, ij, i, j, ix, iy, iz, jx, jy, jz, maxrys
      dimension pij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension minprim(5), maxprim(5)
      dimension   xv(5,5,5),  yv(5,5,5),  zv(5,5,5)
      dimension  dxv(5,5,5), dyv(5,5,5), dzv(5,5,5)
      dimension ddxv(5,5,5),ddyv(5,5,5),ddzv(5,5,5)
      dimension w2(5),w4(5)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-hlfspd-'/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data four   /4.0d+00/
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
      data maxprim / 1, 4,10,20,35/
c
      tol =rln10*itol
c
c     ----- calculate -helfey- term -----
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      i1=1
      i2=i_nprim
      lit = Li + 1
      mini=minprim(lit)
      maxi=maxprim(lit)
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      j1=1
      j2=j_nprim
      ljt = Lj + 1
      minj=minprim(ljt)
      maxj=maxprim(ljt)
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call hnd_hnderr(3,errmsg)
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
        call errquit('hnd_hlf_ij2: illegal -lit- ', lit)
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
        call errquit('hnd_stv_int: illegal -ljt- ', ljt)
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
     2       180,220,220,220,220,220,220,220,220,220,
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
  350 continue
c
      ij=ij+1
      pij(ij)=dum2
  360 continue
c
c     ----- hellmann-feynman term -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
c
c     ----- kat -----
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
      u2=uu
      u4=uu*uu
      ww=w(iroot)*znuc
      w2(iroot)=ww*u2*two
      w4(iroot)=ww*u4*four
      tt=one/(aa+uu)
      tx= sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      do 410 j=1,ljt
      nj=j
      do 410 i=1,lit
      ni=i

      call hnd_sxyz

      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint

      call hnd_dvxyz

      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint

c
c if nder = 2 only?
c
      call hnd_ddvxyz

      ddxv(i,j,iroot)=xint
      ddyv(i,j,iroot)=yint
      ddzv(i,j,iroot)=zint
  410 continue
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
      dumx=zero
      dumy=zero
      dumz=zero
      dumxx=zero
      dumyy=zero
      dumzz=zero
      dumxy=zero
      dumxz=zero
      dumyz=zero
      do 430 iroot=1,nroots
      dumx =dumx+
     1    dxv(ix,jx,iroot)*  yv(iy,jy,iroot)*  zv(iz,jz,iroot)*w2(iroot)
      dumy =dumy+
     1     xv(ix,jx,iroot)* dyv(iy,jy,iroot)*  zv(iz,jz,iroot)*w2(iroot)
      dumz =dumz+
     1     xv(ix,jx,iroot)*  yv(iy,jy,iroot)* dzv(iz,jz,iroot)*w2(iroot)
      dum  =
     1     xv(ix,jx,iroot)*  yv(iy,jy,iroot)*  zv(iz,jz,iroot)*w2(iroot)
      dumxx=dumxx-dum+
     1   ddxv(ix,jx,iroot)*  yv(iy,jy,iroot)*  zv(iz,jz,iroot)*w4(iroot)
      dumyy=dumyy-dum+
     1     xv(ix,jx,iroot)*ddyv(iy,jy,iroot)*  zv(iz,jz,iroot)*w4(iroot)
      dumzz=dumzz-dum+
     1     xv(ix,jx,iroot)*  yv(iy,jy,iroot)*ddzv(iz,jz,iroot)*w4(iroot)
      dumxy=dumxy+
     1    dxv(ix,jx,iroot)* dyv(iy,jy,iroot)*  zv(iz,jz,iroot)*w4(iroot)
      dumxz=dumxz+
     1    dxv(ix,jx,iroot)*  yv(iy,jy,iroot)* dzv(iz,jz,iroot)*w4(iroot)
      dumyz=dumyz+
     1     xv(ix,jx,iroot)* dyv(iy,jy,iroot)* dzv(iz,jz,iroot)*w4(iroot)
  430 continue
      ij=ij+1
      if (nder.eq.1) then
        dumx =dumx *(pi212*aa1*pij(ij))
        dumy =dumy *(pi212*aa1*pij(ij))
        dumz =dumz *(pi212*aa1*pij(ij))
        dvij(ij,1)=dvij(ij,1)+dumx
        dvij(ij,2)=dvij(ij,2)+dumy
        dvij(ij,3)=dvij(ij,3)+dumz
      elseif (nder.eq.2) then
        dumxx=dumxx*(pi212*aa1*pij(ij))
        dumyy=dumyy*(pi212*aa1*pij(ij))
        dumzz=dumzz*(pi212*aa1*pij(ij))
        dumxy=dumxy*(pi212*aa1*pij(ij))
        dumxz=dumxz*(pi212*aa1*pij(ij))
        dumyz=dumyz*(pi212*aa1*pij(ij))
        ddvij(ij,1,1)=ddvij(ij,1,1)+dumxx    
        ddvij(ij,1,2)=ddvij(ij,1,2)+dumxy
        ddvij(ij,1,3)=ddvij(ij,1,3)+dumxz
        ddvij(ij,2,1)=ddvij(ij,2,1)+dumxy
        ddvij(ij,2,2)=ddvij(ij,2,2)+dumyy
        ddvij(ij,2,3)=ddvij(ij,2,3)+dumyz
        ddvij(ij,3,1)=ddvij(ij,3,1)+dumxz
        ddvij(ij,3,2)=ddvij(ij,3,2)+dumyz
        ddvij(ij,3,3)=ddvij(ij,3,3)+dumzz
      endif
      dum=(pi212*aa1)*pij(ij)
c
  440 continue
  450 continue
c
 6000 continue
 7000 continue
c
      return
 9997 format(' in -hlf- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots = ',3i3)
      end
