c
c $Id: hnd_sd_ij.F,v 1.1 1999-04-10 00:15:16 windus Exp $
c
c Taken from HONDO
c
      subroutine hnd_sd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj,
     1  disij,djsij,didisij,didjsij,nder,nint)
c
c  This is the routine that actually computes the 1rst and 2nd
c  derivatives of the overlap.  It is assumed that the buffers
c  for the integrals have been zeroed before entering this routine
c  (disij,djsij,didisij,didjsij).
c
      implicit none
#include "hnd_tol.fh"
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nder     ! [input] 1=1rst der; 2=2nd der
      integer nint     ! [input] number of base integrals
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision disij ! [output] 1rst derivative wrt i
      double precision djsij ! [output] 1rst derivative wrt j
      double precision didisij ! [output] 2nd derivative wrt i,i
      double precision didjsij ! [output] 2nd derivative wrt i,j

      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      double precision xint, yint, zint, t, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision ijx, ijy, ijz, xs, ys, zs
      double precision dxsdi, dysdi, dzsdi, dxsdj, dysdj, dzsdj
      double precision dxsdii, dysdii, dzsdii, dxsdij, dysdij, dzsdij
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      double precision dumx, dumy, dumz, dumxx, dumyy, dumzz
      double precision dumxy, dumxz, dumyz, dumyx, dumzx, dumzy
      integer ni, nj, minprim, maxprim
      integer i1, i2, lit, mini, maxi, litder, j1, j2, ljt, minj, maxj
      integer ljtder, ig, jg, ij, i, j, ix, iy, iz, jx, jy, jz
      dimension didisij(nint,3,3),didjsij(nint,3,3)
      dimension disij(nint,3),djsij(nint,3)
      dimension pij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension  xs   (7,6), ys   (7,6), zs   (7,6)
      dimension dxsdi (5,5),dysdi (5,5),dzsdi (5,5)
      dimension dxsdj (5,5),dysdj (5,5),dzsdj (5,5)
      dimension dxsdii(5,5),dysdii(5,5),dzsdii(5,5)
      dimension dxsdij(5,5),dysdij(5,5),dzsdij(5,5)
      dimension minprim(5), maxprim(5)
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
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
      litder=lit+nder
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
      ljtder=ljt+nder
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
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
        call errquit('hnd_stv_int: illegal -lit- ', lit)
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
c     ----- overlap derivatives -----
c
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljtder
      nj=j
      do 370 i=1,litder
      ni=i
      call hnd_sxyz
      xs(i,j)=xint*t
      ys(i,j)=yint*t
      zs(i,j)=zint*t
  370 continue
c
      call hnd_dderij2(dxsdii,dysdii,dzsdii,dxsdij,dysdij,dzsdij,
     1            dxsdi ,dysdi ,dzsdi ,dxsdj ,dysdj ,dzsdj ,
     2            xs,ys,zs,lit,ljt,ai,aj)
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      do 380 j=minj,maxj
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
c
      ij=ij+1
c
      if (nder.eq.1) then
c
c do derivative wrt the first atom
c
        dumx =dxsdi (ix,jx)* ys   (iy,jy)* zs   (iz,jz)
        dumy = xs   (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
        dumz = xs   (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
c
c switching the ij and xyz breaks unit stride, but it works
c for the moment - tlw
c
        disij(ij,1)=disij(ij,1)+ dumx*pij(ij)
        disij(ij,2)=disij(ij,2)+ dumy*pij(ij)
        disij(ij,3)=disij(ij,3)+ dumz*pij(ij)
c
c do derivative wrt the second atom
c
        dumx =dxsdj (ix,jx)* ys   (iy,jy)* zs   (iz,jz)
        dumy = xs   (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
        dumz = xs   (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
        djsij(ij,1)=djsij(ij,1)+ dumx*pij(ij)
        djsij(ij,2)=djsij(ij,2)+ dumy*pij(ij)
        djsij(ij,3)=djsij(ij,3)+ dumz*pij(ij)
      endif
c
      if (nder.eq.2) then
      dumxx=dxsdii(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
      dumyy= xs   (ix,jx)*dysdii(iy,jy)* zs   (iz,jz)
      dumzz= xs   (ix,jx)* ys   (iy,jy)*dzsdii(iz,jz)
      dumxy=dxsdi (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
      dumxz=dxsdi (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
      dumyz= xs   (ix,jx)*dysdi (iy,jy)*dzsdi (iz,jz)
c     didisij(1,1,ij)=didisij(1,1,ij)+ dumxx*pij(ij) 
c     didisij(1,2,ij)=didisij(1,2,ij)+ dumxy*pij(ij) 
c     didisij(1,3,ij)=didisij(1,3,ij)+ dumxz*pij(ij) 
c     didisij(2,1,ij)=didisij(2,1,ij)+ dumxy*pij(ij) 
c     didisij(2,2,ij)=didisij(2,2,ij)+ dumyy*pij(ij) 
c     didisij(2,3,ij)=didisij(2,3,ij)+ dumyz*pij(ij) 
c     didisij(3,1,ij)=didisij(3,1,ij)+ dumxz*pij(ij) 
c     didisij(3,2,ij)=didisij(3,2,ij)+ dumyz*pij(ij) 
c     didisij(3,3,ij)=didisij(3,3,ij)+ dumzz*pij(ij) 
      didisij(ij,1,1)=didisij(ij,1,1)+ dumxx*pij(ij) 
      didisij(ij,1,2)=didisij(ij,1,2)+ dumxy*pij(ij) 
      didisij(ij,1,3)=didisij(ij,1,3)+ dumxz*pij(ij) 
      didisij(ij,2,1)=didisij(ij,2,1)+ dumxy*pij(ij) 
      didisij(ij,2,2)=didisij(ij,2,2)+ dumyy*pij(ij) 
      didisij(ij,2,3)=didisij(ij,2,3)+ dumyz*pij(ij) 
      didisij(ij,3,1)=didisij(ij,3,1)+ dumxz*pij(ij) 
      didisij(ij,3,2)=didisij(ij,3,2)+ dumyz*pij(ij) 
      didisij(ij,3,3)=didisij(ij,3,3)+ dumzz*pij(ij) 
c
      dumxx=dxsdij(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
      dumxy=dxsdi (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
      dumxz=dxsdi (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
      dumyx=dxsdj (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
      dumyy= xs   (ix,jx)*dysdij(iy,jy)* zs   (iz,jz)
      dumyz= xs   (ix,jx)*dysdi (iy,jy)*dzsdj (iz,jz)
      dumzx=dxsdj (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
      dumzy= xs   (ix,jx)*dysdj (iy,jy)*dzsdi (iz,jz)
      dumzz= xs   (ix,jx)* ys   (iy,jy)*dzsdij(iz,jz)
c     didjsij(1,1,ij)=didjsij(1,1,ij)+ dumxx*pij(ij) 
c     didjsij(1,2,ij)=didjsij(1,2,ij)+ dumxy*pij(ij) 
c     didjsij(1,3,ij)=didjsij(1,3,ij)+ dumxz*pij(ij) 
c     didjsij(2,1,ij)=didjsij(2,1,ij)+ dumyx*pij(ij) 
c     didjsij(2,2,ij)=didjsij(2,2,ij)+ dumyy*pij(ij) 
c     didjsij(2,3,ij)=didjsij(2,3,ij)+ dumyz*pij(ij) 
c     didjsij(3,1,ij)=didjsij(3,1,ij)+ dumzx*pij(ij) 
c     didjsij(3,2,ij)=didjsij(3,2,ij)+ dumzy*pij(ij) 
c     didjsij(3,3,ij)=didjsij(3,3,ij)+ dumzz*pij(ij) 
      didjsij(ij,1,1)=didjsij(ij,1,1)+ dumxx*pij(ij) 
      didjsij(ij,1,2)=didjsij(ij,1,2)+ dumxy*pij(ij) 
      didjsij(ij,1,3)=didjsij(ij,1,3)+ dumxz*pij(ij) 
      didjsij(ij,2,1)=didjsij(ij,2,1)+ dumyx*pij(ij) 
      didjsij(ij,2,2)=didjsij(ij,2,2)+ dumyy*pij(ij) 
      didjsij(ij,2,3)=didjsij(ij,2,3)+ dumyz*pij(ij) 
      didjsij(ij,3,1)=didjsij(ij,3,1)+ dumzx*pij(ij) 
      didjsij(ij,3,2)=didjsij(ij,3,2)+ dumzy*pij(ij) 
      didjsij(ij,3,3)=didjsij(ij,3,3)+ dumzz*pij(ij) 
      endif

  380 continue
  390 continue
c
 6000 continue
 7000 continue
c
      return
      end
