c
c $Id: hnd_tvd_ij.F,v 1.1 1999-04-10 00:15:18 windus Exp $
c
c Taken from HONDO
c
      subroutine hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj, xyz,
     2  ditij,djtij,diditij,didjtij,divij,djvij,didivij,didjvij,
     3  nder,nint,zan,nat,doT,doV,iatom,jatom)
      implicit none
#include "stdio.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
      logical doT      ! [input] Need T contributions
      logical doV      ! [input] Need V contributions
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nder     ! [input] 1=1rst der; 2=2nd der
      integer nint     ! [input] number of base integrals
      integer nat      ! [input] number of atoms
      integer iatom    ! [input] lexical number of center i
      integer jatom    ! [input] lexical number of center j
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,nat) ! [input] all atom positions
      double precision ditij ! [output] 1rst derivative wrt i for T
      double precision djtij ! [output] 1rst derivative wrt j for T
      double precision diditij ! [output] 2nd derivative wrt i,i for T
      double precision didjtij ! [output] 2nd derivative wrt i,j for T
      double precision divij ! [output] 1rst derivative wrt i for V
      double precision djvij ! [output] 1rst derivative wrt j for V
      double precision didivij ! [output] 2nd derivative wrt i,i for V
      double precision didjvij ! [output] 2nd derivative wrt i,j for V
      double precision zan(nat) ! [input] nuclear charges
c
      double precision xint, yint, zint, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision ijx, ijy, ijz, xs, ys, zs
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision pi212, xt, yt, zt, dxsdii, dysdii, dzsdii
      double precision dxsdij, dysdij, dzsdij, dxsdi, dysdi, dzsdi
      double precision dxsdj, dysdj, dzsdj, dxtdii, dytdii, dztdii
      double precision dxtdij, dytdij, dztdij, dxtdi, dytdi, dztdi
      double precision dxtdj, dytdj, dztdj, dumx, dumy, dumz
      double precision dumxx, dumyy, dumzz, dumyx, dumzx, dumzy
      double precision dumxy, dumxz, dumyz, aax, aay, aaz
      double precision uu, ww, tt, xv, yv, zv
      double precision dxvdii, dyvdii, dzvdii, dxvdij, dyvdij, dzvdij
      double precision dxvdi, dyvdi, dzvdi, dxvdj, dyvdj, dzvdj
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision tx, aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      integer ni, nj, maxrys, minprim, maxprim, i1, i2, lit
      integer mini, maxi, litder, j1, j2, ljt, minj, maxj 
      integer ljtder, ljtmod, ig, jg, ij, i, j, ix, iy, iz
      integer jx, jy, jz, kat, znuc, iroot, ijtemp
c
      character*8 errmsg
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      dimension ditij(nint,3),djtij(nint,3)
      dimension divij(nint,3),djvij(nint,3)
c     dimension diditij(3,3,225),didjtij(3,3,225)
c     dimension didivij(3,3,225),didjvij(3,3,225)
      dimension diditij(nint,3,3),didjtij(nint,3,3)
      dimension didivij(nint,3,3),didjvij(nint,3,3)
      dimension pij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension     xs(7,9)  ,    ys(7,9)  ,    zs(7,9)
      dimension  dxsdi(5,5)  , dysdi(5,5)  , dzsdi(5,5)
      dimension  dxsdj(5,5)  , dysdj(5,5)  , dzsdj(5,5)
      dimension dxsdii(5,5)  ,dysdii(5,5)  ,dzsdii(5,5)
      dimension dxsdij(5,5)  ,dysdij(5,5)  ,dzsdij(5,5)
      dimension     xt(7,7)  ,    yt(7,7)  ,    zt(7,7)
      dimension  dxtdi(5,5)  , dytdi(5,5)  , dztdi(5,5)
      dimension  dxtdj(5,5)  , dytdj(5,5)  , dztdj(5,5)
      dimension dxtdii(5,5)  ,dytdii(5,5)  ,dztdii(5,5)
      dimension dxtdij(5,5)  ,dytdij(5,5)  ,dztdij(5,5)
      dimension     xv(7,7,5),    yv(7,7,5),    zv(7,7,5)
      dimension  dxvdi(5,5,5), dyvdi(5,5,5), dzvdi(5,5,5)
      dimension  dxvdj(5,5,5), dyvdj(5,5,5), dzvdj(5,5,5)
      dimension dxvdii(5,5,5),dyvdii(5,5,5),dzvdii(5,5,5)
      dimension dxvdij(5,5,5),dyvdij(5,5,5),dzvdij(5,5,5)
      dimension minprim(5), maxprim(5)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','- tvd  -'/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
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
      ljtmod=ljtder+2
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_tvd_ij: maxrys is to large!',555)
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
        call errquit('hnd_tvd_ij: illegal -lit- ', lit)
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
        call errquit('hnd_tvd_ij: illegal -ljt- ', ljt)
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
c     ----- kinetic energy derivatives -----
c
      if (doT) then
      tx= sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljtmod
      nj=j
      do 370 i=1,litder
      ni=i
      call hnd_sxyz
      xs(i,j)=xint*tx
      ys(i,j)=yint*tx
      zs(i,j)=zint*tx
  370 continue
c
      call hnd_ddtxyz(xt,yt,zt,xs,ys,zs,litder,ljtder,aj)
c
      call hnd_dderij2(dxsdii,dysdii,dzsdii,dxsdij,dysdij,
     1            dzsdij,
     1            dxsdi ,dysdi ,dzsdi ,dxsdj ,dysdj ,dzsdj ,
     1            xs,ys,zs,lit,ljt,ai,aj)
      call hnd_dderij2(dxtdii,dytdii,dztdii,dxtdij,dytdij,
     1            dztdij,
     1            dxtdi ,dytdi ,dztdi ,dxtdj ,dytdj ,dztdj ,
     1            xt,yt,zt,lit,ljt,ai,aj)
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
      ij=ij+1
c
      if (nder.eq.1) then
      dumx =dxtdi (ix,jx)* ys   (iy,jy)* zs   (iz,jz)
     1     +dxsdi (ix,jx)* yt   (iy,jy)* zs   (iz,jz)
     2     +dxsdi (ix,jx)* ys   (iy,jy)* zt   (iz,jz)
      dumy = xt   (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
     1     + xs   (ix,jx)*dytdi (iy,jy)* zs   (iz,jz)
     2     + xs   (ix,jx)*dysdi (iy,jy)* zt   (iz,jz)
      dumz = xt   (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
     1     + xs   (ix,jx)* yt   (iy,jy)*dzsdi (iz,jz)
     2     + xs   (ix,jx)* ys   (iy,jy)*dztdi (iz,jz)
      ditij(ij,1)=ditij(ij,1)+ dumx*pij(ij)
      ditij(ij,2)=ditij(ij,2)+ dumy*pij(ij)
      ditij(ij,3)=ditij(ij,3)+ dumz*pij(ij)
c
c do derivative wrt to the second atom
c
      dumx =dxtdj (ix,jx)* ys   (iy,jy)* zs   (iz,jz)
     1     +dxsdj (ix,jx)* yt   (iy,jy)* zs   (iz,jz)
     2     +dxsdj (ix,jx)* ys   (iy,jy)* zt   (iz,jz)
      dumy = xt   (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
     1     + xs   (ix,jx)*dytdj (iy,jy)* zs   (iz,jz)
     2     + xs   (ix,jx)*dysdj (iy,jy)* zt   (iz,jz)
      dumz = xt   (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
     1     + xs   (ix,jx)* yt   (iy,jy)*dzsdj (iz,jz)
     2     + xs   (ix,jx)* ys   (iy,jy)*dztdj (iz,jz)
      djtij(ij,1)=djtij(ij,1)+ dumx*pij(ij)
      djtij(ij,2)=djtij(ij,2)+ dumy*pij(ij)
      djtij(ij,3)=djtij(ij,3)+ dumz*pij(ij)
c
      elseif (nder.eq.2) then
      dumxx=dxtdii(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
     1     +dxsdii(ix,jx)* yt   (iy,jy)* zs   (iz,jz)
     2     +dxsdii(ix,jx)* ys   (iy,jy)* zt   (iz,jz)
      dumyy= xt   (ix,jx)*dysdii(iy,jy)* zs   (iz,jz)
     1     + xs   (ix,jx)*dytdii(iy,jy)* zs   (iz,jz)
     2     + xs   (ix,jx)*dysdii(iy,jy)* zt   (iz,jz)
      dumzz= xt   (ix,jx)* ys   (iy,jy)*dzsdii(iz,jz)
     1     + xs   (ix,jx)* yt   (iy,jy)*dzsdii(iz,jz)
     2     + xs   (ix,jx)* ys   (iy,jy)*dztdii(iz,jz)
      dumyx=dxtdi (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
     1     +dxsdi (ix,jx)*dytdi (iy,jy)* zs   (iz,jz)
     2     +dxsdi (ix,jx)*dysdi (iy,jy)* zt   (iz,jz)
      dumzx=dxtdi (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
     1     +dxsdi (ix,jx)* yt   (iy,jy)*dzsdi (iz,jz)
     2     +dxsdi (ix,jx)* ys   (iy,jy)*dztdi (iz,jz)
      dumzy= xt   (ix,jx)*dysdi (iy,jy)*dzsdi (iz,jz)
     1     + xs   (ix,jx)*dytdi (iy,jy)*dzsdi (iz,jz)
     2     + xs   (ix,jx)*dysdi (iy,jy)*dztdi (iz,jz)
c     diditij(1,1,ij)=diditij(1,1,ij)+ dumxx*pij(ij)
c     diditij(1,2,ij)=diditij(1,2,ij)+ dumyx*pij(ij)
c     diditij(1,3,ij)=diditij(1,3,ij)+ dumzx*pij(ij)
c     diditij(2,1,ij)=diditij(2,1,ij)+ dumyx*pij(ij)
c     diditij(2,2,ij)=diditij(2,2,ij)+ dumyy*pij(ij)
c     diditij(2,3,ij)=diditij(2,3,ij)+ dumzy*pij(ij)
c     diditij(3,1,ij)=diditij(3,1,ij)+ dumzx*pij(ij)
c     diditij(3,2,ij)=diditij(3,2,ij)+ dumzy*pij(ij)
c     diditij(3,3,ij)=diditij(3,3,ij)+ dumzz*pij(ij)
      diditij(ij,1,1)=diditij(ij,1,1)+ dumxx*pij(ij)
      diditij(ij,1,2)=diditij(ij,1,2)+ dumyx*pij(ij)
      diditij(ij,1,3)=diditij(ij,1,3)+ dumzx*pij(ij)
      diditij(ij,2,1)=diditij(ij,2,1)+ dumyx*pij(ij)
      diditij(ij,2,2)=diditij(ij,2,2)+ dumyy*pij(ij)
      diditij(ij,2,3)=diditij(ij,2,3)+ dumzy*pij(ij)
      diditij(ij,3,1)=diditij(ij,3,1)+ dumzx*pij(ij)
      diditij(ij,3,2)=diditij(ij,3,2)+ dumzy*pij(ij)
      diditij(ij,3,3)=diditij(ij,3,3)+ dumzz*pij(ij)
      dumxx=dxtdij(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
     1     +dxsdij(ix,jx)* yt   (iy,jy)* zs   (iz,jz)
     2     +dxsdij(ix,jx)* ys   (iy,jy)* zt   (iz,jz)
      dumxy=dxtdi (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
     1     +dxsdi (ix,jx)*dytdj (iy,jy)* zs   (iz,jz)
     2     +dxsdi (ix,jx)*dysdj (iy,jy)* zt   (iz,jz)
      dumxz=dxtdi (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
     1     +dxsdi (ix,jx)* yt   (iy,jy)*dzsdj (iz,jz)
     2     +dxsdi (ix,jx)* ys   (iy,jy)*dztdj (iz,jz)
      dumyx=dxtdj (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
     1     +dxsdj (ix,jx)*dytdi (iy,jy)* zs   (iz,jz)
     2     +dxsdj (ix,jx)*dysdi (iy,jy)* zt   (iz,jz)
      dumyy= xt   (ix,jx)*dysdij(iy,jy)* zs   (iz,jz)
     1     + xs   (ix,jx)*dytdij(iy,jy)* zs   (iz,jz)
     2     + xs   (ix,jx)*dysdij(iy,jy)* zt   (iz,jz)
      dumyz= xt   (ix,jx)*dysdi (iy,jy)*dzsdj (iz,jz)
     1     + xs   (ix,jx)*dytdi (iy,jy)*dzsdj (iz,jz)
     2     + xs   (ix,jx)*dysdi (iy,jy)*dztdj (iz,jz)
      dumzx=dxtdj (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
     1     +dxsdj (ix,jx)* yt   (iy,jy)*dzsdi (iz,jz)
     2     +dxsdj (ix,jx)* ys   (iy,jy)*dztdi (iz,jz)
      dumzy= xt   (ix,jx)*dysdj (iy,jy)*dzsdi (iz,jz)
     1     + xs   (ix,jx)*dytdj (iy,jy)*dzsdi (iz,jz)
     2     + xs   (ix,jx)*dysdj (iy,jy)*dztdi (iz,jz)
      dumzz= xt   (ix,jx)* ys   (iy,jy)*dzsdij(iz,jz)
     1     + xs   (ix,jx)* yt   (iy,jy)*dzsdij(iz,jz)
     2     + xs   (ix,jx)* ys   (iy,jy)*dztdij(iz,jz)
c     didjtij(1,1,ij)=didjtij(1,1,ij)+ dumxx*pij(ij)
c     didjtij(1,2,ij)=didjtij(1,2,ij)+ dumxy*pij(ij)
c     didjtij(1,3,ij)=didjtij(1,3,ij)+ dumxz*pij(ij)
c     didjtij(2,1,ij)=didjtij(2,1,ij)+ dumyx*pij(ij)
c     didjtij(2,2,ij)=didjtij(2,2,ij)+ dumyy*pij(ij)
c     didjtij(2,3,ij)=didjtij(2,3,ij)+ dumyz*pij(ij)
c     didjtij(3,1,ij)=didjtij(3,1,ij)+ dumzx*pij(ij)
c     didjtij(3,2,ij)=didjtij(3,2,ij)+ dumzy*pij(ij)
c     didjtij(3,3,ij)=didjtij(3,3,ij)+ dumzz*pij(ij)
      didjtij(ij,1,1)=didjtij(ij,1,1)+ dumxx*pij(ij)
      didjtij(ij,1,2)=didjtij(ij,1,2)+ dumxy*pij(ij)
      didjtij(ij,1,3)=didjtij(ij,1,3)+ dumxz*pij(ij)
      didjtij(ij,2,1)=didjtij(ij,2,1)+ dumyx*pij(ij)
      didjtij(ij,2,2)=didjtij(ij,2,2)+ dumyy*pij(ij)
      didjtij(ij,2,3)=didjtij(ij,2,3)+ dumyz*pij(ij)
      didjtij(ij,3,1)=didjtij(ij,3,1)+ dumzx*pij(ij)
      didjtij(ij,3,2)=didjtij(ij,3,2)+ dumzy*pij(ij)
      didjtij(ij,3,3)=didjtij(ij,3,3)+ dumzz*pij(ij)

      else ! if nder.gt.2
        write(luout,9996) nder
        call errquit('hnd_tvd_ij: derivative order is too high',555)
      endif

  380 continue
  390 continue
      endif ! doT
c
c     ----- nuclear attraction derivatives -----
c
      if (doV) then
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 kat=1,nat
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
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint*ww
  410 continue
c
      call hnd_dderij2(dxvdii(1,1,iroot),dyvdii(1,1,iroot),
     1            dzvdii(1,1,iroot),
     1            dxvdij(1,1,iroot),dyvdij(1,1,iroot),dzvdij(1,1,iroot),
     2             dxvdi(1,1,iroot), dyvdi(1,1,iroot), dzvdi(1,1,iroot),
     3             dxvdj(1,1,iroot), dyvdj(1,1,iroot), dzvdj(1,1,iroot),
     4                xv(1,1,iroot),    yv(1,1,iroot),    zv(1,1,iroot),
     5            lit,ljt,ai,aj)
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
c
      if (nder.eq.1) then
c
c first do the derivative wrt to the first atom
c
      ijtemp = ij
      dumx=zero
      dumy=zero
      dumz=zero
      do 430 iroot=1,nroots
      dumx=dumx+ 
     1  dxvdi(ix,jx,iroot)*  yv  (iy,jy,iroot)*  zv  (iz,jz,iroot)
      dumy=dumy+ 
     1   xv  (ix,jx,iroot)* dyvdi(iy,jy,iroot)*  zv  (iz,jz,iroot)
      dumz=dumz+ 
     1   xv  (ix,jx,iroot)*  yv  (iy,jy,iroot)* dzvdi(iz,jz,iroot)
  430 continue
      dumx=dumx*aa1*pi212
      dumy=dumy*aa1*pi212
      dumz=dumz*aa1*pi212
      ij=ij+1
      divij(ij,1)=divij(ij,1)+dumx*pij(ij)
      divij(ij,2)=divij(ij,2)+dumy*pij(ij)
      divij(ij,3)=divij(ij,3)+dumz*pij(ij)
c
c first do the derivative wrt to the second atom
c
      ij = ijtemp
      dumx=zero
      dumy=zero
      dumz=zero
      do 431 iroot=1,nroots
      dumx=dumx+ 
     1  dxvdj(ix,jx,iroot)*  yv  (iy,jy,iroot)*  zv  (iz,jz,iroot)
      dumy=dumy+ 
     1   xv  (ix,jx,iroot)* dyvdj(iy,jy,iroot)*  zv  (iz,jz,iroot)
      dumz=dumz+ 
     1   xv  (ix,jx,iroot)*  yv  (iy,jy,iroot)* dzvdj(iz,jz,iroot)
  431 continue
      dumx=dumx*aa1*pi212
      dumy=dumy*aa1*pi212
      dumz=dumz*aa1*pi212
      ij=ij+1
      djvij(ij,1)=djvij(ij,1)+dumx*pij(ij)
      djvij(ij,2)=djvij(ij,2)+dumy*pij(ij)
      djvij(ij,3)=djvij(ij,3)+dumz*pij(ij)
c
      elseif (nder.eq.2) then
      dumxx=zero
      dumyy=zero
      dumzz=zero
      dumyx=zero
      dumzx=zero
      dumzy=zero
      do 432 iroot=1,nroots
      dumxx=dumxx
     1     +dxvdii(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumyy=dumyy
     1     + xv   (ix,jx,iroot)*dyvdii(iy,jy,iroot)* zv   (iz,jz,iroot)
      dumzz=dumzz
     1     + xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdii(iz,jz,iroot)
      dumyx=dumyx
     1     +dxvdi (ix,jx,iroot)*dyvdi (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumzx=dumzx
     1     +dxvdi (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdi (iz,jz,iroot)
      dumzy=dumzy
     1     + xv   (ix,jx,iroot)*dyvdi (iy,jy,iroot)*dzvdi (iz,jz,iroot)
  432 continue
      dumxx=dumxx*aa1*pi212
      dumyy=dumyy*aa1*pi212
      dumzz=dumzz*aa1*pi212
      dumyx=dumyx*aa1*pi212
      dumzx=dumzx*aa1*pi212
      dumzy=dumzy*aa1*pi212
c     didivij(1,1,ij)=didivij(1,1,ij)+dumxx*pij(ij)
c     didivij(1,2,ij)=didivij(1,2,ij)+dumyx*pij(ij)
c     didivij(1,3,ij)=didivij(1,3,ij)+dumzx*pij(ij)
c     didivij(2,1,ij)=didivij(2,1,ij)+dumyx*pij(ij)
c     didivij(2,2,ij)=didivij(2,2,ij)+dumyy*pij(ij)
c     didivij(2,3,ij)=didivij(2,3,ij)+dumzy*pij(ij)
c     didivij(3,1,ij)=didivij(3,1,ij)+dumzx*pij(ij)
c     didivij(3,2,ij)=didivij(3,2,ij)+dumzy*pij(ij)
c     didivij(3,3,ij)=didivij(3,3,ij)+dumzz*pij(ij)
      didivij(ij,1,1)=didivij(ij,1,1)+dumxx*pij(ij)
      didivij(ij,1,2)=didivij(ij,1,2)+dumyx*pij(ij)
      didivij(ij,1,3)=didivij(ij,1,3)+dumzx*pij(ij)
      didivij(ij,2,1)=didivij(ij,2,1)+dumyx*pij(ij)
      didivij(ij,2,2)=didivij(ij,2,2)+dumyy*pij(ij)
      didivij(ij,2,3)=didivij(ij,2,3)+dumzy*pij(ij)
      didivij(ij,3,1)=didivij(ij,3,1)+dumzx*pij(ij)
      didivij(ij,3,2)=didivij(ij,3,2)+dumzy*pij(ij)
      didivij(ij,3,3)=didivij(ij,3,3)+dumzz*pij(ij)
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
      do 434 iroot=1,nroots
      dumxx=dumxx
     1     +dxvdij(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumxy=dumxy
     1     +dxvdi (ix,jx,iroot)*dyvdj (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumxz=dumxz
     1     +dxvdi (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdj (iz,jz,iroot)
      dumyx=dumyx
     1     +dxvdj (ix,jx,iroot)*dyvdi (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumyy=dumyy
     1     + xv   (ix,jx,iroot)*dyvdij(iy,jy,iroot)* zv   (iz,jz,iroot)
      dumyz=dumyz
     1     + xv   (ix,jx,iroot)*dyvdi (iy,jy,iroot)*dzvdj (iz,jz,iroot)
      dumzx=dumzx
     1     +dxvdj (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdi (iz,jz,iroot)
      dumzy=dumzy
     1     + xv   (ix,jx,iroot)*dyvdj (iy,jy,iroot)*dzvdi (iz,jz,iroot)
      dumzz=dumzz
     1     + xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdij(iz,jz,iroot)
  434 continue
      dumxx=dumxx*aa1*pi212
      dumxy=dumxy*aa1*pi212
      dumxz=dumxz*aa1*pi212
      dumyx=dumyx*aa1*pi212
      dumyy=dumyy*aa1*pi212
      dumyz=dumyz*aa1*pi212
      dumzx=dumzx*aa1*pi212
      dumzy=dumzy*aa1*pi212
      dumzz=dumzz*aa1*pi212
c     didjvij(1,1,ij)=didjvij(1,1,ij)+dumxx*pij(ij) 
c     didjvij(1,2,ij)=didjvij(1,2,ij)+dumxy*pij(ij) 
c     didjvij(1,3,ij)=didjvij(1,3,ij)+dumxz*pij(ij) 
c     didjvij(2,1,ij)=didjvij(2,1,ij)+dumyx*pij(ij) 
c     didjvij(2,2,ij)=didjvij(2,2,ij)+dumyy*pij(ij) 
c     didjvij(2,3,ij)=didjvij(2,3,ij)+dumyz*pij(ij) 
c     didjvij(3,1,ij)=didjvij(3,1,ij)+dumzx*pij(ij) 
c     didjvij(3,2,ij)=didjvij(3,2,ij)+dumzy*pij(ij) 
c     didjvij(3,3,ij)=didjvij(3,3,ij)+dumzz*pij(ij) 
      didjvij(ij,1,1)=didjvij(ij,1,1)+dumxx*pij(ij) 
      didjvij(ij,1,2)=didjvij(ij,1,2)+dumxy*pij(ij) 
      didjvij(ij,1,3)=didjvij(ij,1,3)+dumxz*pij(ij) 
      didjvij(ij,2,1)=didjvij(ij,2,1)+dumyx*pij(ij) 
      didjvij(ij,2,2)=didjvij(ij,2,2)+dumyy*pij(ij) 
      didjvij(ij,2,3)=didjvij(ij,2,3)+dumyz*pij(ij) 
      didjvij(ij,3,1)=didjvij(ij,3,1)+dumzx*pij(ij) 
      didjvij(ij,3,2)=didjvij(ij,3,2)+dumzy*pij(ij) 
      didjvij(ij,3,3)=didjvij(ij,3,3)+dumzz*pij(ij) 

      else ! if nder.gt.2
        write(luout,9996) nder
        call errquit('hnd_tvd_ij: derivative order is too high',555)
      endif

  440 continue
  450 continue
  500 continue
      endif ! doV
c
 6000 continue
 7000 continue
c
      return
 9996 format(' in -tvd- , derivatives of order ',i3,
     1       ' are not implemented')
 9997 format(' in -tvd- , the rys quadrature is not implememented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots= ',3i3)
      end
