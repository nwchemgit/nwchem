c
c $Id$
c
c Modified from HONDO
c
      subroutine hnd_sd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj,
     1  disij,djsij,didisij,didjsij,djdjsij,nder,nint,scr,lscr)
c
c This is a wrapper routine, setting up scratch blocks used in actual
c integral routine
c
      implicit none
c
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nder     ! [input] 1=1rst der; 2=2nd der
      integer nint     ! [input] number of base integrals
      integer lscr     ! [input] length of scratch space
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision disij(nint,3)! [output] 1rst derivative wrt i
      double precision djsij(nint,3) ! [output] 1rst derivative wrt j
      double precision didisij(nint,3,3) ! [output] 2nd derivative wrt i,i
      double precision didjsij(nint,3,3) ! [output] 2nd derivative wrt i,j
      double precision djdjsij(nint,3,3) ! [output] 2nd derivative wrt j,j
      double precision scr(lscr) ! [input] scratch buffer
c
# include "hnd_pointers.fh"
c
c Use scratch for temps in integral calculation
c Scratch needs are dependent on nder:
c nder=1 (3*(Li+3)*(Lj+2) + 6*(Li+1)*(Lj+1)) + 9 zero length blocks
c nder=2 (3*(Li+3)*(Lj+2) +15*(Li+1)*(Lj+1))
c
c The pointers are defined in hnd_pointers.fh
c
      call hnd_sd_ij21(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1     xyzj,expj,coefj, j_nprim, j_ngen, Lj,
     2     disij,djsij,didisij,didjsij,djdjsij,nder,nint,
     3     scr(tvp(1)) ,scr(tvp(2)) ,scr(tvp(3)) ,scr(tvp(7)),
     4     scr(tvp(8)) ,scr(tvp(9)) ,scr(tvp(10)),scr(tvp(11)),
     5     scr(tvp(12)),scr(tvp(28)),scr(tvp(29)),scr(tvp(30)),
     6     scr(tvp(31)),scr(tvp(32)),scr(tvp(33)),scr(tvp(34)),
     7     scr(tvp(35)),scr(tvp(36)))
c
      return
      end
c
      subroutine hnd_sd_ij21(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj,
     2  disij,djsij,didisij,didjsij,djdjsij,nder,nint,
     3  xs,ys,zs,dxsdi,dysdi,dzsdi,dxsdj,dysdj,dzsdj,dxsdii,dysdii,
     4  dzsdii,dxsdjj,dysdjj,dzsdjj,dxsdij,dysdij,dzsdij)
c
c  This is the routine that actually computes the 1rst and 2nd
c  derivatives of the overlap.  It is assumed that the buffers
c  for the integrals have been zeroed before entering this routine
c  (disij,djsij,didisij,didjsij,djdjsij).
c
      implicit none
#include "hnd_tol.fh"
#include "errquit.fh"
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
      double precision djdjsij ! [output] 2nd derivative wrt j,j

      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      double precision xint, yint, zint, t, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision ijx, ijy, ijz, xs, ys, zs
      double precision dxsdi, dysdi, dzsdi, dxsdj, dysdj, dzsdj
      double precision dxsdii, dysdii, dzsdii 
      double precision dxsdjj, dysdjj, dzsdjj 
      double precision dxsdij, dysdij, dzsdij
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      double precision dumx, dumy, dumz, dumxx, dumyy, dumzz
      double precision dumxy, dumxz, dumyz, dumyx, dumzx, dumzy
      integer ni, nj, Nxyz
      integer lit, maxi, litder, ljt, maxj
      integer ljtder, ig, jg, ij, i, j, ix, iy, iz, jx, jy, jz
      dimension didisij(nint,3,3),didjsij(nint,3,3),djdjsij(nint,3,3)
      dimension disij(nint,3),djsij(nint,3)
      dimension xs(Li+3,Lj+2), ys(Li+3,Lj+2), zs(Li+3,Lj+2)
      dimension dxsdi (Li+1,Lj+1),dysdi (Li+1,Lj+1),dzsdi (Li+1,Lj+1)
      dimension dxsdj (Li+1,Lj+1),dysdj (Li+1,Lj+1),dzsdj (Li+1,Lj+1)
      dimension dxsdii(Li+1,Lj+1),dysdii(Li+1,Lj+1),dzsdii(Li+1,Lj+1)
      dimension dxsdjj(Li+1,Lj+1),dysdjj(Li+1,Lj+1),dzsdjj(Li+1,Lj+1)
      dimension dxsdij(Li+1,Lj+1),dysdij(Li+1,Lj+1),dzsdij(Li+1,Lj+1)
      dimension Nxyz(3)
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
c
      tol =rln10*itol
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      lit = Li + 1
      maxi=lit*(lit+1)/2
c
      litder=lit+nder
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      ljt = Lj + 1
      maxj=ljt*(ljt+1)/2
c
      ljtder=ljt+nder
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- i primitive -----
c
      do 7000 ig=1,i_nprim
      ai=expi(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=coefi(ig,i_ngen)
c
c     ----- j primitive -----
c
      do 6000 jg=1,j_nprim
      aj=expj(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=coefj(jg,j_ngen)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      pij=csi*csj*fac
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
     1            dxsdjj,dysdjj,dzsdjj,
     2            dxsdi ,dysdi ,dzsdi ,dxsdj ,dysdj ,dzsdj ,
     3            xs,ys,zs,lit,ljt,ai,aj,nder)
c
      ij=0
      do 390 i=1,maxi
      call getNxyz(Li,i,Nxyz)
      ix = Nxyz(1) + 1
      iy = Nxyz(2) + 1
      iz = Nxyz(3) + 1
      do 380 j=1,maxj
      call getNxyz(Lj,j,Nxyz)
      jx = Nxyz(1) + 1
      jy = Nxyz(2) + 1
      jz = Nxyz(3) + 1
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
        disij(ij,1)=disij(ij,1)+ dumx*pij
        disij(ij,2)=disij(ij,2)+ dumy*pij
        disij(ij,3)=disij(ij,3)+ dumz*pij
c
c do derivative wrt the second atom
c
        dumx =dxsdj (ix,jx)* ys   (iy,jy)* zs   (iz,jz)
        dumy = xs   (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
        dumz = xs   (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
        djsij(ij,1)=djsij(ij,1)+ dumx*pij
        djsij(ij,2)=djsij(ij,2)+ dumy*pij
        djsij(ij,3)=djsij(ij,3)+ dumz*pij
      endif
c
      if (nder.eq.2) then
      dumxx=dxsdii(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
      dumyy= xs   (ix,jx)*dysdii(iy,jy)* zs   (iz,jz)
      dumzz= xs   (ix,jx)* ys   (iy,jy)*dzsdii(iz,jz)
      dumxy=dxsdi (ix,jx)*dysdi (iy,jy)* zs   (iz,jz)
      dumxz=dxsdi (ix,jx)* ys   (iy,jy)*dzsdi (iz,jz)
      dumyz= xs   (ix,jx)*dysdi (iy,jy)*dzsdi (iz,jz)
      didisij(ij,1,1)=didisij(ij,1,1)+ dumxx*pij 
      didisij(ij,1,2)=didisij(ij,1,2)+ dumxy*pij 
      didisij(ij,1,3)=didisij(ij,1,3)+ dumxz*pij 
      didisij(ij,2,1)=didisij(ij,2,1)+ dumxy*pij 
      didisij(ij,2,2)=didisij(ij,2,2)+ dumyy*pij 
      didisij(ij,2,3)=didisij(ij,2,3)+ dumyz*pij 
      didisij(ij,3,1)=didisij(ij,3,1)+ dumxz*pij 
      didisij(ij,3,2)=didisij(ij,3,2)+ dumyz*pij 
      didisij(ij,3,3)=didisij(ij,3,3)+ dumzz*pij 
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
      didjsij(ij,1,1)=didjsij(ij,1,1)+ dumxx*pij 
      didjsij(ij,1,2)=didjsij(ij,1,2)+ dumxy*pij 
      didjsij(ij,1,3)=didjsij(ij,1,3)+ dumxz*pij 
      didjsij(ij,2,1)=didjsij(ij,2,1)+ dumyx*pij 
      didjsij(ij,2,2)=didjsij(ij,2,2)+ dumyy*pij 
      didjsij(ij,2,3)=didjsij(ij,2,3)+ dumyz*pij 
      didjsij(ij,3,1)=didjsij(ij,3,1)+ dumzx*pij 
      didjsij(ij,3,2)=didjsij(ij,3,2)+ dumzy*pij 
      didjsij(ij,3,3)=didjsij(ij,3,3)+ dumzz*pij 

      dumxx=dxsdjj(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
      dumyy= xs   (ix,jx)*dysdjj(iy,jy)* zs   (iz,jz)
      dumzz= xs   (ix,jx)* ys   (iy,jy)*dzsdjj(iz,jz)
      dumxy=dxsdj (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
      dumxz=dxsdj (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
      dumyz= xs   (ix,jx)*dysdj (iy,jy)*dzsdj (iz,jz)
      djdjsij(ij,1,1)=djdjsij(ij,1,1)+ dumxx*pij 
      djdjsij(ij,1,2)=djdjsij(ij,1,2)+ dumxy*pij 
      djdjsij(ij,1,3)=djdjsij(ij,1,3)+ dumxz*pij 
      djdjsij(ij,2,1)=djdjsij(ij,2,1)+ dumxy*pij 
      djdjsij(ij,2,2)=djdjsij(ij,2,2)+ dumyy*pij 
      djdjsij(ij,2,3)=djdjsij(ij,2,3)+ dumyz*pij 
      djdjsij(ij,3,1)=djdjsij(ij,3,1)+ dumxz*pij 
      djdjsij(ij,3,2)=djdjsij(ij,3,2)+ dumyz*pij 
      djdjsij(ij,3,3)=djdjsij(ij,3,3)+ dumzz*pij 
c
      endif ! nder=2

  380 continue
  390 continue
c
 6000 continue
 7000 continue
c
      return
      end
