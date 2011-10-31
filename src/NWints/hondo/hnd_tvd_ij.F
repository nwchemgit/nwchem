c
c $Id$
c
c Modified from HONDO
c
      subroutine hnd_tvd_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj, xyz,
     2  ditij,djtij,diditij,didjtij,djdjtij,divij,djvij,
     3  didivij,didjvij,djdjvij,nder,nint,zan,nat,doT,doV,
     4  iatom,jatom,scr,lscr)
c
c This is a wrapper routine, setting up scratch blocks used in actual
c integral routine
c
      implicit none
c
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
      integer lscr     ! [input] length of scratch space
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,nat) ! [input] all atom positions
      double precision ditij(nint,3) ! [output] 1rst derivative wrt i for T
      double precision djtij(nint,3) ! [output] 1rst derivative wrt j for T
      double precision diditij(nint,3,3) ! [output] 2nd derivative wrt i,i for T
      double precision didjtij(nint,3,3) ! [output] 2nd derivative wrt i,j for T
      double precision djdjtij(nint,3,3) ! [output] 2nd derivative wrt j,j for T
      double precision divij(nint,3) ! [output] 1rst derivative wrt i for V
      double precision djvij(nint,3) ! [output] 1rst derivative wrt j for V
      double precision didivij(nint,3,3) ! [output] 2nd derivative wrt i,i for V
      double precision didjvij(nint,3,3) ! [output] 2nd derivative wrt i,j for V
      double precision djdjvij(nint,3,3) ! [output] 2nd derivative wrt j,j for V
      double precision zan(nat) ! [input] nuclear charges
      double precision scr(lscr) ! [input] scratch buffer
c
#include "hnd_pointers.fh"     
c
c Use scratch for temps in integral calculation
c Scratch needs are dependent on nder and doT and doV:
c doT:
c  nder=1: 3*(Li+3)*(2*Lj+8)+12*(Li+1)*(Lj+1)
c  nder=2: 3*(Li+3)*(2*Lj+8)+30*(Li+1)*(Lj+1)
c doV:
c  nder=1: (3*(Li+3)*(Lj+3)+6*(Li+1)*(Lj+1)) * ((Li+Lj+2)/2+1)
c  nder=2: (3*(Li+3)*(Lj+3)+6*(Li+1)*(Lj+1)) * ((Li+Lj+2)/2+1) 
c                          +9*(Li+1)*(Lj+1)  * ((Li+Lj+2)/2+1)
c
c The pointers are defined in hnd_pointers.fh
c
      call hnd_tvd_ij21(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     *  xyzj,expj,coefj, j_nprim, j_ngen, Lj, xyz,
     *  ditij,djtij,diditij,didjtij,djdjtij,divij,djvij,
     *  didivij,didjvij,djdjvij,nder,nint,zan,nat,doT,doV,
     *  iatom,jatom,
     *  scr(tvp(1)) ,scr(tvp(2)) ,scr(tvp(3)) ,scr(tvp(4)),
     *  scr(tvp(5)) ,scr(tvp(6)) ,scr(tvp(7)) ,scr(tvp(8)),
     *  scr(tvp(9)) ,scr(tvp(10)),scr(tvp(11)),scr(tvp(12)),
     *  scr(tvp(13)),scr(tvp(14)),scr(tvp(15)),scr(tvp(16)),
     *  scr(tvp(17)),scr(tvp(18)),scr(tvp(19)),scr(tvp(20)),
     *  scr(tvp(21)),scr(tvp(22)),scr(tvp(23)),scr(tvp(24)),
     *  scr(tvp(25)),scr(tvp(26)),scr(tvp(27)),scr(tvp(28)),
     *  scr(tvp(29)),scr(tvp(30)),scr(tvp(31)),scr(tvp(32)),
     *  scr(tvp(33)),scr(tvp(34)),scr(tvp(35)),scr(tvp(36)),
     *  scr(tvp(37)),scr(tvp(38)),scr(tvp(39)),scr(tvp(40)),
     *  scr(tvp(41)),scr(tvp(42)),scr(tvp(43)),scr(tvp(44)),
     *  scr(tvp(45)),scr(tvp(46)),scr(tvp(47)),scr(tvp(48)),
     *  scr(tvp(49)),scr(tvp(50)),scr(tvp(51)),scr(tvp(52)),
     *  scr(tvp(53)),scr(tvp(54)))
c
      return
      end
c
      subroutine hnd_tvd_ij21(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     *  xyzj,expj,coefj, j_nprim, j_ngen, Lj, xyz,
     *  ditij,djtij,diditij,didjtij,djdjtij,divij,djvij,
     *  didivij,didjvij,djdjvij,nder,nint,zan,nat,doT,doV,
     *  iatom,jatom,
     *  xs,ys,zs,xt,yt,zt,dxsdi,dysdi,dzsdi,dxsdj,dysdj,dzsdj,dxtdi,
     *  dytdi,dztdi,dxtdj,dytdj,dztdj,xv,yv,zv,dxvdi,dyvdi,dzvdi,dxvdj,
     *  dyvdj,dzvdj,dxsdii,dysdii,dzsdii,dxsdjj,dysdjj,dzsdjj,dxsdij,
     *  dysdij,dzsdij,dxtdii,dytdii,dztdii,dxtdij,dytdij,dztdij,dxtdjj,
     *  dytdjj,dztdjj,dxvdii,dyvdii,dzvdii,dxvdjj,dyvdjj,dzvdjj,dxvdij,
     *  dyvdij,dzvdij) 
c
      implicit none
#include "stdio.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "errquit.fh"
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
      double precision ditij(nint,3) ! [output] 1rst derivative wrt i for T
      double precision djtij(nint,3) ! [output] 1rst derivative wrt j for T
      double precision diditij(nint,3,3) ! [output] 2nd derivative wrt i,i for T
      double precision didjtij(nint,3,3) ! [output] 2nd derivative wrt i,j for T
      double precision djdjtij(nint,3,3) ! [output] 2nd derivative wrt j,j for T
      double precision divij(nint,3) ! [output] 1rst derivative wrt i for V
      double precision djvij(nint,3) ! [output] 1rst derivative wrt j for V
      double precision didivij(nint,3,3) ! [output] 2nd derivative wrt i,i for V
      double precision didjvij(nint,3,3) ! [output] 2nd derivative wrt i,j for V
      double precision djdjvij(nint,3,3) ! [output] 2nd derivative wrt j,j for V
      double precision zan(nat) ! [input] nuclear charges
c
      double precision xint, yint, zint, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision xs, ys, zs
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision pi212, xt, yt, zt, dxsdii, dysdii, dzsdii
      double precision dxsdjj, dysdjj, dzsdjj
      double precision dxsdij, dysdij, dzsdij, dxsdi, dysdi, dzsdi
      double precision dxsdj, dysdj, dzsdj, dxtdii, dytdii, dztdii
      double precision dxtdjj, dytdjj, dztdjj
      double precision dxtdij, dytdij, dztdij, dxtdi, dytdi, dztdi
      double precision dxtdj, dytdj, dztdj, dumx, dumy, dumz
      double precision dumxx, dumyy, dumzz, dumyx, dumzx, dumzy
      double precision dumxy, dumxz, dumyz, aax, aay, aaz
      double precision uu, ww, tt, xv, yv, zv
      double precision dxvdii, dyvdii, dzvdii, dxvdij, dyvdij, dzvdij
      double precision dxvdjj, dyvdjj, dzvdjj
      double precision dxvdi, dyvdi, dzvdi, dxvdj, dyvdj, dzvdj
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision tx, aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      integer ni, nj, lit, Nxyz
      integer maxi, litder, ljt, maxj 
      integer ljtder, ljtmod, ig, jg, ij, i, j, ix, iy, iz
      integer jx, jy, jz, kat, znuc, iroot
c
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      dimension xs(Li+3,Lj+5),ys(Li+3,Lj+5),zs(Li+3,Lj+5)
      dimension dxsdi(Li+1,Lj+1),dysdi(Li+1,Lj+1),dzsdi(Li+1,Lj+1)
      dimension dxsdj(Li+1,Lj+1),dysdj(Li+1,Lj+1),dzsdj(Li+1,Lj+1)
      dimension dxsdii(Li+1,Lj+1),dysdii(Li+1,Lj+1),dzsdii(Li+1,Lj+1)
      dimension dxsdjj(Li+1,Lj+1),dysdjj(Li+1,Lj+1),dzsdjj(Li+1,Lj+1)
      dimension dxsdij(Li+1,Lj+1),dysdij(Li+1,Lj+1),dzsdij(Li+1,Lj+1)
      dimension xt(Li+3,Lj+3),yt(Li+3,Lj+3),zt(Li+3,Lj+3)
      dimension dxtdi(Li+1,Lj+1),dytdi(Li+1,Lj+1),dztdi(Li+1,Lj+1)
      dimension dxtdj(Li+1,Lj+1),dytdj(Li+1,Lj+1),dztdj(Li+1,Lj+1)
      dimension dxtdii(Li+1,Lj+1),dytdii(Li+1,Lj+1),dztdii(Li+1,Lj+1)
      dimension dxtdij(Li+1,Lj+1),dytdij(Li+1,Lj+1),dztdij(Li+1,Lj+1)
      dimension dxtdjj(Li+1,Lj+1),dytdjj(Li+1,Lj+1),dztdjj(Li+1,Lj+1)
      dimension xv(Li+3,Lj+3,*),yv(Li+3,Lj+3,*),zv(Li+3,Lj+3,*)
      dimension dxvdi(Li+1,Lj+1,*),dyvdi(Li+1,Lj+1,*),dzvdi(Li+1,Lj+1,*)
      dimension dxvdj(Li+1,Lj+1,*),dyvdj(Li+1,Lj+1,*),dzvdj(Li+1,Lj+1,*)
      dimension dxvdii(Li+1,Lj+1,*),dyvdii(Li+1,Lj+1,*)
      dimension dzvdii(Li+1,Lj+1,*),dxvdjj(Li+1,Lj+1,*)
      dimension dyvdjj(Li+1,Lj+1,*),dzvdjj(Li+1,Lj+1,*)
      dimension dxvdij(Li+1,Lj+1,*),dyvdij(Li+1,Lj+1,*)
      dimension dzvdij(Li+1,Lj+1,*)
      dimension Nxyz(3)
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pi212  /1.1283791670955d+00/
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
      ljtmod=ljtder+2
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_tvd_ij: maxrys is to large!',555, INT_ERR)
      endif
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
      pij=fac*csi*csj
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
     1            dzsdij,dxsdjj,dysdjj,dzsdjj,
     1            dxsdi ,dysdi ,dzsdi ,dxsdj ,dysdj ,dzsdj ,
     1            xs,ys,zs,lit,ljt,ai,aj,nder)
c
      call hnd_dderij2(dxtdii,dytdii,dztdii,dxtdij,dytdij,
     1            dztdij,dxtdjj,dytdjj,dztdjj,
     1            dxtdi ,dytdi ,dztdi ,dxtdj ,dytdj ,dztdj ,
     1            xt,yt,zt,lit,ljt,ai,aj,nder)
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
      ditij(ij,1)=ditij(ij,1)+ dumx*pij
      ditij(ij,2)=ditij(ij,2)+ dumy*pij
      ditij(ij,3)=ditij(ij,3)+ dumz*pij
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
      djtij(ij,1)=djtij(ij,1)+ dumx*pij
      djtij(ij,2)=djtij(ij,2)+ dumy*pij
      djtij(ij,3)=djtij(ij,3)+ dumz*pij
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
      diditij(ij,1,1)=diditij(ij,1,1)+ dumxx*pij
      diditij(ij,1,2)=diditij(ij,1,2)+ dumyx*pij
      diditij(ij,1,3)=diditij(ij,1,3)+ dumzx*pij
      diditij(ij,2,1)=diditij(ij,2,1)+ dumyx*pij
      diditij(ij,2,2)=diditij(ij,2,2)+ dumyy*pij
      diditij(ij,2,3)=diditij(ij,2,3)+ dumzy*pij
      diditij(ij,3,1)=diditij(ij,3,1)+ dumzx*pij
      diditij(ij,3,2)=diditij(ij,3,2)+ dumzy*pij
      diditij(ij,3,3)=diditij(ij,3,3)+ dumzz*pij
c
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
      didjtij(ij,1,1)=didjtij(ij,1,1)+ dumxx*pij
      didjtij(ij,1,2)=didjtij(ij,1,2)+ dumxy*pij
      didjtij(ij,1,3)=didjtij(ij,1,3)+ dumxz*pij
      didjtij(ij,2,1)=didjtij(ij,2,1)+ dumyx*pij
      didjtij(ij,2,2)=didjtij(ij,2,2)+ dumyy*pij
      didjtij(ij,2,3)=didjtij(ij,2,3)+ dumyz*pij
      didjtij(ij,3,1)=didjtij(ij,3,1)+ dumzx*pij
      didjtij(ij,3,2)=didjtij(ij,3,2)+ dumzy*pij
      didjtij(ij,3,3)=didjtij(ij,3,3)+ dumzz*pij
c
      dumxx=dxtdjj(ix,jx)* ys   (iy,jy)* zs   (iz,jz)
     1     +dxsdjj(ix,jx)* yt   (iy,jy)* zs   (iz,jz)
     2     +dxsdjj(ix,jx)* ys   (iy,jy)* zt   (iz,jz)
      dumyy= xt   (ix,jx)*dysdjj(iy,jy)* zs   (iz,jz)
     1     + xs   (ix,jx)*dytdjj(iy,jy)* zs   (iz,jz)
     2     + xs   (ix,jx)*dysdjj(iy,jy)* zt   (iz,jz)
      dumzz= xt   (ix,jx)* ys   (iy,jy)*dzsdjj(iz,jz)
     1     + xs   (ix,jx)* yt   (iy,jy)*dzsdjj(iz,jz)
     2     + xs   (ix,jx)* ys   (iy,jy)*dztdjj(iz,jz)
      dumyx=dxtdj (ix,jx)*dysdj (iy,jy)* zs   (iz,jz)
     1     +dxsdj (ix,jx)*dytdj (iy,jy)* zs   (iz,jz)
     2     +dxsdj (ix,jx)*dysdj (iy,jy)* zt   (iz,jz)
      dumzx=dxtdj (ix,jx)* ys   (iy,jy)*dzsdj (iz,jz)
     1     +dxsdj (ix,jx)* yt   (iy,jy)*dzsdj (iz,jz)
     2     +dxsdj (ix,jx)* ys   (iy,jy)*dztdj (iz,jz)
      dumzy= xt   (ix,jx)*dysdj (iy,jy)*dzsdj (iz,jz)
     1     + xs   (ix,jx)*dytdj (iy,jy)*dzsdj (iz,jz)
     2     + xs   (ix,jx)*dysdj (iy,jy)*dztdj (iz,jz)
      djdjtij(ij,1,1)=djdjtij(ij,1,1)+ dumxx*pij
      djdjtij(ij,1,2)=djdjtij(ij,1,2)+ dumyx*pij
      djdjtij(ij,1,3)=djdjtij(ij,1,3)+ dumzx*pij
      djdjtij(ij,2,1)=djdjtij(ij,2,1)+ dumyx*pij
      djdjtij(ij,2,2)=djdjtij(ij,2,2)+ dumyy*pij
      djdjtij(ij,2,3)=djdjtij(ij,2,3)+ dumzy*pij
      djdjtij(ij,3,1)=djdjtij(ij,3,1)+ dumzx*pij
      djdjtij(ij,3,2)=djdjtij(ij,3,2)+ dumzy*pij
      djdjtij(ij,3,3)=djdjtij(ij,3,3)+ dumzz*pij
c
      else ! if nder.gt.2
        write(luout,9996) nder
        call errquit('hnd_tvd_ij: derivative order is too high',555,
     &       INT_ERR)
      endif

  380 continue
  390 continue
      endif ! doT
c
c     ----- nuclear attraction derivatives -----
c
      if (doV) then
      pij=pij*aa1*pi212
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 kat=1,nat
      znuc=-zan(kat)
      cx=xyz(1,kat)
      cy=xyz(2,kat)
      cz=xyz(3,kat)
      yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      call hnd_droot
      do 420 iroot=1,nroots
      uu=u9(iroot)*aa
      ww=w9(iroot)*znuc
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
     1            dxvdjj(1,1,iroot),dyvdjj(1,1,iroot),dzvdjj(1,1,iroot),
     2             dxvdi(1,1,iroot), dyvdi(1,1,iroot), dzvdi(1,1,iroot),
     3             dxvdj(1,1,iroot), dyvdj(1,1,iroot), dzvdj(1,1,iroot),
     4                xv(1,1,iroot),    yv(1,1,iroot),    zv(1,1,iroot),
     5            lit,ljt,ai,aj,nder)
c
  420 continue
c
      ij=0
      do 450 i=1,maxi
      call getNxyz(Li,i,Nxyz)
      ix = Nxyz(1) + 1
      iy = Nxyz(2) + 1
      iz = Nxyz(3) + 1
      do 440 j=1,maxj
      call getNxyz(Lj,j,Nxyz)
      jx = Nxyz(1) + 1
      jy = Nxyz(2) + 1
      jz = Nxyz(3) + 1
      ij=ij+1
c
      if (nder.eq.1) then
c
c first do the derivative wrt to the first atom
c
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
      divij(ij,1)=divij(ij,1)+dumx*pij
      divij(ij,2)=divij(ij,2)+dumy*pij
      divij(ij,3)=divij(ij,3)+dumz*pij
c
c first do the derivative wrt to the second atom
c
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
      djvij(ij,1)=djvij(ij,1)+dumx*pij
      djvij(ij,2)=djvij(ij,2)+dumy*pij
      djvij(ij,3)=djvij(ij,3)+dumz*pij
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
      didivij(ij,1,1)=didivij(ij,1,1)+dumxx*pij
      didivij(ij,1,2)=didivij(ij,1,2)+dumyx*pij
      didivij(ij,1,3)=didivij(ij,1,3)+dumzx*pij
      didivij(ij,2,1)=didivij(ij,2,1)+dumyx*pij
      didivij(ij,2,2)=didivij(ij,2,2)+dumyy*pij
      didivij(ij,2,3)=didivij(ij,2,3)+dumzy*pij
      didivij(ij,3,1)=didivij(ij,3,1)+dumzx*pij
      didivij(ij,3,2)=didivij(ij,3,2)+dumzy*pij
      didivij(ij,3,3)=didivij(ij,3,3)+dumzz*pij
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
      didjvij(ij,1,1)=didjvij(ij,1,1)+dumxx*pij 
      didjvij(ij,1,2)=didjvij(ij,1,2)+dumxy*pij 
      didjvij(ij,1,3)=didjvij(ij,1,3)+dumxz*pij 
      didjvij(ij,2,1)=didjvij(ij,2,1)+dumyx*pij 
      didjvij(ij,2,2)=didjvij(ij,2,2)+dumyy*pij 
      didjvij(ij,2,3)=didjvij(ij,2,3)+dumyz*pij 
      didjvij(ij,3,1)=didjvij(ij,3,1)+dumzx*pij 
      didjvij(ij,3,2)=didjvij(ij,3,2)+dumzy*pij 
      didjvij(ij,3,3)=didjvij(ij,3,3)+dumzz*pij 
c
      dumxx=zero
      dumyy=zero
      dumzz=zero
      dumyx=zero
      dumzx=zero
      dumzy=zero
      do 436 iroot=1,nroots
      dumxx=dumxx
     1     +dxvdjj(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumyy=dumyy
     1     + xv   (ix,jx,iroot)*dyvdjj(iy,jy,iroot)* zv   (iz,jz,iroot)
      dumzz=dumzz
     1     + xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdjj(iz,jz,iroot)
      dumyx=dumyx
     1     +dxvdj (ix,jx,iroot)*dyvdj (iy,jy,iroot)* zv   (iz,jz,iroot)
      dumzx=dumzx
     1     +dxvdj (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdj (iz,jz,iroot)
      dumzy=dumzy
     1     + xv   (ix,jx,iroot)*dyvdj (iy,jy,iroot)*dzvdj (iz,jz,iroot)
  436 continue
      djdjvij(ij,1,1)=djdjvij(ij,1,1)+dumxx*pij
      djdjvij(ij,1,2)=djdjvij(ij,1,2)+dumyx*pij
      djdjvij(ij,1,3)=djdjvij(ij,1,3)+dumzx*pij
      djdjvij(ij,2,1)=djdjvij(ij,2,1)+dumyx*pij
      djdjvij(ij,2,2)=djdjvij(ij,2,2)+dumyy*pij
      djdjvij(ij,2,3)=djdjvij(ij,2,3)+dumzy*pij
      djdjvij(ij,3,1)=djdjvij(ij,3,1)+dumzx*pij
      djdjvij(ij,3,2)=djdjvij(ij,3,2)+dumzy*pij
      djdjvij(ij,3,3)=djdjvij(ij,3,3)+dumzz*pij
c
      else ! if nder.gt.2
        write(luout,9996) nder
        call errquit('hnd_tvd_ij: derivative order is too high',555,
     &       INT_ERR)
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
