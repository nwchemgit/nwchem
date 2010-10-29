c
c $Id$
c
c Modified from HONDO
c
      subroutine hnd_hlf_ij2(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj, kat, dvij, ddvij,
     2  zan, xyz, nder, nint, scr,lscr)
c
c This is a wrapper routine, setting up scratch blocks used in actual
c integral routine
c
      implicit none
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nder     ! [input] 1=1rst der; 2=2nd der
      integer nint     ! [input] number of base integrals
      integer kat      ! [input] lexical number of an atom
      integer lscr     ! [input] length of scratch space
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,*)  ! [input] all atom positions
      double precision dvij(nint,3) ! [output] 1rst. der integrals
      double precision ddvij(nint,3,3) ! [output] 2nd der integrals
      double precision zan(*)   ! [input] nuclear charges
      double precision scr(lscr) ! [input] scratch buffer
c
# include "hnd_pointers.fh"
c
c Use scratch for temps in integral calculation
c Scratch needs are dependent on nder:
c nder=1 (6 * (Li+1)*(Lj+1)*((Li+Lj+2)/2+1) + 3 0.0d0 length blocks
c nder=2 (9 * (Li+1)*(Lj+1)*((Li+Lj+2)/2+1)
c
c The pointers are defined in hnd_pointers.fh
c
      call hnd_hlf_ij21(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     2     xyzj,expj,coefj, j_nprim, j_ngen, Lj, kat, dvij, ddvij,
     3     zan, xyz, nder, nint,
     4     scr(tvp(22)),scr(tvp(23)),scr(tvp(24)),
     5     scr(tvp(25)),scr(tvp(26)),scr(tvp(27)),
     6     scr(tvp(46)),scr(tvp(47)),scr(tvp(48)))
c
      return
      end
c
      subroutine hnd_hlf_ij21(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj, j_nprim, j_ngen, Lj, kat, dvij, ddvij,
     2  zan, xyz, nder, nint,xv,yv,zv,dxv,dyv,dzv,ddxv,ddyv,ddzv)
c
c  This is the routine that actually computes the 1rst and 2nd
c  derivatives for the Helman-Feynman term.  It is assumed that the buffers
c  for the integrals have been 0.0d0ed before entering this routine.
c
      implicit none
#include "stdio.fh"
#include "hnd_rys.fh"
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
      integer kat      ! [input] lexical number of an atom
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision xyz(3,*)  ! [input] all atom positions
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
      double precision rln10, tol, rr, ai, aj, arri
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      double precision dumx, dumy, dumz, dumxx, dumyy, dumzz
      double precision dumxy, dumxz, dumyz
      double precision pi212, aax, aay, aaz, znuc
      double precision uu, u2, u4, ww, w2, w4, tt, xv, yv, zv
      double precision dxv, dyv, dzv, ddxv, ddyv, ddzv
      integer ni, nj, iroot, Nxyz
      integer lit,maxi, ljt, maxj
      integer ig, jg, ij, i, j, ix, iy, iz, jx, jy, jz
      dimension   xv(Li+1,Lj+1,*),  yv(Li+1,Lj+1,*),  zv(Li+1,Lj+1,*)
      dimension  dxv(Li+1,Lj+1,*), dyv(Li+1,Lj+1,*), dzv(Li+1,Lj+1,*)
      dimension ddxv(Li+1,Lj+1,*),ddyv(Li+1,Lj+1,*),ddzv(Li+1,Lj+1,*)
      dimension Nxyz(3)
      dimension w2(maxrys),w4(maxrys)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-hlfspd-'/
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
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
      lit = Li + 1
      maxi=lit*(lit+1)/2
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      ljt = Lj + 1
      maxj=ljt*(ljt+1)/2
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
      aa1=1.0d0/aa
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
      pij=fac*csi*csj*pi212*aa1
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
      yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      call hnd_droot
      do 420 iroot=1,nroots
      uu=u9(iroot)*aa
      u2=uu
      u4=uu*uu
      ww=w9(iroot)*znuc
      w2(iroot)=ww*u2*2.0d0
      w4(iroot)=ww*u4*4.0d0
      tt=1.0d0/(aa+uu)
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

      call hnd_dervxyz(1)

      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint

      if (nder .gt. 1) then 
         call hnd_dervxyz(2)

         ddxv(i,j,iroot)=xint
         ddyv(i,j,iroot)=yint
         ddzv(i,j,iroot)=zint
      endif
  410 continue
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
      if (nder.eq.1) then
        dumx=0.0d0
        dumy=0.0d0
        dumz=0.0d0
        do 430 iroot=1,nroots
           dumx = dumx + dxv(ix,jx,iroot)*  yv(iy,jy,iroot)*  
     1             zv(iz,jz,iroot)*w2(iroot)
           dumy = dumy + xv(ix,jx,iroot)*  dyv(iy,jy,iroot)*  
     1             zv(iz,jz,iroot)*w2(iroot)
           dumz = dumz + xv(ix,jx,iroot)*   yv(iy,jy,iroot)* 
     1            dzv(iz,jz,iroot)*w2(iroot)
  430   continue
        dumx =dumx *pij
        dumy =dumy *pij
        dumz =dumz *pij
        dvij(ij,1)=dvij(ij,1)+dumx
        dvij(ij,2)=dvij(ij,2)+dumy
        dvij(ij,3)=dvij(ij,3)+dumz
      elseif (nder.eq.2) then
        dumxx=0.0d0
        dumyy=0.0d0
        dumzz=0.0d0
        dumxy=0.0d0
        dumxz=0.0d0
        dumyz=0.0d0
        do 431 iroot=1,nroots
           dum=xv(ix,jx,iroot)*yv(iy,jy,iroot)*zv(iz,jz,iroot)*w2(iroot)
           dumxx=dumxx-dum+ ddxv(ix,jx,iroot)*  yv(iy,jy,iroot)*  
     1             zv(iz,jz,iroot)*w4(iroot)
           dumyy=dumyy-dum+   xv(ix,jx,iroot)*ddyv(iy,jy,iroot)*  
     1             zv(iz,jz,iroot)*w4(iroot)
           dumzz=dumzz-dum+   xv(ix,jx,iroot)*  yv(iy,jy,iroot)*
     1           ddzv(iz,jz,iroot)*w4(iroot)
           dumxy=dumxy+ dxv(ix,jx,iroot)* dyv(iy,jy,iroot)*  
     1             zv(iz,jz,iroot)*w4(iroot)
           dumxz=dumxz+ dxv(ix,jx,iroot)*  yv(iy,jy,iroot)* 
     1            dzv(iz,jz,iroot)*w4(iroot)
           dumyz=dumyz+  xv(ix,jx,iroot)* dyv(iy,jy,iroot)*
     1            dzv(iz,jz,iroot)*w4(iroot)
  431   continue
        dumxx=dumxx*pij
        dumyy=dumyy*pij
        dumzz=dumzz*pij
        dumxy=dumxy*pij
        dumxz=dumxz*pij
        dumyz=dumyz*pij
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
