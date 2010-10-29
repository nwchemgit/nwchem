c
c $Id$
c
c Modified from HONDO integral code
c
      subroutine hnd_hlfd_ij2(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1  xyzj,expj,coefj,j_nprim,j_ngen,Lj,kat,didvij,djdvij,zan,xyz,
     2  nder,nint,scr,lscr)
c
c This is a wrapper routine, setting up scratch blocks used in actual
c integral routine
c
      implicit none
c
      integer i_nprim   ! [input] num. prims on function i
      integer i_ngen    ! [input] num general conts on func. i
      integer Li        ! [input] angular momentum of func. i
      integer j_nprim   ! [input] num. prims on function j
      integer j_ngen    ! [input] num general conts on func. j
      integer Lj        ! [input] angular momentum of func. j
      integer nder      ! [input] should always be 2=2nd der
      integer nint      ! [input] number of base integrals
      integer kat       ! [input] lexical number of an atom
      integer lscr      ! [input] length of scratch space
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
      double precision scr(lscr) ! [input] scratch buffer
c
#include "hnd_pointers.fh"
c
c Use scratch for temps in integral calculation
c Scratch needs are 
c (6*(Li+2)*(Lj+2)+12*(Li+1)*(Lj+1))*((Li+Lj+nder)/2+1)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_hlfd_ij21(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1  coefj,j_nprim,j_ngen,Lj,kat,didvij,djdvij,zan,xyz,nder,nint,
     2  scr(hdp(1)) ,scr(hdp(2)) ,scr(hdp(3)) ,scr(hdp(4)),scr(hdp(5)),
     3  scr(hdp(6)) ,scr(hdp(7)) ,scr(hdp(8)) ,scr(hdp(9)),scr(hdp(10)),
     4  scr(hdp(11)),scr(hdp(12)),scr(hdp(13)),scr(hdp(14)),
     5  scr(hdp(15)),scr(hdp(16)),scr(hdp(17)),scr(hdp(18)))
c
      return
      end
c
      subroutine hnd_hlfd_ij21(xyzi,expi,coefi,i_nprim,i_ngen,Li,
     1  xyzj,expj,coefj,j_nprim,j_ngen,Lj,kat,didvij,djdvij,zan,xyz,
     2  nder,nint,xv,yv,zv,dxv,dyv,dzv,dxvdi,dyvdi,dzvdi,dxvdj,dyvdj, 
     3  dzvdj,ddxvdi,ddyvdi,ddzvdi,ddxvdj,ddyvdj,ddzvdj)
c
c  This is a routine that computes additional cross 2nd derivative Helman-
c  Feynman/orbital terms.  It is assumed that the buffers for the integrals
c  have been zeroed before entering this routine.
c
      implicit none
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
#include "errquit.fh"
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
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                             ni,nj,cx,cy,cz
      double precision xint, yint, zint, tx, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz
      double precision rln10, zero, one, tol, rr, ai, aj, arri
      double precision axi, ayi, azi, csi, cpi, cdi, cfi, cgi
      double precision aa, aa1, dum, fac, csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2, pij
      double precision dumxx, dumyy, dumzz, dumxy, dumxz, dumyz
      double precision dumyx, dumzx, dumzy
      double precision two, pi212, aax, aay, aaz, znuc
      double precision uu, ww, w2, tt, xv, yv, zv
      double precision dxv, dyv, dzv
      double precision dxvdi, dyvdi, dzvdi
      double precision dxvdj, dyvdj, dzvdj
      double precision ddxvdi, ddyvdi, ddzvdi
      double precision ddxvdj, ddyvdj, ddzvdj
      integer ni, nj, iroot, Nxyz
      integer lit, maxi, ljt,maxj
      integer ig, jg, ij, i, j, ix, iy, iz, jx, jy, jz
      integer litder, ljtder
      dimension w2(maxrys), Nxyz(3)
      dimension    xv(Li+2,Lj+2,*),   yv(Li+2,Lj+2,*),   zv(Li+2,Lj+2,*)
      dimension   dxv(Li+2,Lj+2,*),  dyv(Li+2,Lj+2,*),  dzv(Li+2,Lj+2,*)
      dimension dxvdi(Li+1,Lj+1,*),dyvdi(Li+1,Lj+1,*),dzvdi(Li+1,Lj+1,*)
      dimension dxvdj(Li+1,Lj+1,*),dyvdj(Li+1,Lj+1,*),dzvdj(Li+1,Lj+1,*)
      dimension ddxvdi(Li+1,Lj+1,*),ddyvdi(Li+1,Lj+1,*)
      dimension ddzvdi(Li+1,Lj+1,*),ddxvdj(Li+1,Lj+1,*)
      dimension ddyvdj(Li+1,Lj+1,*),ddzvdj(Li+1,Lj+1,*)
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
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
      lit=Li + 1
      maxi=lit*(lit+1)/2
c
      litder=lit+nder-1
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      ljt=Lj + 1
      maxj=ljt*(ljt+1)/2
c
      ljtder=ljt+nder-1
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_hlfd_ij2: maxrys is too large!',555, INT_ERR)
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
      pij=fac*csi*csj*aa1*pi212
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
      yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      call hnd_droot
      do 420 iroot=1,nroots
      uu=u9(iroot)*aa
      ww=w9(iroot)*znuc
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
      call hnd_dervxyz(1)
      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint
  410 continue
c
      call hnd_deriaj( dxvdi(1,1,iroot), dyvdi(1,1,iroot),
     1              dzvdi(1,1,iroot), dxvdj(1,1,iroot),
     1              dyvdj(1,1,iroot), dzvdj(1,1,iroot),
     1              xv(1,1,iroot),    yv(1,1,iroot),   zv(1,1,iroot),
     2          lit,ljt,ai,aj)
      call hnd_deriaj(ddxvdi(1,1,iroot),ddyvdi(1,1,iroot),
     1             ddzvdi(1,1,iroot), ddxvdj(1,1,iroot),
     1             ddyvdj(1,1,iroot), ddzvdj(1,1,iroot),
     1             dxv(1,1,iroot),   dyv(1,1,iroot),  dzv(1,1,iroot),
     2          lit,ljt,ai,aj)
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
     1  (ddxvdi(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxy=dumxy+w2(iroot)*
     1  (dxvdi (ix,jx,iroot)*dyv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxz=dumxz+w2(iroot)*
     1  (dxvdi (ix,jx,iroot)* yv   (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumyx=dumyx+w2(iroot)*
     1  (dxv   (ix,jx,iroot)*dyvdi (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyy=dumyy+w2(iroot)*
     1  ( xv   (ix,jx,iroot)*ddyvdi(iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyz=dumyz+w2(iroot)*
     1  ( xv   (ix,jx,iroot)*dyvdi (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumzx=dumzx+w2(iroot)*
     1  (dxv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdi (iz,jz,iroot))
      dumzy=dumzy+w2(iroot)*
     1  ( xv   (ix,jx,iroot)*dyv   (iy,jy,iroot)*dzvdi (iz,jz,iroot))
      dumzz=dumzz+w2(iroot)*
     1  ( xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*ddzvdi(iz,jz,iroot))
  430 continue
      ij=ij+1
      didvij(ij,1,1)=didvij(ij,1,1)+dumxx*pij
      didvij(ij,1,2)=didvij(ij,1,2)+dumxy*pij
      didvij(ij,1,3)=didvij(ij,1,3)+dumxz*pij
      didvij(ij,2,1)=didvij(ij,2,1)+dumyx*pij
      didvij(ij,2,2)=didvij(ij,2,2)+dumyy*pij
      didvij(ij,2,3)=didvij(ij,2,3)+dumyz*pij
      didvij(ij,3,1)=didvij(ij,3,1)+dumzx*pij
      didvij(ij,3,2)=didvij(ij,3,2)+dumzy*pij
      didvij(ij,3,3)=didvij(ij,3,3)+dumzz*pij
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
     1  (ddxvdj(ix,jx,iroot)* yv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxy=dumxy+w2(iroot)*
     1  (dxvdj (ix,jx,iroot)*dyv   (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumxz=dumxz+w2(iroot)*
     1  (dxvdj (ix,jx,iroot)* yv   (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumyx=dumyx+w2(iroot)*
     1  (dxv   (ix,jx,iroot)*dyvdj (iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyy=dumyy+w2(iroot)*
     1  ( xv   (ix,jx,iroot)*ddyvdj(iy,jy,iroot)* zv   (iz,jz,iroot))
      dumyz=dumyz+w2(iroot)*
     1  ( xv   (ix,jx,iroot)*dyvdj (iy,jy,iroot)*dzv   (iz,jz,iroot))
      dumzx=dumzx+w2(iroot)*
     1  (dxv   (ix,jx,iroot)* yv   (iy,jy,iroot)*dzvdj (iz,jz,iroot))
      dumzy=dumzy+w2(iroot)*
     1  ( xv   (ix,jx,iroot)*dyv   (iy,jy,iroot)*dzvdj (iz,jz,iroot))
      dumzz=dumzz+w2(iroot)*
     1  ( xv   (ix,jx,iroot)* yv   (iy,jy,iroot)*ddzvdj(iz,jz,iroot))
  435 continue
      djdvij(ij,1,1)=djdvij(ij,1,1)+dumxx*pij
      djdvij(ij,1,2)=djdvij(ij,1,2)+dumxy*pij
      djdvij(ij,1,3)=djdvij(ij,1,3)+dumxz*pij
      djdvij(ij,2,1)=djdvij(ij,2,1)+dumyx*pij
      djdvij(ij,2,2)=djdvij(ij,2,2)+dumyy*pij
      djdvij(ij,2,3)=djdvij(ij,2,3)+dumyz*pij
      djdvij(ij,3,1)=djdvij(ij,3,1)+dumzx*pij
      djdvij(ij,3,2)=djdvij(ij,3,2)+dumzy*pij
      djdvij(ij,3,3)=djdvij(ij,3,3)+dumzz*pij
c
  440 continue
  450 continue
c
 6000 continue
 7000 continue
c
      return
 9997 format(' in -hlfd- , the rys quadrature is not implememented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots= ',3i3)
      end
