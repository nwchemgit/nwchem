c
c $Id$
c
c  Modified from the HONDO dipole integral code
c
      subroutine hnd_dipder(xyzi, expi, coefi,i_nprim,i_ngen, Li,
     *  xyzj, expj, coefj, j_nprim, j_ngen, Lj,
     *  center, buf, lbuf, nint, scr, lscr)
c
c This routine calculates the dipole derivative integrals.  This is a
c wrapper to the work routine hnd_ddipint.
c
      implicit none
c
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer lscr     ! [input] size of scratch array
      integer lbuf     ! [input] size of any integral
      integer nint     ! [input] size of any integral buffer
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision center(3)  ! [input] center for expansion
      double precision scr(lscr)  ! [scratch] scratch buffer
      double precision buf(lbuf)  ! [output] dipole derivative integrals
c
      double precision zero
      data zero   /0.0d+00/
c
# include "hnd_pointers.fh"
c
c Zero out the integral 
c
      call dcopy(nint*18,zero,0,buf,1)
c
c Use scratch for temps in integral calculation
c Scratch needs are 6*(Li+2)*(Lj+2)+12*(Li+1)*(Lj+1)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_ddipint(xyzi, expi, coefi,i_nprim,i_ngen, Li,
     *  xyzj, expj, coefj, j_nprim, j_ngen, Lj,
     *  center, buf, buf((nint*9)+1),nint,
     *  scr(tvp(1)),scr(tvp(2)),scr(tvp(3)),scr(tvp(4)),scr(tvp(5)),
     *  scr(tvp(6)),scr(tvp(7)),scr(tvp(8)),scr(tvp(9)),
     *  scr(tvp(10)),scr(tvp(11)),scr(tvp(12)),scr(tvp(13)),
     *  scr(tvp(14)),scr(tvp(15)),scr(tvp(16)),scr(tvp(17)),
     *  scr(tvp(18)))
c
      return
      end
c
      subroutine hnd_ddipint(xyzi, expi, coefi,i_nprim,i_ngen, Li,
     *  xyzj, expj, coefj, j_nprim, j_ngen, Lj, center, didij, djdij, 
     *  nint, xin, yin,zin, xxin, yyin, zzin, dxsdi, dysdi, dzsdi, 
     *  dxxdi, dyydi, dzzdi, dxsdj, dysdj, dzsdj, dxxdj, dyydj, dzzdj)
c
c This is the main work routine for the dipole derivative integrals.
c It is assumed the didij and djdij have been zeroed before entering.
c
      implicit none
#include "hnd_tol.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer i_nprim  ! [input] num. prims on function i
      integer i_ngen   ! [input] num general conts on func. i
      integer Li       ! [input] angular momentum of func. i
      integer j_nprim  ! [input] num. prims on function j
      integer j_ngen   ! [input] num general conts on func. j
      integer Lj       ! [input] angular momentum of func. j
      integer nint     ! [input] number of base integrals
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim)  ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen)  ! [input] j coeffs
      double precision center(3)  ! [input] center for expansion
      double precision didij(nint,9) ! [output] dip.der. integrals wrt i
      double precision djdij(nint,9) ! [ouptut] dip.der. integrals wrt j
c
      integer ni, nj
      integer ijx, ijy, ijz, ix, iy, iz, jx, jy, jz
      integer nder, lit, ljt, litder, ljtder
      integer maxi, maxj
      integer ig, jg, ij, i, j
      double precision zero, one
      double precision rln10, tol
      double precision rr, ai, arri, axi, ayi, azi
      double precision aj, aa, aa1, dum, fac
      double precision csi, cpi, cdi, cfi, cgi
      double precision csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2
      double precision dij
      double precision xin, yin,zin, xxin, yyin, zzin
      double precision dxsdi, dysdi, dzsdi, dxxdi, dyydi, dzzdi
      double precision dxsdj, dysdj, dzsdj, dxxdj, dyydj, dzzdj
      double precision dumxdx, dumydx, dumzdx
      double precision dumxdy, dumydy, dumzdy
      double precision dumxdz, dumydz, dumzdz
      double precision xint, yint, zint, xintx, yinty, zintz
      double precision t, xc, yc, zc, x0, y0, z0
      double precision xi, yi, zi, xj, yj, zj
      double precision t1
c
      common/hnd_xyzdip/xint,yint,zint,xintx,yinty,zintz,t,xc,yc,zc,
     1              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension   xin(Li+2,Lj+2),  yin(Li+2,Lj+2),  zin(Li+2,Lj+2)
      dimension  xxin(Li+2,Lj+2), yyin(Li+2,Lj+2), zzin(Li+2,Lj+2)
      dimension dxsdi(Li+1,Lj+1),dysdi(Li+1,Lj+1),dzsdi(Li+1,Lj+1)
      dimension dxsdj(Li+1,Lj+1),dysdj(Li+1,Lj+1),dzsdj(Li+1,Lj+1)
      dimension dxxdi(Li+1,Lj+1),dyydi(Li+1,Lj+1),dzzdi(Li+1,Lj+1)
      dimension dxxdj(Li+1,Lj+1),dyydj(Li+1,Lj+1),dzzdj(Li+1,Lj+1)
      dimension Nxyz(3)
      data zero  /0.0d+00/
      data one   /1.0e+00/
      data rln10 /2.30258d+00/
c
      tol=rln10*itol
      nder = 1
c
      xc=center(1)
      yc=center(2)
      zc=center(3)
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
      nroots=(lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_dipder: angular momentum too high',555,
     &       BASIS_ERR)
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
      if(dum.gt.tol) go to 6000  ! the integral is zero
      fac= exp(-dum)
      csj=coefj(jg,j_ngen)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      dij=fac*csi*csj
c
c     ----- dipole moment integrals -----
c
      t = sqrt(aa)
      t1=one/t
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljtder
      nj=j
      do 370 i=1,litder
      ni=i
      call hnd_dipxyz
       xin(i,j)=xint*t1
       yin(i,j)=yint*t1
       zin(i,j)=zint*t1
      xxin(i,j)=xintx*t1
      yyin(i,j)=yinty*t1
      zzin(i,j)=zintz*t1
  370 continue
c
      call hnd_deriaj(dxsdi,dysdi,dzsdi,dxsdj,dysdj,dzsdj,
     *                xin,yin,zin,
     *                lit,ljt,ai,aj)
      call hnd_deriaj(dxxdi,dyydi,dzzdi,dxxdj,dyydj,dzzdj,
     *                xxin,yyin,zzin,
     *                lit,ljt,ai,aj)
      ij=0
      do 390 i=1,maxi
      call getNxyz(Li,i,Nxyz)
      ix = Nxyz(1) + 1
      iy = Nxyz(2) + 1
      iz = Nxyz(3) + 1
      do 390 j=1,maxj
      call getNxyz(Lj,j,Nxyz)
      jx = Nxyz(1) + 1
      jy = Nxyz(2) + 1
      jz = Nxyz(3) + 1
      ij=ij+1
c
c  First do derivative wrt the first atom
c
      dumxdx = dxxdi(ix,jx)*  yin(iy,jy)*  zin(iz,jz)
      dumydx = dxsdi(ix,jx)* yyin(iy,jy)*  zin(iz,jz)
      dumzdx = dxsdi(ix,jx)*  yin(iy,jy)* zzin(iz,jz)
      dumxdy =  xxin(ix,jx)*dysdi(iy,jy)*  zin(iz,jz)
      dumydy =   xin(ix,jx)*dyydi(iy,jy)*  zin(iz,jz)
      dumzdy =   xin(ix,jx)*dysdi(iy,jy)* zzin(iz,jz)
      dumxdz =  xxin(ix,jx)*  yin(iy,jy)*dzsdi(iz,jz)
      dumydz =   xin(ix,jx)* yyin(iy,jy)*dzsdi(iz,jz)
      dumzdz =   xin(ix,jx)*  yin(iy,jy)*dzzdi(iz,jz)
c
      didij(ij,1) = didij(ij,1) + dij*dumxdx
      didij(ij,2) = didij(ij,2) + dij*dumydx
      didij(ij,3) = didij(ij,3) + dij*dumzdx
      didij(ij,4) = didij(ij,4) + dij*dumxdy
      didij(ij,5) = didij(ij,5) + dij*dumydy
      didij(ij,6) = didij(ij,6) + dij*dumzdy
      didij(ij,7) = didij(ij,7) + dij*dumxdz
      didij(ij,8) = didij(ij,8) + dij*dumydz
      didij(ij,9) = didij(ij,9) + dij*dumzdz
c
c  Next do derivative wrt the second atom
c
      dumxdx = dxxdj(ix,jx)*  yin(iy,jy)*  zin(iz,jz)
      dumydx = dxsdj(ix,jx)* yyin(iy,jy)*  zin(iz,jz)
      dumzdx = dxsdj(ix,jx)*  yin(iy,jy)* zzin(iz,jz)
      dumxdy =  xxin(ix,jx)*dysdj(iy,jy)*  zin(iz,jz)
      dumydy =   xin(ix,jx)*dyydj(iy,jy)*  zin(iz,jz)
      dumzdy =   xin(ix,jx)*dysdj(iy,jy)* zzin(iz,jz)
      dumxdz =  xxin(ix,jx)*  yin(iy,jy)*dzsdj(iz,jz)
      dumydz =   xin(ix,jx)* yyin(iy,jy)*dzsdj(iz,jz)
      dumzdz =   xin(ix,jx)*  yin(iy,jy)*dzzdj(iz,jz)
c
      djdij(ij,1) = djdij(ij,1) + dij*dumxdx
      djdij(ij,2) = djdij(ij,2) + dij*dumydx
      djdij(ij,3) = djdij(ij,3) + dij*dumzdx
      djdij(ij,4) = djdij(ij,4) + dij*dumxdy
      djdij(ij,5) = djdij(ij,5) + dij*dumydy
      djdij(ij,6) = djdij(ij,6) + dij*dumzdy
      djdij(ij,7) = djdij(ij,7) + dij*dumxdz
      djdij(ij,8) = djdij(ij,8) + dij*dumydz
      djdij(ij,9) = djdij(ij,9) + dij*dumzdz
  390 continue
c
 6000 continue
 7000 continue
      return
c
 9997 format(' in -hnd_stvint- the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i2,/,
     2       ' lit,ljt,nroots= ',3i3)
      end
