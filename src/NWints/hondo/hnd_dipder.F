c
c $Id: hnd_dipder.F,v 1.1 2002-01-20 05:52:52 windus Exp $
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
      integer index, i
      integer mini, maxi, minj, maxj, minprim, maxprim
      dimension minprim(5), maxprim(5)
      double precision zero
      data zero   /0.0d+00/
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 4,10,20,35/
c
c Zero out the integral 
c
      call dfill(nint*18,zero,buf,1)
c
      call hnd_ddipint(xyzi, expi, coefi,i_nprim,i_ngen, Li,
     *  xyzj, expj, coefj, j_nprim, j_ngen, Lj,
     *  center, buf, buf((nint*9)+1),nint)
c
c
c     ----- perform transformation to NWChem order if necessary ----
c
      if ((Li.gt.1).or.(Lj.gt.1)) then
        mini = minprim(Li+1)
        maxi = maxprim(Li+1)
        minj = minprim(Lj+1)
        maxj = maxprim(Lj+1)
        index = 1
        do i=1,6
          call hnd_to_nw_order(scr,buf(index),minj,mini,
     1                         maxj-minj+1,maxi-mini+1)
          call dcopy (nint,scr,1,buf(index),1)
          index = index + nint
        enddo
      endif
c
      return
      end
c
      subroutine hnd_ddipint(xyzi, expi, coefi,i_nprim,i_ngen, Li,
     *  xyzj, expj, coefj, j_nprim, j_ngen, Lj,
     *  center, didij, djdij, nint)
c
c This is the main work routine for the dipole derivative integrals.
c It is assumed the didij and djdij have been zeroed before entering.
c
      implicit none
#include "hnd_tol.fh"
#include "stdio.fh"
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
      integer ni, nj, minprim, maxprim
      integer ijx, ijy, ijz, ix, iy, iz, jx, jy, jz
      integer nder, i1, i2, j1, j2, lit, ljt, litder, ljtder
      integer mini, maxi, minj, maxj, nroots, maxrys
      integer ig, jg, ij, i, j
      double precision zero, one
      double precision rln10, tol
      double precision rr, ai, arri, axi, ayi, azi
      double precision aj, aa, aa1, dum, fac
      double precision csi, cpi, cdi, cfi, cgi
      double precision csj, cpj, cdj, cfj, cgj
      double precision ax, ay, az, dum1, dum2
      double precision dij, xin, yin,zin, xxin, yyin, zzin
      double precision t1
      double precision dxsdi, dysdi, dzsdi, dxxdi, dyydi, dzzdi
      double precision dxsdj, dysdj, dzsdj, dxxdj, dyydj, dzzdj
      double precision dumxdx, dumydx, dumzdx
      double precision dumxdy, dumydy, dumzdy
      double precision dumxdz, dumydz, dumzdz
      double precision xint, yint, zint, xintx, yinty, zintz
      double precision t, xc, yc, zc, x0, y0, z0
      double precision xi, yi, zi, xj, yj, zj
c
      common/hnd_xyzdip/xint,yint,zint,xintx,yinty,zintz,t,xc,yc,zc,
     1              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension dij(225)
      dimension  xin(7,6), yin(7,6), zin(7,6)
      dimension xxin(7,6),yyin(7,6),zzin(7,6)
      dimension dxsdi (5,5),dysdi (5,5),dzsdi (5,5)
      dimension dxsdj (5,5),dysdj (5,5),dzsdj (5,5)
      dimension dxxdi (5,5),dyydi (5,5),dzzdi (5,5)
      dimension dxxdj (5,5),dyydj (5,5),dzzdj (5,5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension minprim(5), maxprim(5)
      data maxrys /5/
      data zero  /0.0d+00/
      data one   /1.0e+00/
      data rln10 /2.30258d+00/
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
      nroots=(lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_dipder: angular momentum too high',555)
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
        call errquit('hnd_dipder: illegal -lit- ', lit)
      endif
c
c     ----- j primitive -----
c
      do 6000 jg=j1,j2
      aj=expj(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000  ! the integral is zero
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
  360 dij(ij)=dum2
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
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      do 390 j=minj,maxj
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
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
      didij(ij,1) = didij(ij,1) + dij(ij)*dumxdx
      didij(ij,2) = didij(ij,2) + dij(ij)*dumydx
      didij(ij,3) = didij(ij,3) + dij(ij)*dumzdx
      didij(ij,4) = didij(ij,4) + dij(ij)*dumxdy
      didij(ij,5) = didij(ij,5) + dij(ij)*dumydy
      didij(ij,6) = didij(ij,6) + dij(ij)*dumzdy
      didij(ij,7) = didij(ij,7) + dij(ij)*dumxdz
      didij(ij,8) = didij(ij,8) + dij(ij)*dumydz
      didij(ij,9) = didij(ij,9) + dij(ij)*dumzdz
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
      djdij(ij,1) = djdij(ij,1) + dij(ij)*dumxdx
      djdij(ij,2) = djdij(ij,2) + dij(ij)*dumydx
      djdij(ij,3) = djdij(ij,3) + dij(ij)*dumzdx
      djdij(ij,4) = djdij(ij,4) + dij(ij)*dumxdy
      djdij(ij,5) = djdij(ij,5) + dij(ij)*dumydy
      djdij(ij,6) = djdij(ij,6) + dij(ij)*dumzdy
      djdij(ij,7) = djdij(ij,7) + dij(ij)*dumxdz
      djdij(ij,8) = djdij(ij,8) + dij(ij)*dumydz
      djdij(ij,9) = djdij(ij,9) + dij(ij)*dumzdz
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
