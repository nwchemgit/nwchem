c
c $Id$
c
c Modified from HONDO
c
      subroutine hnd_vstat(xyzi,expi,coefi,i_nprim,i_ngen, Li, 
     1     xyzj,expj,coefj, j_nprim, j_ngen, Lj, 
     1     xyz,zan,nat,dens,s,ti,vi,lstv,doS,doT,doV,scr,lscr)
c
c This is a wrapper routine, setting up scratch blocks used in actual
c integral routine
c
      implicit none
c
      integer i_nprim           ! [input] num. prims on function i
      integer i_ngen            ! [input] num general conts on func. i
      integer Li                ! [input] angular momentum of func. i
      integer j_nprim           ! [input] num. prims on function j
      integer j_ngen            ! [input] num general conts on func. j
      integer Lj                ! [input] angular momentum of func. j
      integer nat               ! [input] number of atoms
      integer lscr              ! [input] size of scratch array
      integer lstv              ! [input] size of any integral buffer
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim) ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen) ! [input] j coeffs
      double precision xyz(3,nat) ! [input] all atom positions
      double precision zan(nat) ! [input] charges on all atoms
      double precision dens(*)  ! [input] (transition)density 
      double precision scr(lscr) ! [scratch] scratch buffers
      double precision s(lstv)  ! [output] overlap integrals
      double precision ti(lstv) ! [output] kinetic energy integrals
      double precision vi(lstv) ! [output] potential integrals
      logical doS               ! [input] compute overlap (True/False)
      logical doT               ! [input] compute kinetic (True/False)
      logical doV               ! [input] compute potential (True/False)
c
      integer ss,st,tt,tv,vv, dd
c
c Use scratch for temps in integral calculation
c Scratch needs are dependent on integral types:
c doS : 3*(Li+1)*(Lj+3)
c doT : 3*(Li+1)*(Lj+3)+3*(Li+1)*(Lj+1)
c doV : 3*(Li+1)*(Lj+1)*((Li+Lj)/2+1)
c
      ss = 0
      tt = 0
      vv = 0
      st = 1
      tv = 1
      dd = (Li+1)*(Li+2)*(Lj+1)*(Lj+2)/4
      if (doS.or.doT) then
         ss = dd + (Li+1) * (Lj+3)
         st = st + 3 * ss
         if (doV) tv = st
         if (doT) then
            tt = (Li+1) * (Lj+1)
            if (doV) tv = tv + 3 * tt
         endif
      endif
      if (doV) vv = (Li+1) * (Lj+1) * ((Li+Lj)/2+1)
c
      call hnd_vstat1(xyzi,expi,coefi,i_nprim,i_ngen, Li, 
     1     xyzj,expj,coefj, j_nprim, j_ngen, Lj, 
     1     xyz,zan,nat,dens,s,ti,vi,lstv,doS,doT,doV,scr(1),
     3     scr(1+dd),scr(1+ss),scr(1+2*ss),scr(st),scr(st+tt),
     4     scr(st+2*tt),scr(tv),scr(tv+vv),scr(tv+2*vv))
c
      return
      end
c
      subroutine hnd_vstat1(xyzi,expi,coefi,i_nprim,i_ngen, Li, 
     1     xyzj,expj,coefj, j_nprim, j_ngen, Lj, 
     1     xyz,zan,nat,dens,s,ti,vi,lstv,doS,doT,doV,dij,
     1     xs,ys,zs,xt,yt,zt,xv,yv,zv)
c     
      implicit none
#include "stdio.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "errquit.fh"
      integer i_nprim           ! [input] num. prims on function i
      integer i_ngen            ! [input] num general conts on func. i
      integer Li                ! [input] angular momentum of func. i
      integer j_nprim           ! [input] num. prims on function j
      integer j_ngen            ! [input] num general conts on func. j
      integer Lj                ! [input] angular momentum of func. j
      integer nat               ! [input] number of atoms
      integer lscr              ! [input] size of scratch array
      integer lstv              ! [input] size of any integral buffer
      double precision xyzi(3)  ! [input] position of center i
      double precision expi(i_nprim) ! [input] exponents on i
      double precision coefi(i_nprim,i_ngen) ! [input] i coeffs
      double precision xyzj(3)  ! [input] position of center j
      double precision expj(j_nprim) ! [input] exponents on j
      double precision coefj(j_nprim,j_ngen) ! [input] j coeffs
      double precision xyz(3,nat) ! [input] all atom positions
      double precision zan(nat) ! [input] charges on all atoms
      double precision dens(*)  ! [input] (transition)density 
      double precision s(lstv)  ! [output] overlap integrals
      double precision ti(lstv) ! [output] kinetic energy integrals
      double precision vi(lstv) ! [output] potential integrals
      logical doS               ! [input] compute overlap (True/False)
      logical doT               ! [input] compute kinetic (True/False)
      logical doV               ! [input] compute potential (True/False)
      double precision tol, aa, aa1, rr, ai, arri, axi, ayi, azi
      double precision csi, cpi, cdi, cfi, cgi, aj, fac
      double precision csj, cpj, cdj, cfj, cgj, ax, ay, az
      double precision xs, ys, zs, xt, yt, zt, dum, dum1, dum2
      double precision xint, yint, zint, t, x0, y0, z0, xi, yi, zi 
      double precision xj, yj, zj, ccx, ccy, ccz 
      double precision rln10, zero, one, pi212
      double precision aax, aay, aaz, znuc, cx, cy, cz, uu, ww, tt
      double precision xv, yv, zv, dij
      integer lit, ljt, maxi, maxj, ljtmod, jg
      integer ijx, ijy, ijz, ij, i, j
      integer ni, nj, ig, ic, iroot
      integer ix, iy, iz, jx, jy, jz, Nxyz(3)
      logical some
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1     ni,nj,ccx,ccy,ccz
      dimension xs(Li+1,Lj+3)  ,ys(Li+1,Lj+3)  ,zs(Li+1,Lj+3)
      dimension xt(Li+1,Lj+1)  ,yt(Li+1,Lj+1)  ,zt(Li+1,Lj+1)
      dimension xv(Li+1,Lj+1,*),yv(Li+1,Lj+1,*),zv(Li+1,Lj+1,*)
      dimension dij(*)
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pi212  /1.1283791670955d+00/
c    
c      write(*,*)"basis function information" 
c      write(*,*)"xyz=", xyz, xyzi,xyzj
c      write(*,*)"expi,coefi,i_nprim,i_ngen, Li", 
c     &     expi, coefi,i_nprim,i_ngen, Li
c      write(*,*)"expj,coefj,j_nprim,j_ngen, Lj", 
c     &     expj, coefj,j_nprim,j_ngen, Lj
      
c      write(*,*)"vints in hnd_vstat, before"
c      write(*,"(4f12.8)")((xyz(ig,jg),ig=1,3),vi(jg),jg=1,nat)
c
      tol=rln10*itol
c     
c     ----- calculate -s-, -t-, and -v- integrals -----
c     
      some = .false.
      if(some) write(luout,9999)
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
      ljtmod=ljt+2
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_stvint: angular momentum too high',555,
     &       INT_ERR)
      endif
c
      if (doS) call dcopy(lstv,zero,0,s,1)
      if (doT) call dcopy(lstv,zero,0,ti,1)
c
c      if (doV) call dcopy(lstv,zero,0,vi,1)
c
c     do not fill vi with zero's since the contributions are summed over 
c     all basis pairs
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
            if(dum.gt.tol) go to 6000 ! the integral is zero
            fac= exp(-dum)
            csj=coefj(jg,j_ngen)
            ax=(axi+aj*xj)*aa1
            ay=(ayi+aj*yj)*aa1
            az=(azi+aj*zj)*aa1
c     
c     ----- density factor -----
c     
            ij=0
            dum2=fac*csi*csj
            do 3600 i=1,maxi
               do 360 j=1,maxj
                  ij=ij+1
                  dij(ij)=dum2*dens(ij)
 360           continue
 3600       continue
c     
c     ----- overlap and kinetic energy -----
c     
            if (doS.or.doT) then
            t = sqrt(aa1)
            x0=ax
            y0=ay
            z0=az
            do 3700 j=1,ljtmod
               nj=j
               do 370 i=1,lit
                  ni=i
                  call hnd_sxyz
                  xs(i,j)=xint*t
                  ys(i,j)=yint*t
                  zs(i,j)=zint*t
 370           continue
 3700       continue
c
            if (doT) call hnd_txyz(xt,yt,zt,xs,ys,zs,lit,ljt,aj)
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
                  dum =xs(ix,jx)*ys(iy,jy)*zs(iz,jz)
                  if (doS) s(ij)= s(ij)+ dum*dij(ij)
                  if (doT) then
                     dum =xt(ix,jx)*ys(iy,jy)*zs(iz,jz)
     1                    +xs(ix,jx)*yt(iy,jy)*zs(iz,jz)
     2                    +xs(ix,jx)*ys(iy,jy)*zt(iz,jz)
                     ti(ij)=ti(ij)+(dum*dij(ij))
                  endif
 380           continue
 390        continue
            endif  ! of overlap and kinetic energy integral
c     
c     ----- nuclear attraction -----
c     
            if (doV) then
               aax=aa*ax
               aay=aa*ay
               aaz=aa*az
               do 500 ic=1,nat
                  znuc=-zan(ic)
                  cx=xyz(1,ic)
                  cy=xyz(2,ic)
                  cz=xyz(3,ic)
                  yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                  call hnd_droot
                  do 420 iroot=1,nroots
                     uu=u9(iroot)*aa
                     ww=w9(iroot)*znuc
                     tt=one/(aa+uu)
                     t = sqrt(tt)
                     x0=(aax+uu*cx)*tt
                     y0=(aay+uu*cy)*tt
                     z0=(aaz+uu*cz)*tt
                     do 4100 j=1,ljt
                        nj=j
                        do 410 i=1,lit
                           ni=i
                           call hnd_sxyz
                           xv(i,j,iroot)=xint
                           yv(i,j,iroot)=yint
                           zv(i,j,iroot)=zint*ww
 410                    continue
 4100                continue
 420              continue
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
                        dum=zero
                        do 430 iroot=1,nroots
c                           write(*,*)"ic,ix,iy,iz,jx,jy,jz,xv,yy,zv",
c     &                          ic,ix,iy,iz,jx,jy,jz,
c     &                          xv(ix,jx,iroot),yv(iy,jy,iroot),
c     &                          zv(iz,jz,iroot)
                           dum=dum+xv(ix,jx,iroot)
     &                          *yv(iy,jy,iroot)*zv(iz,jz,iroot)
 430                    continue
                        dum=dum*(aa1*pi212)
                        ij=ij+1
c                        write(*,"('den(ij)=', f10.4,'   ai=',f10.4,
c     &                       '   aj=',f10.4
c     &                       '   ci=',f10.4,'   aj=',f10.4,'   soo=',
c     &                       f10.4)")
c     &                       dens(ij),ai,aj,dum1/fac,dum2/dum1,dij(ij)
c                        write(*,"('i=',i5,'   j=',i5,'   ij=',i5,
c     &                       '   dum=',f12.8,'   dij=',f12.8,
c     &                       '   dij*dum=',f12.8)")
c     &                       i,j,ij,dum,dij(ij),dum*dij(ij)
                        vi(ic)=vi(ic)+dum*dij(ij)
 440                 continue
 450              continue
c
 500           continue
            endif               ! of nuclear part
c     
 6000    continue
 7000 continue
c
      if(some) write(luout,9998)
      return
 9999 format(/,10x,20(1h-),/,10x,'1 electron integrals',
     2       /,10x,20(1h-))
 9998 format(' ...... end of one-electron integrals ......')
 9997 format(' in -hnd_stvint- the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i2,/,
     2       ' lit,ljt,nroots= ',3i3)
      end






