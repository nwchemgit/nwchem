c
c $Id: hnd_vstat.F,v 1.1 2002-01-12 06:47:30 zzhang Exp $
c
c Taken from HONDO
c
      subroutine hnd_vstat(xyzi,expi,coefi,i_nprim,i_ngen, Li, 
     1     xyzj,expj,coefj, j_nprim, j_ngen, Lj, 
     1     xyz,zan,nat,dens,s,ti,vi,lstv,doS,doT,doV,dij,lscr)
c     
      implicit none
#include "stdio.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
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
      double precision dij(lscr) ! [scratch] scratch buffers
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
      double precision xv, yv, zv
      double precision vzhang(6,6,2)
      integer lit, ljt, mini, minj, maxi, maxj, ljtmod, jg
      integer ijx, ijy, ijz, minprim, maxprim, ij, i, j
      integer ni, nj, maxrys, i1, i2, j1, j2, ig, ic, iroot
      integer ix, iy, iz, jx, jy, jz, idim, jdim
      logical some
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1     ni,nj,ccx,ccy,ccz
      dimension ijx(35),ijy(35),ijz(35)
      dimension xs(5,7)  ,ys(5,7)  ,zs(5,7)
      dimension xt(5,5)  ,yt(5,5)  ,zt(5,5)
      dimension xv(5,5,5),yv(5,5,5),zv(5,5,5)
      dimension minprim(5), maxprim(5)
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pi212  /1.1283791670955d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1     4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2     5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3     3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1     1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2     1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3     1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1     1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2     1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3     3, 3, 2, 2, 3/
      data minprim / 1, 2, 5,11,21/
      data maxprim / 1, 4,10,20,35/
     
c      write(*,*)"basis function information" 
c      write(*,*)"xyz=", xyz, xyzi,xyzj
c      write(*,*)"expi,coefi,i_nprim,i_ngen, Li", 
c     &     expi, coefi,i_nprim,i_ngen, Li
c      write(*,*)"expj,coefj,j_nprim,j_ngen, Lj", 
c     &     expj, coefj,j_nprim,j_ngen, Lj
c      call dfill(72,0d0,vzhang,1)
      
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
      i1=1
      i2=i_nprim
      lit = Li + 1
      mini=minprim(lit)
      maxi=maxprim(lit)
      idim=maxi-mini+1
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
      jdim=maxj-minj+1
      ljtmod=ljt+2
c
c     Transform the density from nw order to hnd order
c
      if ((lit.gt.2).or.(ljt.gt.2)) then
         if (doV) then
            call nw_to_hnd_order(dij,dens,minj,mini,maxj-minj+1,
     1           maxi-mini+1)
            call dcopy (idim*jdim,dij,1,dens,1)
         endif
      endif
c
c      write(*,*)"dens in hdn"
c      call output(dens, 1, jdim, 1, idim, jdim,
c     &              idim, 1)

      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_stvint: angular momentum too high',555)
      endif
c
      if (doS) call dfill(lstv,zero,s,1)
      if (doT) call dfill(lstv,zero,ti,1)
c
c      if (doV) call dfill(lstv,zero,vi,1)
c
c     do not fill vi with zero's since the contributions are summed over 
c     all basis pairs
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
            if(dum.gt.tol) go to 6000 ! the integral is zero
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
            do 3600 i=mini,maxi
               go to (110,120,220,220,130,220,220,220,220,220,
     1              150,220,220,220,220,220,220,220,220,220,
     2              180,220,220,220,220,220,220,220,220,220,
     3              220,220,220,220,220),i
 110           dum1=csi*fac
               go to 220
 120           dum1=cpi*fac
               go to 220
 130           dum1=cdi*fac
               go to 220
 150           dum1=cfi*fac
               go to 220
 180           dum1=cgi*fac
 220           continue
c     
               do 360 j=minj,maxj
                  go to (230,250,350,350,260,350,350,350,350,350,
     1                 280,350,350,350,350,350,350,350,350,350,
     2                 310,350,350,350,350,350,350,350,350,350,
     3                 350,350,350,350,350),j
 230              dum2=dum1*csj
                  go to 350
 250              dum2=dum1*cpj
                  go to 350
 260              dum2=dum1*cdj
                  go to 350
 280              dum2=dum1*cfj
                  go to 350
 310              dum2=dum1*cgj
 350              continue
c     
                  ij=ij+1
c
c
c     dij(ij)=dum2
c     dens(ij)=1
                  dij(ij)=dum2*dens(ij)
 360           continue
 3600       continue
c     
c     ----- overlap and kinetic energy -----
c     
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
            call hnd_txyz(xt,yt,zt,xs,ys,zs,lit,ljt,aj)
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
                  xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                  if(nroots.le.3) call hnd_rt123
                  if(nroots.eq.4) call hnd_root4
                  if(nroots.eq.5) call hnd_root5
                  do 420 iroot=1,nroots
                     uu=u(iroot)*aa
                     ww=w(iroot)*znuc
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
                  do 450 i=mini,maxi
                     ix=ijx(i)
                     iy=ijy(i)
                     iz=ijz(i)
                     do 440 j=minj,maxj
                        jx=ijx(j)
                        jy=ijy(j)
                        jz=ijz(j)
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
c                        vzhang(i-mini+1,j-minj+1,ic)=
c     &                       vzhang(i-mini+1,j-minj+1,ic)+dum*dij(ij)
 440                 continue
 450              continue
c
 500           continue
            endif               ! of nuclear part
c     
 6000    continue
 7000 continue
      do 550 i=mini,maxi
         ix=ijx(i)
         iy=ijy(i)
         iz=ijz(i)
         do 540 j=minj,maxj
            jx=ijx(j)
            jy=ijy(j)
            jz=ijz(j)
c            write(*,"(6i1,2f12.5)")ix-1,iy-1,iz-1,jx-1,jy-1,jz-1,
c     &           vzhang(i-mini+1,j-minj+1,1),
c     &           vzhang(i-mini+1,j-minj+1,2)
 540     continue
 550  continue
c      write(*,*)"vints in hnd_vstat,after"
c      write(*,"(4f12.8)")((xyz(ig,jg),ig=1,3),vi(jg),jg=1,nat)
c
      goto 101
c     
c     ----- perform transformation to NWChem order if necessary ----
c     
      if ((lit.gt.2).or.(ljt.gt.2)) then
         if (doS) then
            call hnd_to_nw_order(dij,s,minj,mini,maxj-minj+1,
     1           maxi-mini+1)
            call dcopy (lstv,dij,1,s,1)
         endif
         if (doT) then
            call hnd_to_nw_order(dij,ti,minj,mini,maxj-minj+1,
     1           maxi-mini+1)
            call dcopy (lstv,dij,1,ti,1)
         endif
         if (doV) then
            call hnd_to_nw_order(dij,vi,minj,mini,maxj-minj+1,
     1           maxi-mini+1)
            call dcopy (lstv,dij,1,vi,1)
         endif
      endif
c
 101  continue
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






