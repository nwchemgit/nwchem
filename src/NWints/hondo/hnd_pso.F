      subroutine hnd_pso(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,nint,e11,scr,lscr,catms,nat)
c
c $Id$
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr)
      logical para,dia
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension e11(nint,3,*)
      dimension catms(3,nat)
c
c Use scratch for temps in integral calculation
c Scratch needs are
c 3*(Li+1)*(Lj+2)*((Li+Lj+1)/2+1)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_pso1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,
     2 scr(dsot(1)),scr(dsot(2)),scr(dsot(3)))
c
      return
      end
c
      subroutine hnd_pso1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1 expj,coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,
     4 v1,dv1,v1d)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension e11(nint,3,*)
      dimension catms(3,nat)
      dimension   v1(3,Li+1,Lj+2,*)    ! blocks for vx1,vy1,vz1
      dimension  dv1(3,Li+1,Lj+2,*)    ! blocks for dvx1,dvy1,dvz1
      dimension  v1d(3,Li+1,Lj+2,*)    ! blocks for vx1d,vy1d,vz1d
      logical para,dia
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*3*nat,0.0d0,0,e11,1)
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      lit = Li + 1
      maxi = lit*(lit+1)/2
c
c     ----- jshell -----
c
      xj=xyzj(1)
      yj=xyzj(2)
      zj=xyzj(3)
      ljt = Lj + 1
      maxj = ljt*(ljt+1)/2
      ljtmod=ljt+1
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      nroots=(lit+ljt+2-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9999) maxrys,lit,ljt,nroots
         call errquit('hnd_pso: need higher Rys rt',nroots,INT_ERR)
      endif
c
c     ----- i primitive -----
c
      do ig=1, i_nprim
      ai=expi(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=coefi(ig,i_ngen)
c
c     ----- j primitive -----
c
      do jg=1,j_nprim
      aj=expj(jg)
      aa=ai+aj
      aa1=1.0d0/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) goto 1000
      fac= exp(-dum)
      csj=coefj(jg,j_ngen)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      cij=csi*csj*fac*pi212*aa1
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do kat=1,nat
         cx=catms(1,kat)
         cy=catms(2,kat)
         cz=catms(3,kat)
         xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
         yy=xx
         call hnd_droot
         do ir=1,nroots
            uu=u9(ir)*aa
            ww=w9(ir)
            vv=ww
            ww=ww*(uu+uu)
            tt=1.0d0/(aa+uu)
            t = sqrt(tt)
            x0=(aax+uu*cx)*tt
            y0=(aay+uu*cy)*tt
            z0=(aaz+uu*cz)*tt
            do j=1,ljtmod
               nj=j
               do i=1,lit
                  ni=i
                  call hnd_sxyz
c
c     ----- for x/r**3 -----
c
                  v1(1,i,j,ir)=xint         ! vx1
                  v1(2,i,j,ir)=yint         ! vy1
                  v1(3,i,j,ir)=zint*ww      ! vz1
                  call hnd_dervxyz(1)       
                  dv1(1,i,j,ir)=xint        ! dvx1
                  dv1(2,i,j,ir)=yint        ! dvy1
                  dv1(3,i,j,ir)=zint*ww     ! dvz1
               enddo ! i-loop
            enddo    ! j-loop
         enddo       ! ir-loop
c
         do ir=1,nroots
c
c     ----- d/dx ... operators -----
c
         do i=1,lit
            v1d(1,i,1,ir)= (-(aj+aj)*v1(1,i,2,ir))             ! vx1d and vx1
            v1d(2,i,1,ir)= (-(aj+aj)*v1(2,i,2,ir))             ! vy1d and vy1
            v1d(3,i,1,ir)= (-(aj+aj)*v1(3,i,2,ir))             ! vz1d and vz1
            do j=2,ljt
               v1d(1,i,j,ir)= (dble(j-1)*v1(1,i,j-1,ir)-       ! vx1d and vx1
     &                           (aj+aj)*v1(1,i,j+1,ir))
               v1d(2,i,j,ir)= (dble(j-1)*v1(2,i,j-1,ir)-       ! vy1d and vy1
     &                           (aj+aj)*v1(2,i,j+1,ir))
               v1d(3,i,j,ir)= (dble(j-1)*v1(3,i,j-1,ir)-       ! vz1d and vz1
     &                           (aj+aj)*v1(3,i,j+1,ir))
            enddo ! j-loop
         enddo    ! i-loop
c
         enddo    ! ir-loop
c
         ij=0
         do j=1,maxj
            call getNxyz(Lj,j,Nxyz)
            jx = Nxyz(1) + 1
            jy = Nxyz(2) + 1
            jz = Nxyz(3) + 1
            do i=1,maxi
               call getNxyz(Li,i,Nxyz)
               ix = Nxyz(1) + 1
               iy = Nxyz(2) + 1
               iz = Nxyz(3) + 1
               tmx = 0.0d0
               tmy = 0.0d0
               tmz = 0.0d0
               do ir=1,nroots
c
c     ----- for basic paramagnetic spin-orbit integrals -----
c
c               t10 = t10 +  vx1 * dvy1 * vz1d -  vx1 * vy1d * dvz1
c               t11 = t11 + vx1d *  vy1 * vz1d - dvx1 *  vy1 * vz1d
c               t12 = t12 + dvx1 * vy1d *  vz1 - vx1d * dvy1 *  vz1
c
                tmx = tmx
     1                +v1(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
     2                -v1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*dv1(3,iz,jz,ir)
                tmy = tmy
     1                +v1d(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1(3,iz,jz,ir)
     2                -dv1(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
                tmz = tmz
     1                +dv1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*v1(3,iz,jz,ir)
     2                -v1d(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1(3,iz,jz,ir)
c
               enddo     
               ij=ij+1
               e11(ij,1,kat)=e11(ij,1,kat)+tmx*cij
               e11(ij,2,kat)=e11(ij,2,kat)+tmy*cij
               e11(ij,3,kat)=e11(ij,3,kat)+tmz*cij
c
            enddo  ! j-loop final summation
         enddo     ! i-loop final summation
c
      enddo        ! kat-loop
c
 1000 continue
      enddo        ! jprim loop
      enddo        ! iprim loop
c
      return
 9999 format(' in -pso- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,' lit,ljt,nroots= ',3i3)
 9993 format(' something wrong with translational',
     1       ' invariance in -pso- . stop. ')
      end
