      subroutine hnd_giaoh11(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,scr,lscr, 
     & para,dia)
c
c $Id$
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr)
      logical para,dia
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension e11(nint,3,3,*)
      dimension catms(3,nat)
c
c Use scratch for temps in integral calculation
c Scratch needs are
c 11*3*(Li+2)*(Lj+2)*((Li+Lj+3)/2+1)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_giaoh111(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,para,dia,
     2 scr(gh01(1)),scr(gh01(2)),scr(gh01(3)),scr(gh01(4)),
     3 scr(gh01(5)),scr(gh01(6)),scr(gh01(7)),scr(gh01(8)),
     4 scr(gh01(9)),scr(gh01(10)),scr(gh01(11)))
c
      return
      end
c
      subroutine hnd_giaoh111(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1 expj,coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,para,dia,
     4 v0,dv0,v0d,v1,dv1,v11,dv11,v1d,v12,dv12,v11d)
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
      dimension e11(nint,3,3,*)
      dimension tm(21),catms(3,nat)
      dimension   v0(3,Li+2,Lj+2,*)    ! blocks for vx0,vy0,vz0
      dimension  dv0(3,Li+2,Lj+2,*)    ! blocks for dvx0,dvy0,dvz0
      dimension  v0d(3,Li+2,Lj+2,*)    ! blocks for vx0d,vy0d,vz0d
      dimension   v1(3,Li+2,Lj+2,*)    ! blocks for vx1,vy1,vz1
      dimension  dv1(3,Li+2,Lj+2,*)    ! blocks for dvx1,dvy1,dvz1
      dimension  v11(3,Li+2,Lj+2,*)    ! blocks for vx1x,vy1y,vz1z
      dimension dv11(3,Li+2,Lj+2,*)    ! blocks for dvx1x,dvy1y,dvz1z
      dimension  v1d(3,Li+2,Lj+2,*)    ! blocks for vx1d,vy1d,vz1d
      dimension  v12(3,Li+2,Lj+2,*)    ! blocks for xvx1,yvy1,zvz1
      dimension dv12(3,Li+2,Lj+2,*)    ! blocks for xdvx1,ydvy1,zdvz1
      dimension v11d(3,Li+2,Lj+2,*)    ! blocks for xvx1d,yvy1d,zvz1d
      logical para,dia
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*9*nat,0.0d0,0,e11,1)
c
c     ----- ishell -----
c
      xi=xyzi(1)
      yi=xyzi(2)
      zi=xyzi(3)
      lit = Li + 1
      maxi = lit*(lit+1)/2
      litmod=lit+1
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
      qijx=yi*zj-zi*yj
      qijy=zi*xj-xi*zj
      qijz=xi*yj-yi*xj
      tijx=xi   -   xj
      tijy=yi   -   yj
      tijz=zi   -   zj
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      nroots=(lit+ljt+3-2)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9999) maxrys,lit,ljt,nroots
         call errquit('hnd_giaoh11: need higher Rys rt',nroots,INT_ERR)
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
c     ----- -h11- dia + paramagnetic terms -----
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
               do i=1,litmod
                  ni=i
                  call hnd_sxyz
c
c     ----- for 1/r -----
c
                  if(kat.eq.1) then
                      v0(1,i,j,ir)=xint     ! vx0
                      v0(2,i,j,ir)=yint     ! vy0
                      v0(3,i,j,ir)=zint*vv  ! vz0
                  endif
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
            do j=1,ljt
               do i=1,lit
                    v11(1,i,j,ir)= v1(1,i,j+1,ir)   ! vx1x
                    v11(2,i,j,ir)= v1(2,i,j+1,ir)   ! vy1y
                    v11(3,i,j,ir)= v1(3,i,j+1,ir)   ! vz1z
                   dv11(1,i,j,ir)=dv1(1,i,j+1,ir)   ! dvx1x
                   dv11(2,i,j,ir)=dv1(2,i,j+1,ir)   ! dvy1y
                   dv11(3,i,j,ir)=dv1(3,i,j+1,ir)   ! dvz1z
               enddo ! i-loop
            enddo    ! j-loop
            do j=1,ljtmod
               do i=1,lit
                   v12(1,i,j,ir)= v1(1,i+1,j,ir)    ! xvx1
                   v12(2,i,j,ir)= v1(2,i+1,j,ir)    ! yvy1
                   v12(3,i,j,ir)= v1(3,i+1,j,ir)    ! zvz1
                  dv12(1,i,j,ir)=dv1(1,i+1,j,ir)    ! xdvx1
                  dv12(2,i,j,ir)=dv1(2,i+1,j,ir)    ! ydvy1
                  dv12(3,i,j,ir)=dv1(3,i+1,j,ir)    ! zdvz1
               enddo ! i-loop
            enddo    ! j-loop
         enddo       ! ir-loop
c
         do ir=1,nroots
c
c     ----- for 1/r derivatives wrt. centers -i- and -j- -----
c
         if (kat.eq.1) then
c
c     ----- derivatives with respect to xj ... -----
c
            do i=1,lit
               v0d(1,i,1,ir)=-(-(aj+aj)*v0(1,i,2,ir))          ! vx0d and vx0
               v0d(2,i,1,ir)=-(-(aj+aj)*v0(2,i,2,ir))          ! vy0d and vy0
               v0d(3,i,1,ir)=-(-(aj+aj)*v0(3,i,2,ir))          ! vz0d and vz0
               do j=2,ljt
                  v0d(1,i,j,ir)=-(dble(j-1)*v0(1,i,j-1,ir)-    ! vx0d and vx0
     &                              (aj+aj)*v0(1,i,j+1,ir))
                  v0d(2,i,j,ir)=-(dble(j-1)*v0(2,i,j-1,ir)-    ! vy0d and vy0
     &                              (aj+aj)*v0(2,i,j+1,ir))
                  v0d(3,i,j,ir)=-(dble(j-1)*v0(3,i,j-1,ir)-    ! vz0d and vz0
     &                              (aj+aj)*v0(3,i,j+1,ir))
               enddo ! j-loop
            enddo    ! i-loop
c
c     ----- derivatives with respect to xi ... -----
c
            do j=1,ljt
               dv0(1,1,j,ir)=-(-(ai+ai)*v0(1,2,j,ir))          ! dvx0 and vx0
               dv0(2,1,j,ir)=-(-(ai+ai)*v0(2,2,j,ir))          ! dvy0 and vy0
               dv0(3,1,j,ir)=-(-(ai+ai)*v0(3,2,j,ir))          ! dvz0 and vz0
               do i=2,lit
                  dv0(1,i,j,ir)=-(dble(i-1)*v0(1,i-1,j,ir)-    ! dvx0 and vx0
     &                              (ai+ai)*v0(1,i+1,j,ir))
                  dv0(2,i,j,ir)=-(dble(i-1)*v0(2,i-1,j,ir)-    ! dvy0 and vy0
     &                              (ai+ai)*v0(2,i+1,j,ir))
                  dv0(3,i,j,ir)=-(dble(i-1)*v0(3,i-1,j,ir)-    ! vz0d and vz0
     &                              (ai+ai)*v0(3,i+1,j,ir))
               enddo ! i-loop
            enddo    ! j-loop
         endif
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
         do i=1,lit
            v11d(1,i,1,ir)= (-(aj+aj)*v12(1,i,2,ir))           ! xvx1d and xvx1
            v11d(2,i,1,ir)= (-(aj+aj)*v12(2,i,2,ir))           ! yvy1d and yvy1
            v11d(3,i,1,ir)= (-(aj+aj)*v12(3,i,2,ir))           ! zvz1d and zvz1
            do j=2,ljt
               v11d(1,i,j,ir)= (dble(j-1)*v12(1,i,j-1,ir)-     ! xvx1d and xvx1
     &                            (aj+aj)*v12(1,i,j+1,ir))
               v11d(2,i,j,ir)= (dble(j-1)*v12(2,i,j-1,ir)-     ! yvy1d and yvy1
     &                            (aj+aj)*v12(2,i,j+1,ir))
               v11d(3,i,j,ir)= (dble(j-1)*v12(3,i,j-1,ir)-     ! zvz1d and zvz1
     &                            (aj+aj)*v12(3,i,j+1,ir))
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
               transx=0.0d0
               transy=0.0d0
               transz=0.0d0
               call dcopy(21,0.0d0,0,tm,1)
               do ir=1,nroots
c
                if(kat.eq.1) then
c
c     ----- translation invariance for nuclear derivatives of 1/r -----
c
                  transx=transx
     1                   +dv0(1,ix,jx,ir)*v0(2,iy,jy,ir)*v0(3,iz,jz,ir)
     2                   +dv1(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1(3,iz,jz,ir)
     3                   +v0d(1,ix,jx,ir)*v0(2,iy,jy,ir)*v0(3,iz,jz,ir)
                  transy=transy
     1                   +v0(1,ix,jx,ir)*dv0(2,iy,jy,ir)*v0(3,iz,jz,ir)
     2                   +v1(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1(3,iz,jz,ir)
     3                   +v0(1,ix,jx,ir)*v0d(2,iy,jy,ir)*v0(3,iz,jz,ir)
                  transz=transz
     1                   +v0(1,ix,jx,ir)*v0(2,iy,jy,ir)*dv0(3,iz,jz,ir)
     2                   +v1(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1(3,iz,jz,ir)
     3                   +v0(1,ix,jx,ir)*v0(2,iy,jy,ir)*v0d(3,iz,jz,ir)
                endif
c
c     ----- for h(1,1) diamagnetic -----
c
c       xx = xx + dvx1x *   vy1 *   vz1
c       xy = xy +  vx1x *  dvy1 *   vz1
c       xz = xz +  vx1x *   vy1 *  dvz1
c       yx = yx +  dvx1 *  vy1y *   vz1
c       yy = yy +   vx1 * dvy1y *   vz1
c       yz = yz +   vx1 *  vy1y *   vz1
c       zx = zx +  dvx1 *   vy1 *  dvz1
c       zy = zy +   vx1 *  dvy1 *  vz1z
c       zz = zz +   vx1 *   vy1 * dvz1z
c
                if (dia) then
              tm(1)=tm(1)+dv11(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1(3,iz,jz,ir)
              tm(2)=tm(2)+v11(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1(3,iz,jz,ir)
              tm(3)=tm(3)+v11(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1(3,iz,jz,ir)
              tm(4)=tm(4)+dv1(1,ix,jx,ir)*v11(2,iy,jy,ir)*v1(3,iz,jz,ir)
              tm(5)=tm(5)+v1(1,ix,jx,ir)*dv11(2,iy,jy,ir)*v1(3,iz,jz,ir)
              tm(6)=tm(6)+v1(1,ix,jx,ir)*v11(2,iy,jy,ir)*dv1(3,iz,jz,ir)
              tm(7)=tm(7)+dv1(1,ix,jx,ir)*v1(2,iy,jy,ir)*v11(3,iz,jz,ir)
              tm(8)=tm(8)+v1(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v11(3,iz,jz,ir)
              tm(9)=tm(9)+v1(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv11(3,iz,jz,ir)
                endif
c
c     ----- for h(1,1) paramagnetic -----
c
c               t10 = t10 +  vx1 * dvy1 * vz1d -  vx1 * vy1d * dvz1
c               t11 = t11 + vx1d *  vy1 * vz1d - dvx1 *  vy1 * vz1d
c               t12 = t12 + dvx1 * vy1d *  vz1 - vx1d * dvy1 *  vz1
c
                if (para) then
                tm(10)=tm(10)
     1                +v1(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
     2                -v1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*dv1(3,iz,jz,ir)
                tm(11)=tm(11)
     1                +v1d(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1(3,iz,jz,ir)
     2                -dv1(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
                tm(12)=tm(12)
     1                +dv1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*v1(3,iz,jz,ir)
     2                -v1d(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1(3,iz,jz,ir)
c
c       t13 = t10 +  xvx1 *  dvy1 *  vz1d -  xvx1 *  vy1d *  dvz1
c       t14 = t11 + xvx1d *   vy1 *  dvz1 - xdvx1 *   vy1 *  vz1d
c       t15 = t12 + xdvx1 *  vy1d *   vz1 - xvx1d *  dvy1 *   vz1
c       t16 = t12 +   vx1 * ydvy1 *  vz1d -   vx1 * yvy1d *  dvz1
c       t17 = t12 +  vx1d *  yvy1 *  dvz1 -  dvx1 *  yvy1 *  vz1d
c       t18 = t12 +  dvx1 * yvy1d *   vz1 -  vx1d * ydvy1 *   vz1
c       t19 = t12 +   vx1 *  dvy1 * zvz1d -   vx1 *  vy1d * zdvz1
c       t20 = t12 +  vx1d *   vy1 * zdvz1 -  dvx1 *   vy1 * zvz1d
c       t21 = t12 +  dvx1 *  vy1d *  zvz1 -  vx1d *  dvy1 *  zvz1
c
                tm(13)=tm(13)
     1                +v12(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
     2                -v12(1,ix,jx,ir)*v1d(2,iy,jy,ir)*dv1(3,iz,jz,ir)
                tm(14)=tm(14)
     1                +v11d(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1(3,iz,jz,ir)
     2                -dv12(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
                tm(15)=tm(15)
     1                +dv12(1,ix,jx,ir)*v1d(2,iy,jy,ir)*v1(3,iz,jz,ir)
     2                -v11d(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1(3,iz,jz,ir)
                tm(16)=tm(16)
     1                +v1(1,ix,jx,ir)*dv12(2,iy,jy,ir)*v1d(3,iz,jz,ir)
     2                -v1(1,ix,jx,ir)*v11d(2,iy,jy,ir)*dv1(3,iz,jz,ir)
                tm(17)=tm(17)
     1                +v1d(1,ix,jx,ir)*v12(2,iy,jy,ir)*dv1(3,iz,jz,ir)
     2                -dv1(1,ix,jx,ir)*v12(2,iy,jy,ir)*v1d(3,iz,jz,ir)
                tm(18)=tm(18)
     1                +dv1(1,ix,jx,ir)*v11d(2,iy,jy,ir)*v1(3,iz,jz,ir)
     2                -v1d(1,ix,jx,ir)*dv12(2,iy,jy,ir)*v1(3,iz,jz,ir)
                tm(19)=tm(19)
     1                +v1(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v11d(3,iz,jz,ir)
     2                -v1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*dv12(3,iz,jz,ir)
                tm(20)=tm(20)
     1                +v1d(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv12(3,iz,jz,ir)
     2                -dv1(1,ix,jx,ir)*v1(2,iy,jy,ir)*v11d(3,iz,jz,ir)
                tm(21)=tm(21)
     1                +dv1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*v12(3,iz,jz,ir)
     2                -v1d(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v12(3,iz,jz,ir)
                endif
c
               enddo     
               ij=ij+1
c
c     ----- h(1,1) diamagnetic -----
c
               if (dia) then
                  e11(ij,1,1,kat)=e11(ij,1,1,kat)+(tm(5)+tm(9))*cij
                  e11(ij,1,2,kat)=e11(ij,1,2,kat)-tm(2)*cij
                  e11(ij,1,3,kat)=e11(ij,1,3,kat)-tm(3)*cij
                  e11(ij,2,1,kat)=e11(ij,2,1,kat)-tm(4)*cij
                  e11(ij,2,2,kat)=e11(ij,2,2,kat)+(tm(1)+tm(9))*cij
                  e11(ij,2,3,kat)=e11(ij,2,3,kat)-tm(6)*cij
                  e11(ij,3,1,kat)=e11(ij,3,1,kat)-tm(7)*cij
                  e11(ij,3,2,kat)=e11(ij,3,2,kat)-tm(8)*cij
                  e11(ij,3,3,kat)=e11(ij,3,3,kat)+(tm(1)+tm(5))*cij
               endif
c
c     ----- h(1,1) paramagnetic -----
c
c     ----- a bit unclear ... yx ... instead of ... xy ... -----
c
               if (para) then
                  tm(1)=(qijx*tm(10)+tijy*tm(19)-tijz*tm(16))*cij
                  tm(2)=(qijx*tm(11)+tijy*tm(20)-tijz*tm(17))*cij 
                  tm(3)=(qijx*tm(12)+tijy*tm(21)-tijz*tm(18))*cij
                  tm(4)=(qijy*tm(10)+tijz*tm(13)-tijx*tm(19))*cij
                  tm(5)=(qijy*tm(11)+tijz*tm(14)-tijx*tm(20))*cij
                  tm(6)=(qijy*tm(12)+tijz*tm(15)-tijx*tm(21))*cij
                  tm(7)=(qijz*tm(10)+tijx*tm(16)-tijy*tm(13))*cij
                  tm(8)=(qijz*tm(11)+tijx*tm(17)-tijy*tm(14))*cij
                  tm(9)=(qijz*tm(12)+tijx*tm(18)-tijy*tm(15))*cij
c   
                  e11(ij,1,1,kat)=e11(ij,1,1,kat)+tm(1)
                  e11(ij,1,2,kat)=e11(ij,1,2,kat)+tm(4)
                  e11(ij,1,3,kat)=e11(ij,1,3,kat)+tm(7)
                  e11(ij,2,1,kat)=e11(ij,2,1,kat)+tm(2)
                  e11(ij,2,2,kat)=e11(ij,2,2,kat)+tm(5)
                  e11(ij,2,3,kat)=e11(ij,2,3,kat)+tm(8)
                  e11(ij,3,1,kat)=e11(ij,3,1,kat)+tm(3)
                  e11(ij,3,2,kat)=e11(ij,3,2,kat)+tm(6)
                  e11(ij,3,3,kat)=e11(ij,3,3,kat)+tm(9)
               endif
c
c              if(kat.eq.1) then
c                if(abs(transx*cij).ge.1.0d-08.or.abs(transy*cij)
c    &              .ge.1.0d-08.or.abs(transz*cij).ge.1.0d-08) then
c                  write(luout,9993)
c                  call errquit('hnd_giaoh11: no transl inv',0,INT_ERR)
c                endif
c              endif
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
 9999 format(' in -giaoh11- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,' lit,ljt,nroots= ',3i3)
 9993 format(' something wrong with translational',
     1       ' invariance in -giaoh11- . stop. ')
      end
