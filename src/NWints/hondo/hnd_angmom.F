      subroutine hnd_angmom(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,nint,angmom,scr,lscr)
c
c $Id:
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr)
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension angmom(nint,3)
c
c Use scratch for temps in integral calculation
c Scratch needs are 9*(Li+1)*(Lj+2)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_angmoma(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,nint,angmom,
     2 scr(dip(1)),scr(dip(2)),scr(dip(3)),scr(dip(4)),scr(dip(5)),
     3 scr(dip(6)),scr(dip(7)),scr(dip(8)),scr(dip(9)))
c
      return
      end
c
      subroutine hnd_angmoma(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1                       expj,coefj,j_nprim,j_ngen,Lj,nint,angmom,
     2                       xv,yv,zv,dx,dy,dz,xxv,yyv,zzv)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     &                  ni,nj,cx,cy,cz
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension angmom(nint,3)
      dimension  xv(Li+1,Lj+2),  yv(Li+1,Lj+2),  zv(Li+1,Lj+2)
      dimension  dx(Li+1,Lj+2),  dy(Li+1,Lj+2),  dz(Li+1,Lj+2)
      dimension xxv(Li+1,Lj+2), yyv(Li+1,Lj+2), zzv(Li+1,Lj+2)
c
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*3,0.0d0,0,angmom,1)
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
      cij=csi*csj*fac
c
c     ----- integrals, differentiation -----
c
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do i=1,lit
         ni=i
         do j=1,ljtmod
            nj=j
            call hnd_sxyz
            xv(i,j)=xint*t
            yv(i,j)=yint*t
            zv(i,j)=zint*t
         enddo
         do j=1,ljtmod
            xxv(i,j)=xv(i,j+1)
            yyv(i,j)=yv(i,j+1)
            zzv(i,j)=zv(i,j+1)
         enddo
      enddo
      do i=1,lit
         dx(i,1)=(-(aj+aj)*xv(i,2))
         dy(i,1)=(-(aj+aj)*yv(i,2))
         dz(i,1)=(-(aj+aj)*zv(i,2))
         do j=2,ljt
            dx(i,j)= (dble(j-1)*xv(i,j-1)-(aj+aj)*xv(i,j+1))
            dy(i,j)= (dble(j-1)*yv(i,j-1)-(aj+aj)*yv(i,j+1))
            dz(i,j)= (dble(j-1)*zv(i,j-1)-(aj+aj)*zv(i,j+1))
         enddo
      enddo
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
            dumx =  xv(ix,jx)*yyv(iy,jy)* dz(iz,jz)
     1            - xv(ix,jx)* dy(iy,jy)*zzv(iz,jz)
            dumy =  dx(ix,jx)* yv(iy,jy)*zzv(iz,jz)
     1            -xxv(ix,jx)* yv(iy,jy)* dz(iz,jz)
            dumz = xxv(ix,jx)* dy(iy,jy)* zv(iz,jz)
     1            - dx(ix,jx)*yyv(iy,jy)* zv(iz,jz)
            dumdx=  dx(ix,jx)* yv(iy,jy)* zv(iz,jz)
            dumdy=  xv(ix,jx)* dy(iy,jy)* zv(iz,jz)
            dumdz=  xv(ix,jx)* yv(iy,jy)* dz(iz,jz)
            ij=ij+1
            angmom(ij,1)=angmom(ij,1)-cij*(dumx+yj*dumdz-zj*dumdy)
            angmom(ij,2)=angmom(ij,2)-cij*(dumy+zj*dumdx-xj*dumdz)
            angmom(ij,3)=angmom(ij,3)-cij*(dumz+xj*dumdy-yj*dumdx)
         enddo     
      enddo     
c
 1000 continue
c
      enddo
      enddo
c
      return
      end
