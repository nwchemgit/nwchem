      subroutine hnd_veloc(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,nint,veloc,scr,lscr)
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
      dimension veloc(nint,3)
c
c Use scratch for temps in integral calculation
c Scratch needs are 6*(Li+1)*(Lj+2)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_veloca(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,nint,veloc,
     2 scr(dip(1)),scr(dip(2)),scr(dip(3)),
     3 scr(dip(4)),scr(dip(5)),scr(dip(6)))
c
      return
      end
c
      subroutine hnd_veloca(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1                       expj,coefj,j_nprim,j_ngen,Lj,nint,veloc,
     2                       xs,ys,zs,dx,dy,dz)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,
     1                  xi,yi,zi,xj,yj,zj,ni,nj,cx,cy,cz
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension veloc(nint,3)
      dimension  xs(Li+1,Lj+2), ys(Li+1,Lj+2), zs(Li+1,Lj+2)
      dimension  dx(Li+1,Lj+2), dy(Li+1,Lj+2), dz(Li+1,Lj+2)
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*3,0.0d0,0,veloc,1)
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
c     ----- 
c
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do j=1,ljtmod
         nj=j
         do i=1,lit
            ni=i
            call hnd_sxyz
            xs(i,j)=xint*t
            ys(i,j)=yint*t
            zs(i,j)=zint*t
         enddo
      enddo
c
      do i=1,lit
         dx(i,1)=(-(aj+aj)*xs(i,2))
         dy(i,1)=(-(aj+aj)*ys(i,2))
         dz(i,1)=(-(aj+aj)*zs(i,2))
         do j=2,ljt
            dx(i,j)= (dble(j-1)*xs(i,j-1)-(aj+aj)*xs(i,j+1))
            dy(i,j)= (dble(j-1)*ys(i,j-1)-(aj+aj)*ys(i,j+1))
            dz(i,j)= (dble(j-1)*zs(i,j-1)-(aj+aj)*zs(i,j+1))
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
            dumx= dx(ix,jx)*ys(iy,jy)*zs(iz,jz)
            dumy= xs(ix,jx)*dy(iy,jy)*zs(iz,jz)
            dumz= xs(ix,jx)*ys(iy,jy)*dz(iz,jz)
            ij=ij+1
            veloc(ij,1)=veloc(ij,1)+dumx*cij
            veloc(ij,2)=veloc(ij,2)+dumy*cij
            veloc(ij,3)=veloc(ij,3)+dumz*cij
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
