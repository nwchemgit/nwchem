      subroutine hnd_giaol10(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,nint,l10,scr,lscr)
c
c $Id$
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      double precision l10(nint,3)
      dimension scr(lscr)
c
c Use scratch for temps in integral calculation
c Scratch needs are
c 9*(Li+1)*(Lj+2)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_giaol101(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,nint,l10,scr(gh01(1)),scr(gh01(2)),
     2 scr(gh01(3)),scr(gh01(4)),scr(gh01(5)),scr(gh01(6)),scr(gh01(7)),
     3 scr(gh01(8)),scr(gh01(9)))
c
      return
      end
c
      subroutine hnd_giaol101(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1                 expj,coefj,j_nprim,j_ngen,Lj,nint,l10,
     2                 xs,ys,zs,xsx,ysy,zsz,xsd,ysd,zsd)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_tol.fh"
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      double precision l10(nint,3)
      dimension  xs(Li+1,Lj+2), ys(Li+1,Lj+2), zs(Li+1,Lj+2)
      dimension xsx(Li+1,Lj+2),ysy(Li+1,Lj+2),zsz(Li+1,Lj+2)
      dimension xsd(Li+1,Lj+2),ysd(Li+1,Lj+2),zsd(Li+1,Lj+2)
      data rln10 /2.30258d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*3,0.0d0,0,l10,1)
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
c     ----- integral components -----
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
         do j=1,ljt
            do i=1,lit
               xsx(i,j)=xs(i,j+1)
               ysy(i,j)=ys(i,j+1)
               zsz(i,j)=zs(i,j+1)
            enddo
         enddo
C
         do i=1,lit
            xsd(i,1)= (-(aj+aj)*xs(i,2))
            ysd(i,1)= (-(aj+aj)*ys(i,2))
            zsd(i,1)= (-(aj+aj)*zs(i,2))
            do j=2,ljt
               xsd(i,j)= (dble(j-1)*xs(i,j-1)-(aj+aj)*xs(i,j+1))
               ysd(i,j)= (dble(j-1)*ys(i,j-1)-(aj+aj)*ys(i,j+1))
               zsd(i,j)= (dble(j-1)*zs(i,j-1)-(aj+aj)*zs(i,j+1))
            enddo
         enddo
c
c     ----- l integrals -----
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
            ij=ij+1
            dumx=  xs(ix,jx)*ysy(iy,jy)*zsd(iz,jz)
     1           - xs(ix,jx)*ysd(iy,jy)*zsz(iz,jz)
            dumy= xsd(ix,jx)* ys(iy,jy)*zsz(iz,jz)
     1           -xsx(ix,jx)* ys(iy,jy)*zsd(iz,jz)
            dumz= xsx(ix,jx)*ysd(iy,jy)* zs(iz,jz)
     1           -xsd(ix,jx)*ysy(iy,jy)* zs(iz,jz)
            l10(ij,1)=l10(ij,1)-cij*dumx
            l10(ij,2)=l10(ij,2)-cij*dumy
            l10(ij,3)=l10(ij,3)-cij*dumz
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
