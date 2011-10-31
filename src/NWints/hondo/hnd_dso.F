      subroutine hnd_dso(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj,j_nprim,j_ngen,Lj,nint,dso,scr,lscr,
     2  xyzpt, npairs)
c
c $Id$
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr) 
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension xyzpt(3,npairs,2),dso(nint,9,npairs)
c
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      call hnd_dso1(xyzi,expi,coefi,i_nprim,i_ngen, Li, xyzj,expj,
     1 coefj, j_nprim, j_ngen, Lj, nint, dso, xyzpt, npairs,
     2 scr(dsot(1)),scr(dsot(2)),scr(dsot(3)))
c
      return
      end
c
      subroutine hnd_dso1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1  expj,coefj, j_nprim, j_ngen, Lj, nint, dso, xyzpt, npairs,
     2  xv,yv,zv)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "hnd_wleg.fh"
#include "stdio.fh"
#include "hnd_tol.fh"
      common/hnd_xyzher/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  dx,dy,dz,ni,nj,cx,cy,cz,nc,nd
      dimension w2(maxrys)
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension xyzpt(3,npairs,2),dso(nint,9,npairs)
c
      dimension xv(Li+1,Lj+1,2,2,*)
      dimension yv(Li+1,Lj+1,2,2,*)
      dimension zv(Li+1,Lj+1,2,2,*)
c
      integer cg
c
c     Maximum scratch size needs to be Max(nder*3,1))*(Li+1)*(Lj+1)*((Li+Lj)/2+1)
c     Data format dso will be (nints,9,ipts)
c
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*9*npairs,0.0d0,0,dso,1)
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
c
      rrij=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- Gauss-Legendre center C is a p-function -----
c
c     ----- Number of Rys roots needed -----
c
      nroots=(lit+ljt+1+1-2)/2+1

      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_dso: need higher Rys root',nroots,INT_ERR)
      endif
c
c     ----- Loop over Gauss-Legendre center C -----
c
      do 8000 jpt=1,npairs
      dx=xyzpt(1,jpt,2)
      dy=xyzpt(2,jpt,2)
      dz=xyzpt(3,jpt,2)
      rrid=(xi-dx)**2+(yi-dy)**2+(zi-dz)**2
      rrjd=(xj-dx)**2+(yj-dy)**2+(zj-dz)**2
c
c     ----- i primitive -----
c
      do 7000 ig=1,i_nprim
      ai=expi(ig)
      csi=coefi(ig,i_ngen)
c
c     ----- j primitive -----
c
      do 6000 jg=1,j_nprim
      aj=expj(jg)
      csj=coefj(jg,j_ngen)
c
c     ----- Gauss-Legendre loop for one of the operators
c
      do 5000 cg=1,gleg_points
      ac=0.75d0*(1.0d0+hl(cg))/(1.0d0-hl(cg))
      ac=ac**2
      aa=ai+aj+ac
      aa1=1.0d0/aa
      dum=(ac*(ai*rrid+aj*rrjd)+ai*aj*rrij)*aa1
      if(dum.gt.tol) go to 5000
      fac=exp(-dum)
      ax=ai*xi+aj*xj+ac*dx
      ay=ai*yi+aj*yj+ac*dy
      az=ai*zi+aj*zj+ac*dz
c
c     Prefactor independent of electron position
c
      dij = fac * csi * csj
c
c     Prefactor from Gauss-Legendre Quadrature
c
c     alpha = 0.75*(1+t)/(1-t)
c     alphasquared= 0.75^2*(1+t)^2/(1-t)^2
c     wl    = 2*0.75/(1-t)^2
c     prefactor = 2*alphasquared*wl (* 2/Sqrt(pi))
c               = 4*0.75^3(1+t)^2/(1-t)^4
c
      dij = dij * wl(cg) / (1.0d0-hl(cg))**2 * ac * 3.0d0 * pi212
c
c     Prefactor from Rys Quadrature
c
      dij = dij * pi212 * aa1
c
      cx=xyzpt(1,jpt,1)
      cy=xyzpt(2,jpt,1)
      cz=xyzpt(3,jpt,1)
      yy=aa*((ax*aa1-cx)**2+(ay*aa1-cy)**2+(az*aa1-cz)**2)
      call hnd_droot
      do 420 iroot=1,nroots
         uu=u9(iroot)*aa
         ww=w9(iroot)
         w2(iroot)=ww*uu*2.0d0
         tt=1.0d0/(aa+uu)
         tx= sqrt(tt)
         x0=(ax+uu*cx)*tt
         y0=(ay+uu*cy)*tt
         z0=(az+uu*cz)*tt
         do 410 i=1,lit
            ni=i
            do 410 j=1,ljt
            nj=j
            do 410 nc=1,2
            do 410 nd=1,2
               call hnd_hermxyz
               xv(i,j,nc,nd,iroot)=xint
               yv(i,j,nc,nd,iroot)=yint
               zv(i,j,nc,nd,iroot)=zint
  410    continue
  420 continue
c
c        ----- combining the pieces together -----
c
      ij=0
      do 440 j=1,maxj
         call getNxyz(Lj,j,Nxyz)
         jx = Nxyz(1) + 1
         jy = Nxyz(2) + 1
         jz = Nxyz(3) + 1
         do 450 i=1,maxi
            call getNxyz(Li,i,Nxyz)
            ix = Nxyz(1) + 1
            iy = Nxyz(2) + 1
            iz = Nxyz(3) + 1
            ij=ij+1
            dumxx=0.0d0
            dumxy=0.0d0
            dumxz=0.0d0
            dumyx=0.0d0
            dumyy=0.0d0
            dumyz=0.0d0
            dumzx=0.0d0
            dumzy=0.0d0
            dumzz=0.0d0
            do 431 iroot=1,nroots
            dumxx=dumxx+xv(ix,jx,2,2,iroot)*yv(iy,jy,1,1,iroot)*
     1            zv(iz,jz,1,1,iroot)*w2(iroot)
            dumxy=dumxy+xv(ix,jx,2,1,iroot)*yv(iy,jy,1,2,iroot)*
     1            zv(iz,jz,1,1,iroot)*w2(iroot)
            dumxz=dumxz+xv(ix,jx,2,1,iroot)*yv(iy,jy,1,1,iroot)*
     1            zv(iz,jz,1,2,iroot)*w2(iroot)

            dumyx=dumyx+xv(ix,jx,1,2,iroot)*yv(iy,jy,2,1,iroot)*
     1            zv(iz,jz,1,1,iroot)*w2(iroot)
            dumyy=dumyy+xv(ix,jx,1,1,iroot)*yv(iy,jy,2,2,iroot)*
     1            zv(iz,jz,1,1,iroot)*w2(iroot)
            dumyz=dumyz+xv(ix,jx,1,1,iroot)*yv(iy,jy,2,1,iroot)*
     1            zv(iz,jz,1,2,iroot)*w2(iroot)

            dumzx=dumzx+xv(ix,jx,1,2,iroot)*yv(iy,jy,1,1,iroot)*
     1            zv(iz,jz,2,1,iroot)*w2(iroot)
            dumzy=dumzy+xv(ix,jx,1,1,iroot)*yv(iy,jy,1,2,iroot)*
     1            zv(iz,jz,2,1,iroot)*w2(iroot)
            dumzz=dumzz+xv(ix,jx,1,1,iroot)*yv(iy,jy,1,1,iroot)*
     1            zv(iz,jz,2,2,iroot)*w2(iroot)
  431       continue
            dso(ij,1,jpt) = dso(ij,1,jpt) + (dumyy+dumzz)*dij
            dso(ij,2,jpt) = dso(ij,2,jpt) - dumxy*dij
            dso(ij,3,jpt) = dso(ij,3,jpt) - dumxz*dij
            dso(ij,4,jpt) = dso(ij,4,jpt) - dumyx*dij
            dso(ij,5,jpt) = dso(ij,5,jpt) + (dumxx+dumzz)*dij
            dso(ij,6,jpt) = dso(ij,6,jpt) - dumyz*dij
            dso(ij,7,jpt) = dso(ij,7,jpt) - dumzx*dij
            dso(ij,8,jpt) = dso(ij,8,jpt) - dumzy*dij
            dso(ij,9,jpt) = dso(ij,9,jpt) + (dumxx+dumyy)*dij
  450    continue
  440 continue
c
 5000 continue
 6000 continue
 7000 continue
 8000 continue
c
      return
 9997 format(' in -hnd_dso- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,lct,nroots = ',4i3)
      end
