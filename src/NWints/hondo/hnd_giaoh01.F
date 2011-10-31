      subroutine hnd_giaoh01(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,nint,h01,xyzpt,nat,scr,lscr)
c
c $Id$
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr)
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension h01(nint,3,nat), xyzpt(3,nat)
c
c Use scratch for temps in integral calculation
c Scratch needs are
c 9*(Li+1)*(Lj+2)*((Li+Lj+2)/2+1)
c
c The pointers are defined in hdn_pointers.fh
c
      call hnd_giaoh011(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,nint,h01,xyzpt,nat,
     2 scr(gh01(1)),scr(gh01(2)),scr(gh01(3)),scr(gh01(4)),scr(gh01(5)),
     3 scr(gh01(6)),scr(gh01(7)),scr(gh01(8)),scr(gh01(9)))
c
      return
      end
c
      subroutine hnd_giaoh011(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1 expj,coefj,j_nprim,j_ngen,Lj,nint,h01,xyzpt,nat,
     2 xv,yv,zv,dxv,dyv,dzv,xvx,yvy,zvz)
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
      dimension h01(nint,3,nat), xyzpt(3,nat)
      dimension   xv(Li+1,Lj+2,*),  yv(Li+1,Lj+2,*),  zv(Li+1,Lj+2,*)
      dimension  dxv(Li+1,Lj+2,*), dyv(Li+1,Lj+2,*), dzv(Li+1,Lj+2,*)
      dimension  xvx(Li+1,Lj+2,*), yvy(Li+1,Lj+2,*), zvz(Li+1,Lj+2,*)
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*3*nat,0.0d0,0,h01,1)
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
      nroots=(lit+ljt)/2+1
      if(nroots.gt.maxrys) then
         write(luout,9999) maxrys,lit,ljt,nroots
         call errquit('hnd_giah01: need higher Rys root',nroots,INT_ERR)
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
      cij=csi*csj*fac*aa1*pi212
c
c     ----- -h01- term -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
c
      do ipt=1,nat
         znuc=1.0d0
         cx=xyzpt(1,ipt)
         cy=xyzpt(2,ipt)
         cz=xyzpt(3,ipt)
         yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
         call hnd_droot
         do iroot=1,nroots
            uu=u9(iroot)*aa
            ww=w9(iroot)*znuc
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
                  xv(i,j,iroot)=xint
                  yv(i,j,iroot)=yint
                  zv(i,j,iroot)=zint*ww
                  call hnd_dervxyz(1)
                  dxv(i,j,iroot)=xint
                  dyv(i,j,iroot)=yint
                  dzv(i,j,iroot)=zint*ww
               enddo
            enddo
            do i=1,lit
               xvx(i,1,iroot)=(-(aj+aj)*xv(i,2,iroot))
               yvy(i,1,iroot)=(-(aj+aj)*yv(i,2,iroot))
               zvz(i,1,iroot)=(-(aj+aj)*zv(i,2,iroot))
               do j=2,ljt
                  xvx(i,j,iroot)= (dble(j-1)*xv(i,j-1,iroot)-
     &                             (aj+aj)*xv(i,j+1,iroot))
                  yvy(i,j,iroot)= (dble(j-1)*yv(i,j-1,iroot)-
     &                             (aj+aj)*yv(i,j+1,iroot))
                  zvz(i,j,iroot)= (dble(j-1)*zv(i,j-1,iroot)-
     &                             (aj+aj)*zv(i,j+1,iroot))
               enddo
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
               dumx=0.0d0
               dumy=0.0d0
               dumz=0.0d0
               do ir=1,nroots
                  dumx=dumx+xv(ix,jx,ir)*dyv(iy,jy,ir)*zvz(iz,jz,ir)
     1                     -xv(ix,jx,ir)*yvy(iy,jy,ir)*dzv(iz,jz,ir)
                  dumy=dumy+xvx(ix,jx,ir)* yv(iy,jy,ir)*dzv(iz,jz,ir)
     1                     -dxv(ix,jx,ir)* yv(iy,jy,ir)*zvz(iz,jz,ir)
                  dumz=dumz+dxv(ix,jx,ir)*yvy(iy,jy,ir)* zv(iz,jz,ir)
     1                     -xvx(ix,jx,ir)*dyv(iy,jy,ir)* zv(iz,jz,ir)
               enddo     
               ij=ij+1
               h01(ij,1,ipt)=h01(ij,1,ipt)+dumx*cij
               h01(ij,2,ipt)=h01(ij,2,ipt)+dumy*cij
               h01(ij,3,ipt)=h01(ij,3,ipt)+dumz*cij
            enddo     
         enddo     
      enddo      
c
 1000 continue
c
      enddo
      enddo
c
      return
 9999 format(' in -giah01- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,' lit,ljt,nroots= ',3i3)
      end
