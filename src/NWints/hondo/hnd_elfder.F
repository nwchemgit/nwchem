      subroutine hnd_elfder(xyzi,expi,coefi,i_nprim,i_ngen, Li,
     1  xyzj,expj,coefj,j_nprim,j_ngen,Lj,nder,nint,elfder,scr,lscr,
     2  xyzpt, npt)
c
c $Id$
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension xyzpt(3,npt)
      dimension scr(lscr) 
      dimension elfder(*)
c
c
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      idim = max(nder*3,1)
c
      call hnd_elfder1(xyzi,expi,coefi,i_nprim,i_ngen, Li, xyzj,expj,
     1 coefj, j_nprim, j_ngen, Lj, nder, nint, elfder, xyzpt, npt, idim,
     2 scr(elpt(1)),scr(elpt(2)),scr(elpt(3)),scr(elpt(4)),scr(elpt(5)),
     3 scr(elpt(6)),scr(elpt(7)),scr(elpt(8)),scr(elpt(9)))
c
      return
      end
c
      subroutine hnd_elfder1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1  expj,coefj, j_nprim, j_ngen, Lj, nder, nint, elfder, xyzpt, npt,
     2  idim,xv,yv,zv,dxv,dyv,dzv,ddxv,ddyv,ddzv)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "stdio.fh"
#include "hnd_tol.fh"
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
      dimension w2(maxrys),w4(maxrys)
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension xyzpt(3,npt),elfder(nint,idim,npt)
c
      dimension xv(Li+1,Lj+1,*)
      dimension yv(Li+1,Lj+1,*)
      dimension zv(Li+1,Lj+1,*)
      dimension dxv(Li+1,Lj+1,*)
      dimension dyv(Li+1,Lj+1,*)
      dimension dzv(Li+1,Lj+1,*)
      dimension ddxv(Li+1,Lj+1,*)
      dimension ddyv(Li+1,Lj+1,*)
      dimension ddzv(Li+1,Lj+1,*)
c
c     Routine calculates the following integrals for a given shell:
c
c     nder =-2 : Electronic wave function
c     nder =-1 : Electronic density
c     nder = 0 : Electrostatic potential
c     nder = 1 : Electric field
c     nder = 2 : Electric field gradient
c
c     Maximum scratch size needs to be Max(nder*3,1))*(Li+1)*(Lj+1)*((Li+Lj)/2+1)
c     Data format elfder will be (nints,ipts,Max(nder*3,1))
c
      data rln10  /2.30258d+00/
      data pi212  /1.1283791670955d+00/
c
      tol=rln10*itol
c
c     Zero integral array
c
      call dcopy(nint*max(nder*3,1)*npt,0.0d0,0,elfder,1)
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
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1

      if(nroots.gt.maxrys) then
         write(luout,9997) maxrys,lit,ljt,nroots
         call errquit('hnd_elfgrd: need higher Rys root',nroots,INT_ERR)
      endif
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
      aa1=1.0d0/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=coefj(jg,j_ngen)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      dum1 = csi * fac
      dij = dum1 * csj
c
c     Electronic density integrals and electronic wave
c     function are done differently
c     then the electric field, electric field gradient, and
c     electrostatic potential integrals
c     Switch between the two classes here:
c
      if (nder.ge.0) goto 399
c
c     ----- electronic wave function / density integral -----
c
c     For wave function we simply set j to 0.0d0 and reset the
c     density factor appropriately
c
      if (nder.eq.-2) then
         dij = csi
         ax = xi
         ay = yi
         az = zi
         aa = ai
         ljt = 1
         maxj = 1
      endif
c
      do 300 ipt=1,npt
         x0=xyzpt(1,ipt)
         y0=xyzpt(2,ipt)
         z0=xyzpt(3,ipt)
         dum = aa*((x0-ax)**2+(y0-ay)**2+(z0-az)**2)
         if(dum.gt.tol) go to 300
         fac = exp(-dum)
c
c        ----- density values -----
c
         do 270 j=1,ljt
            nj=j
            do 270 i=1,lit
               ni=i
               call hnd_denxyz
               xv(i,j,1)=xint
               yv(i,j,1)=yint
               zv(i,j,1)=zint
  270    continue
c
c        ----- combining the pieces together -----
c
         ij=0
         do 290 j=1,maxj
            call getNxyz(Lj,j,Nxyz)
            jx = Nxyz(1) + 1
            jy = Nxyz(2) + 1
            jz = Nxyz(3) + 1
            do 290 i=1,maxi
               call getNxyz(Li,i,Nxyz)
               ix = Nxyz(1) + 1
               iy = Nxyz(2) + 1
               iz = Nxyz(3) + 1
               ij=ij+1
               elfder(ij,1,ipt)=elfder(ij,1,ipt)+fac*dij*xv(ix,jx,1)*
     1                          yv(iy,jy,1)*zv(iz,jz,1)
  290    continue
  300 continue
c
      goto 6000
c
c     ----- electric field (gradient) term -----
c
  399 dij = dij * pi212 * aa1
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ipt=1,npt
         znuc=1.0d0
         cx=xyzpt(1,ipt)
         cy=xyzpt(2,ipt)
         cz=xyzpt(3,ipt)
         yy=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
         call hnd_droot
!DEC$ LOOP COUNT MAX=20, MIN=6
         do 420 iroot=1,nroots
            uu=u9(iroot)*aa
            u4=uu*uu
            ww=w9(iroot)*znuc
            w2(iroot)=ww*uu*2.0d0
            w4(iroot)=ww*u4*4.0d0
            w9(iroot)=ww
            tt=1.0d0/(aa+uu)
            tx= sqrt(tt)
            x0=(aax+uu*cx)*tt
            y0=(aay+uu*cy)*tt
            z0=(aaz+uu*cz)*tt
            do 410 i=1,lit
               ni=i
               do 410 j=1,ljt
                  nj=j
                  goto (402,401) nder+1
                  call hnd_dervxyz(2)
                  ddxv(i,j,iroot)=xint
                  ddyv(i,j,iroot)=yint
                  ddzv(i,j,iroot)=zint
  401             call hnd_dervxyz(1)
                  dxv(i,j,iroot)=xint
                  dyv(i,j,iroot)=yint
                  dzv(i,j,iroot)=zint
  402             call hnd_sxyz
                  xv(i,j,iroot)=xint
                  yv(i,j,iroot)=yint
                  zv(i,j,iroot)=zint
  410       continue
  420    continue
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
            if (nder.eq.2) then
               dumxx=0.0d0
               dumyy=0.0d0
               dumzz=0.0d0
               dumxy=0.0d0
               dumxz=0.0d0
               dumyz=0.0d0
!DEC$ LOOP COUNT MAX=20, MIN=6
               do 430 iroot=1,nroots
                  dum=xv(ix,jx,iroot)*yv(iy,jy,iroot)*
     1                zv(iz,jz,iroot)*w2(iroot)
                  dumxx=dumxx - dum + ddxv(ix,jx,iroot)*  
     1                  yv(iy,jy,iroot)*zv(iz,jz,iroot)*w4(iroot)
                  dumyy=dumyy - dum + ddyv(iy,jy,iroot)*
     1                  xv(ix,jx,iroot)*zv(iz,jz,iroot)*w4(iroot)
                  dumzz=dumzz - dum + ddzv(iz,jz,iroot)* 
     1                  xv(ix,jx,iroot)*yv(iy,jy,iroot)*w4(iroot)
                  dumxy=dumxy + dxv(ix,jx,iroot)*dyv(iy,jy,iroot)*
     1                  zv(iz,jz,iroot)*w4(iroot)
                  dumxz=dumxz + dxv(ix,jx,iroot)*dzv(iz,jz,iroot)*
     1                  yv(iy,jy,iroot)*w4(iroot)
                  dumyz=dumyz + dyv(iy,jy,iroot)*dzv(iz,jz,iroot)*
     1                  xv(ix,jx,iroot)*w4(iroot)
  430          continue
               
               elfder(ij,1,ipt) = elfder(ij,1,ipt) + dumxx*dij 
               elfder(ij,2,ipt) = elfder(ij,2,ipt) + dumyy*dij
               elfder(ij,3,ipt) = elfder(ij,3,ipt) + dumzz*dij
               elfder(ij,4,ipt) = elfder(ij,4,ipt) + dumxy*dij
               elfder(ij,5,ipt) = elfder(ij,5,ipt) + dumxz*dij
               elfder(ij,6,ipt) = elfder(ij,6,ipt) + dumyz*dij
            elseif (nder.eq.1) then
               dumx=0.0d0
               dumy=0.0d0
               dumz=0.0d0
!DEC$ LOOP COUNT MAX=20, MIN=6
               do 431 iroot=1,nroots
                  dumx=dumx + dxv(ix,jx,iroot)*yv(iy,jy,iroot)*  
     1                 zv(iz,jz,iroot)*w2(iroot)
                  dumy=dumy + xv(ix,jx,iroot)*dyv(iy,jy,iroot)*  
     1                 zv(iz,jz,iroot)*w2(iroot)
                  dumz=dumz + xv(ix,jx,iroot)*yv(iy,jy,iroot)* 
     1                 dzv(iz,jz,iroot)*w2(iroot)
  431          continue
               elfder(ij,1,ipt) = elfder(ij,1,ipt) + dumx*dij 
               elfder(ij,2,ipt) = elfder(ij,2,ipt) + dumy*dij
               elfder(ij,3,ipt) = elfder(ij,3,ipt) + dumz*dij
            else
               dumx=0.0d0
!DEC$ LOOP COUNT MAX=20, MIN=6
               do 432 iroot=1,nroots
                  dumx=dumx + xv(ix,jx,iroot)*yv(iy,jy,iroot)*
     1                 zv(iz,jz,iroot)*w9(iroot)
  432          continue
               elfder(ij,1,ipt) = elfder(ij,1,ipt) + dumx*dij 
            endif
  450       continue
  440    continue
c
  500 continue
c
 6000 continue
 7000 continue 
c
      return
 9997 format(' in -elfgrd- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots = ',3i3)
      end
