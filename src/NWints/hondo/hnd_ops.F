      subroutine hnd_ops(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     & expj,coefj,j_nprim,j_ngen,Lj,nint,e11,scr,lscr,catms,nat)
c
c $Id$
c
c ... jochen: calculates integrals of the type
c     <i | u/r**3 | d/dv j> with u, v, = x,y,z and u, r being a
c     cartesian coordinate and distance relative to a nucleus.
c     for u.ne.v an antisymmetrized combination of these matrix
c     elements gives the paramagnetic nuclear spin electron orbit
c     (pso) hyperfine integrals as calculated in hnd_pso.F.
c     'ops' doesn't really mean anything but i picked the name
c     as a variant of 'pso'. 
c     
c     ----- Wrapper routine that sets the sizes of scratch blocks -----
c
      implicit double precision (a-h,o-z)
#include "hnd_pointers.fh"
      dimension scr(lscr)
      logical para,dia
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension e11(nint,3,3)
      dimension catms(3,nat)
c
c Use scratch for temps in integral calculation
c Scratch needs are
c 3*(Li+1)*(Lj+2)*((Li+Lj+1)/2+1)
c
c The pointers are defined in hnd_pointers.fh
c
      call hnd_ops1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,expj,
     1 coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,
     2 scr(dsot(1)),scr(dsot(2)),scr(dsot(3)),scr(dsot(4)))
c
      return
      end
c
      subroutine hnd_ops1(xyzi,expi,coefi,i_nprim,i_ngen,Li,xyzj,
     1 expj,coefj,j_nprim,j_ngen,Lj,catms,nat,nint,e11,
     4 v1,dv1,v1d,dv1d)
c
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "stdio.fh"
#include "errquit.fh"
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
      dimension Nxyz(3),xyzi(3),xyzj(3),expi(i_nprim),expj(j_nprim)
      dimension coefi(i_nprim,i_ngen),coefj(j_nprim,j_ngen)
      dimension e11(nint,3,3)
      dimension catms(3,nat)
      dimension    v1(3,Li+1,Lj+2,*)   ! blocks for vx1,vy1,vz1
      dimension   dv1(3,Li+1,Lj+2,*)   ! blocks for dvx1,dvy1,dvz1
      dimension   v1d(3,Li+1,Lj+2,*)   ! blocks for vx1d,vy1d,vz1d
      dimension  dv1d(3,Li+1,Lj+2,*)   ! blocks for dvx1d,dvy1d,dvz1d
      logical para,dia
      data rln10  /2.30258509299405d+00/   ! ln(10)
      data pi212  /1.12837916709551d+00/   ! 2/sqrt(pi)
      dimension tmp(3,3)
      integer dbg
c
      dbg = 0
c
      tol=rln10*itol
c
c     this routine does one atomic center at a time.
c     error exit if nat is not equal to one
      if (nat.ne.1) call errquit('hnd_ops: nat.ne.1',nat,UNKNOWN_ERR)
c
c     Zero integral array
c
      call dcopy(nint*9,0.0d0,0,e11,1)
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
         call errquit('hnd_ops: need higher Rys rt',nroots,INT_ERR)
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

c      do kat=1,nat
c ... jochen: here, we pick one center instead of 
c     going through the whole list. For now it is set to the
c     first center in the list.

      kat = 1

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
c           ----- for x/r**3 -----
c           
            v1(1,i,j,ir)=xint    ! vx1
            v1(2,i,j,ir)=yint    ! vy1
            v1(3,i,j,ir)=zint*ww ! vz1
            call hnd_dervxyz(1)       
            dv1(1,i,j,ir)=xint  ! dvx1
            dv1(2,i,j,ir)=yint  ! dvy1
            dv1(3,i,j,ir)=zint*ww ! dvz1
          enddo                 ! i-loop
        enddo                   ! j-loop
      enddo                     ! ir-loop
c     
c ... jochen 03/15:
c     looking at Dupuis' Comp.Phys.Commun. 134, 2001, 150,
c     v1 is used to calculate the basis function derivative,
c     while dv1 has the potential derivative.
c     This is all one needs for Eq. (2.61). For (2.60), we
c     also need the a basis function derivative term of the
c     potential derivative, which is presumably calculated 
c     like v1d, but from dv1 instead of v1.
c     The number of roots in both Eqs. are the same.
c
      do ir=1,nroots
c       
c       ----- d/dx ... operators -----
c       
        do i=1,lit
          v1d(1,i,1,ir)= (-(aj+aj)*v1(1,i,2,ir)) ! vx1d and vx1
          v1d(2,i,1,ir)= (-(aj+aj)*v1(2,i,2,ir)) ! vy1d and vy1
          v1d(3,i,1,ir)= (-(aj+aj)*v1(3,i,2,ir)) ! vz1d and vz1
          do j=2,ljt
            v1d(1,i,j,ir)= (dble(j-1)*v1(1,i,j-1,ir)- ! vx1d and vx1
     &        (aj+aj)*v1(1,i,j+1,ir))
            v1d(2,i,j,ir)= (dble(j-1)*v1(2,i,j-1,ir)- ! vy1d and vy1
     &        (aj+aj)*v1(2,i,j+1,ir))
            v1d(3,i,j,ir)= (dble(j-1)*v1(3,i,j-1,ir)- ! vz1d and vz1
     &        (aj+aj)*v1(3,i,j+1,ir))
          enddo                 ! j-loop
        enddo                   ! i-loop
        
c ... jochen: new integrals, compared to pso.
c             same as above but with v1->dv1
        do i=1,lit
          dv1d(1,i,1,ir)= (-(aj+aj)*dv1(1,i,2,ir))
          dv1d(2,i,1,ir)= (-(aj+aj)*dv1(2,i,2,ir)) 
          dv1d(3,i,1,ir)= (-(aj+aj)*dv1(3,i,2,ir))
          do j=2,ljt
            dv1d(1,i,j,ir)= (dble(j-1)*dv1(1,i,j-1,ir)- 
     &        (aj+aj)*dv1(1,i,j+1,ir))
            dv1d(2,i,j,ir)= (dble(j-1)*dv1(2,i,j-1,ir)- 
     &        (aj+aj)*dv1(2,i,j+1,ir))
            dv1d(3,i,j,ir)= (dble(j-1)*dv1(3,i,j-1,ir)- 
     &        (aj+aj)*dv1(3,i,j+1,ir))
          enddo                 ! j-loop
        enddo                   ! i-loop
c         
      enddo                     ! ir-loop
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

          tmp(:,:) = 0.0d0 ! will an ancient compiler complain here?

          do ir=1,nroots
          
c           operator = electric-field * derivative (D)
c           dv1 is for the field component, v1d for the
c           derivative, as long as the coordinates don't co-incide.
c           otherwise, we need dv1d

c           diagonal terms:

c           field-x D-x
            tmp(1,1) = tmp(1,1)
     &        + dv1d(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1(3,iz,jz,ir)
c           field-y D-y
            tmp(2,2) = tmp(2,2)
     &        + v1(1,ix,jx,ir)*dv1d(2,iy,jy,ir)*v1(3,iz,jz,ir)
c           field-z D-z
            tmp(3,3) = tmp(3,3)
     &        + v1(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1d(3,iz,jz,ir)

c           off-diagonal terms:

c           field-x D-y
            tmp(1,2) = tmp(1,2)
     &        + dv1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*v1(3,iz,jz,ir)
c           field-y D-x
            tmp(2,1) = tmp(2,1)
     &        + v1d(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1(3,iz,jz,ir)
c           field-x D-z
            tmp(1,3) = tmp(1,3)
     &        + dv1(1,ix,jx,ir)*v1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
c           field-z D-x
            tmp(3,1) = tmp(3,1)
     &        + v1d(1,ix,jx,ir)*v1(2,iy,jy,ir)*dv1(3,iz,jz,ir)
c           field-y D-z
            tmp(2,3) = tmp(2,3)
     &        + v1(1,ix,jx,ir)*dv1(2,iy,jy,ir)*v1d(3,iz,jz,ir)
c           field-z D-y
            tmp(3,2) = tmp(3,2)
     &        + v1(1,ix,jx,ir)*v1d(2,iy,jy,ir)*dv1(3,iz,jz,ir)
           

          enddo     
          ij=ij+1

          do k = 1,3
            do l = 1,3
              e11(ij,k,l)=e11(ij,k,l)+tmp(k,l)*cij
            end do
          end do
c         
        enddo                   ! j-loop final summation
      enddo                     ! i-loop final summation
c     
c     enddo        ! kat-loop
c
 1000 continue
      enddo        ! jprim loop
      enddo        ! iprim loop
c
      return
 9999 format(' in -ops- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,' lit,ljt,nroots= ',3i3)
      end
