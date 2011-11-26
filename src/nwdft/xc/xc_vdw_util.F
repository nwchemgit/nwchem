c
c      Interpolation of C6 to account for coordination number dependence
c
      double precision function c6cn(iat,jat,cni,cnj)
      implicit none
#include "xc_vdw.fh"
      integer iat,jat,i,j
      double precision cni,cnj,top,bottom,dist,c6_ref
      double precision cna,cnb
c
      c6cn=0.0d0
      top=0.0d0
      bottom=0.0d0
      do i=1,maxcn(iat)
         do j=1,maxcn(jat)
            c6_ref=c6ab(iat,jat,i,j,1)
            if (c6_ref.gt.0.0d0) then
               cna=c6ab(iat,jat,i,j,2)
               cnb=c6ab(iat,jat,i,j,3)
               dist=(cna-cni)**2+(cnb-cnj)**2
               top=top+dexp(k3*dist)*c6_ref
               bottom=bottom+dexp(k3*dist)
            endif
         enddo
      enddo
c
      if (bottom.gt.0.0d0) then
         c6cn=top/bottom
      else
         c6cn=0.0d0
      endif
      return
      end
c
c     Analytical gradient of coordination number dependent C6
c
      subroutine c6_grad(grad,iat,jat,kat,x,z,n,cnij,cnijk)
      implicit none
#include "xc_vdw.fh"
      integer n,iat,jat,kat,z(n),i,j
      double precision x(3,n)
      double precision grad(3),cnij(3,n),cnijk(3,n,n),cni,cnj
      double precision crd_nr,cnik(3),cnjk(3)
      external crd_nr
      double precision t1,t2,dt1x,dt1y,dt1z,dt2x,dt2y,dt2z
      double precision tmp1,tmp2,tmp3,tmp4,fac1,fac2
c
      cni=crd_nr(iat,n,x,z)
      cnj=crd_nr(jat,n,x,z)
      if (iat.eq.kat) then
         cnik(1)=cnij(1,iat)
         cnik(2)=cnij(2,iat)
         cnik(3)=cnij(3,iat)
         cnjk(1)=cnijk(1,iat,jat)
         cnjk(2)=cnijk(2,iat,jat)
         cnjk(3)=cnijk(3,iat,jat)
      else
         cnik(1)=cnijk(1,kat,iat)
         cnik(2)=cnijk(2,kat,iat)
         cnik(3)=cnijk(3,kat,iat)
         cnjk(1)=cnijk(1,kat,jat)
         cnjk(2)=cnijk(2,kat,jat)
         cnjk(3)=cnijk(3,kat,jat)
      endif
      t1=0.0d0
      t2=0.0d0
      dt1x=0.0d0
      dt1y=0.0d0
      dt1z=0.0d0
      dt2x=0.0d0
      dt2y=0.0d0
      dt2z=0.0d0
      do i=1,maxcn(z(iat))
        do j=1,maxcn(z(jat))
          tmp1=c6AB(z(iat),z(jat),i,j,3)-cnj
          tmp2=c6AB(z(iat),z(jat),i,j,2)-cni
          tmp3=dexp(k3*(tmp1*tmp1+tmp2*tmp2))
          t1=t1+c6AB(z(iat),z(jat),i,j,1)*tmp3
          t2=t2+tmp3
          fac1=tmp3*k3*2.0d0
          fac2=fac1*c6AB(z(iat),z(jat),i,j,1)
          tmp4=(tmp2*cnik(1)+tmp1*cnjk(1))
          dt1x=dt1x+fac2*tmp4
          dt2x=dt2x+fac1*tmp4
          tmp4=(tmp2*cnik(2)+tmp1*cnjk(2))
          dt1y=dt1y+fac2*tmp4
          dt2y=dt2y+fac1*tmp4
          tmp4=(tmp2*cnik(3)+tmp1*cnjk(3))
          dt1z=dt1z+fac2*tmp4
          dt2z=dt2z+fac1*tmp4
        enddo
      enddo
c
c     Numerical test to avoid division by zero (following Grimme)
c
      if (t2.gt.0) then
        grad(1)=(dt1x*t2-dt2x*t1)/(t2**2)
        grad(2)=(dt1y*t2-dt2y*t1)/(t2**2)
        grad(3)=(dt1z*t2-dt2z*t1)/(t2**2)
      else
        grad(1)=0.0d0
        grad(2)=0.0d0
        grad(3)=0.0d0
      endif
      end
c
      double precision function c6ij_sk(i,j,z)
      implicit none
#include "xc_vdw.fh"
      integer i,j
      integer z(*)
c     
      integer zi,zj
      double precision ni,nj
      external xc_neff
      double precision xc_neff,ci,cj
c
      c6ij_sk=0.0d0
      zi=z(i)
      zj=z(j)
      if (zi.le.0.or.zj.le.0) then 
         c6ij_sk=0.0d0
         return
      end if
c
      ci=c6(zi)
      cj=c6(zj)
c
c     write(6,*) zi,zj,ci,cj
c
      if(ci.le.0.d0.or.cj.le.0.d0) then
         c6ij_sk=0.0d0
         return
      endif
c
      if(ivdw.eq.1)then
c version 1 
c     S. Grimme J Comp Chem 25, 1463 (2004)
        c6ij_sk =2.d0*ci*cj/(ci+cj)
        else if(ivdw.eq.2) then
c version 2
c     S. Grimme J Comp Chem 27, 1787 (2006)
          c6ij_sk =sqrt(ci*cj)
c version 3
c     S. Grimme J Chem Phys 132, 154104 (2010)
        else if(ivdw.eq.0) then
          ni=xc_neff(zi)
          nj=xc_neff(zj)
          c6ij_sk=2d0*(ci*ci*cj*cj*ni*nj)**(1d0/3d0)/
     /       ((nj*nj*ci)**(1d0/3d0)+(ni*ni*cj)**(1d0/3d0))
      endif
      return
      end
c
c      Derivatives of coordination number function
c
      subroutine crd_nr_der(n,xyz,iz,cnij,cnijk)
      implicit none
#include "xc_vdw.fh"
      integer n,iat,i,iz(n)
      double precision cnij(3,n),cnijk(3,n,n),xyz(3,n)
      double precision dx,dy,dz,r,cov_rad,expf,fac2,fac3
c
      do i=1,n
         cnij(1,i)=0.0d0
         cnij(2,i)=0.0d0
         cnij(3,i)=0.0d0
         do iat=1,n
            if (i.ne.iat) then
               dx=xyz(1,iat)-xyz(1,i)
               dy=xyz(2,iat)-xyz(2,i)
               dz=xyz(3,iat)-xyz(3,i)
               r=dsqrt(dx*dx+dy*dy+dz*dz)
               cov_rad=cov_table(iz(iat))+cov_table(iz(i))
               expf=dexp(-k1*((cov_rad/r)-1.0d0))
               fac2=1.0d0/(expf+1.0d0)
               fac3=k1*cov_rad*expf*fac2*fac2/(r*r*r) 
               cnij(1,i)=cnij(1,i)-fac3*dx
               cnij(2,i)=cnij(2,i)-fac3*dy
               cnij(3,i)=cnij(3,i)-fac3*dz
               cnijk(1,iat,i)=fac3*dx
               cnijk(2,iat,i)=fac3*dy
               cnijk(3,iat,i)=fac3*dz
            endif
         enddo
      enddo
      return
      end
c
c      Coordination numbers based on inverse damping function 
c
      double precision function crd_nr(iat,nat,xyz,iz)
      implicit none
#include "xc_vdw.fh"
      integer iat,nat,iz(nat),i
      double precision xyz(3,nat),dx,dy,dz,r,rcov
c       
      crd_nr=0.0d0
      do i=1,nat
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r=dsqrt(dx*dx+dy*dy+dz*dz)
            rcov=cov_table(iz(i))+cov_table(iz(iat))
            crd_nr=crd_nr+1.d0/(1.d0+dexp(-k1*(rcov/r-1.0d0)))
         endif
      enddo
      return
      end
c
c      f1dmp(x1,x2,ff)=(alpha/x2)*ff*ff*exp(-alpha*((x1/x2)-1.0d0))
c
      double precision function f1dmp(x1,x2,ff)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2,ff
      f1dmp=(alpha/x2)*ff*ff*dexp(-alpha*((x1/x2)-1.0d0))
      return
      end
c
c      fdmp3(x1,x2)=1.0d0/(1.0d0+6.0d0*(x1/x2)**-a)
c
      double precision function fdmp3(x1,x2,a)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2,a
      fdmp3=1.0d0/(1.0d0+6.0d0*((x1/x2)**(-a)))
      return
      end
c
c      fdmp(x1,x2)=1.0d0/(1.0d0+exp(-alpha*((x1/x2)-1.0d0)))
c
      double precision function fdmp(x1,x2)
      implicit none
#include "xc_vdw.fh"
      double precision x1,x2
      fdmp=1.0d0/(1.0d0+dexp(-alpha*((x1/x2)-1.0d0)))
      return
      end
c
      subroutine get_scaling_fac(scales6,scales8,scalesr6,scalesr8)
      implicit none
#include "cdft.fh"
#include "xc_vdw.fh"
      double precision scales6,scales8,scalesr6,scalesr8
c
c DFT-D1 version
c
      if( ivdw.eq.1) then
        scales8  = 0.0d0
        scalesr6 = 0.0d0
        scalesr8 = 0.0d0
c BLYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
          scales6 = 1.40d0
c PBE
        else if(xfac(7).eq.1.0d0.and.cfac(12).eq.1.0d0) then
         scales6 = 0.70d0
c BP86
         else if(xfac(3).eq.1.0d0.and.xfac(2).eq.1.0d0.and.
     &              cfac(1).eq.1.0d0.and.cfac(4).eq.1.d0) then
         scales6 = 1.30d0
c default
        else 
         scales6 = 1.00d0
         write(6,*)'WARNING:'
         write(6,*)'you are using HF or an untested DF together'
         write(6,*)'with the dispersion correction!'
         endif
       else if ( ivdw.eq.2) then
c
c DFT-D2 version
c
        scales8  = 0.0d0
        scalesr6 = 0.0d0
        scalesr8 = 0.0d0
c BLYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
           scales6 = 1.20d0
c PBE
        else if(xfac(7).eq.1.0d0.and.cfac(12).eq.1.0d0) then
           scales6 = 0.75d0
c TPSS
        else if(xfac(21).eq.1.0d0.and.cfac(27).eq.1.0d0) then
           scales6 = 1.00d0
c B3LYP
        else if (xccomb(2)) then
           scales6 = 1.05d0
c B97-D
        else if (xccomb(37)) then
           scales6 = 1.25d0
c SSB-D
        else if (xccomb(46)) then
           scales6 = 0.847455d0
c BP86
        else if(xfac(3).eq.1.0d0.and.xfac(2).eq.1.0d0.and.
     &          cfac(1).eq.1.0d0.and.cfac(4).eq.1.d0) then
          scales6 = 1.05d0
c B2-PLYP
        else if(xccomb(2).and.xfac(1).eq.0.53d0.and.
     &          xfac(3).eq.0.47d0.and.cfac(2).eq.0.73d0) then
          scales6 = 0.55d0
c revPBE
c         scales6 = 1.25d0
c PBE0
        else if (xccomb(10)) then
          scales6 = 0.60d0
c PW6B95
c         scales6 = 0.50d0
c B2GP-PLYP
c         scales6 = 0.40d0
         else
c default
           scales6 =1.20d0
           write(6,*)'WARNING:'
           write(6,*)'you are using HF or an untested DF together'
           write(6,*)'with the dispersion correction!'
        endif
      else if (ivdw.eq.3) then
c
c DFT-D3 version
c
        scales6  = 1.0d0
        scalesr8 = 1.0d0
c BLYP
        if(xfac(3).eq.1.0d0.and.cfac(2).eq.1.0d0) then
           scalesr6 = 1.094d0
           scales8  = 1.682d0
c BP86 
        else if(xfac(3).eq.1.0d0.and.xfac(2).eq.1.0d0.and.
     &          cfac(1).eq.1.0d0.and.cfac(4).eq.1.d0) then
           scalesr6 = 1.139d0
           scales8  = 1.683d0
c B97-D
        else if (xccomb(37)) then
           scalesr6 = 0.892d0
           scales8  = 0.909d0
c revPBE
c          scalesr6 = 0.923d0
c          scales8  = 1.010d0
c PBE96
        else if(xccomb(51)) then
           scalesr6 = 1.217d0
           scales8  = 0.722d0
c PBESOL
c          scalesr6 = 1.320d0
c          scales8  = 0.581d0
c RPW86-PBE
c          scalesr6 = 1.224d0
c          scales8  = 0.901d0
c rPBE
c          scalesr6 = 0.872d0
c          scales8  = 0.514d0
c TPSS
        else if(xfac(21).eq.1.0d0.and.cfac(27).eq.1.0d0) then
           scalesr6 = 1.166d0
           scales8  = 1.105d0
c B3LYP
        else if (xccomb(2)) then
           scalesr6 = 1.261d0
           scales8  = 1.703d0
        else if (xccomb(10)) then
c PBE0
           scalesr6 = 1.278d0
           scales8  = 0.928d0
c revPBE38
c          scalesr6 = 1.021d0
c          scales8  = 0.862d0
c PW6B95
        else if (xccomb(25)) then
           scalesr6 = 1.532d0
           scales8  = 1.682d0
c TPSS0
c          scalesr6 = 1.252d0
c          scales8  = 1.242d0
c B2-PLYP
        else if(xccomb(2).and.xfac(1).eq.0.53d0.and.
     &          xfac(3).eq.0.47d0.and.cfac(2).eq.0.73d0) then
           scales6  = 0.640d0
           scalesr6 = 1.094d0
           scales8  = 1.682d0
c PWPB95
c          scales6  = 0.820d0
c          scalesr6 = 1.557d0
c          scales8  = 0.705d0
c B2GP-PLYP
c          scales6  = 0.560d0
c          scalesr6 = 1.586d0
c          scales8  = 0.760d0
c PTPSS
c          scales6  = 0.750d0
c          scalesr6 = 1.541d0
c          scales8  = 0.879d0
c HF   
c          scalesr6 = 1.158d0
c          scales8  = 1.746d0
c mpwlyp
c          scalesr6 = 1.239d0
c          scales8  = 1.098d0
c bpbe
c          scalesr6 = 1.087d0
c          scales8  = 2.033d0
c bhlyp
        else if (xccomb(48)) then
           scalesr6 = 1.370d0
           scales8  = 1.442d0
c tpssh
        else if (xccomb(18)) then
           scalesr6 = 1.223d0
           scales8  = 1.219d0
c pwb6k
        else if (xccomb(26)) then
           scalesr6 = 1.660d0
           scales8  = 0.550d0
c b1b95
        else if (xccomb(23)) then
           scalesr6 = 1.613d0
           scales8  = 1.868d0
c bop
c          scalesr6 = 0.929d0
c          scales8  = 1.975d0
c o-lyp
c          scalesr6 = 0.806d0
c          scales8  = 1.764d0
c o-pbe
c          scalesr6 = 0.837d0
c          scales8  = 1.868d0
c ssb
        else if (xccomb(46)) then
           scalesr6 = 1.215d0
           scales8  = 0.663d0
c revssb
c          scalesr6 = 1.221d0
c          scales8  = 0.560d0
c otpss
c          scalesr6 = 1.128d0
c          scales8  = 1.494d0
c b3pw91
        else if (xccomb(50)) then
           scalesr6 = 1.176d0
           scales8  = 1.775d0
c revpbe0
c          scalesr6 = 0.949d0
c          scales8  = 0.792d0
c pbe38
c          scalesr6 = 1.333d0
c          scales8  = 0.998d0
c mpw1b95
        else if (xccomb(21)) then
           scalesr6 = 1.605d0
           scales8  = 1.118d0
c mpwb1k
        else if (xccomb(22)) then
           scalesr6 = 1.671d0
           scales8  = 1.061d0
c bmk
c          scalesr6 = 1.931d0
c          scales8  = 2.168d0
c cam-b3lyp
c          scalesr6 = 1.378d0
c          scales8  = 1.217d0
c lc-wpbe
c          scalesr6 = 1.355d0
c          scales8  = 1.279d0
c m05
        else if (xccomb(29)) then
           scalesr6 = 1.373d0
           scales8  = 0.595d0
c m052x
        else if (xccomb(30)) then
           scalesr6 = 1.417d0
           scales8  = 0.000d0
c m06l
        else if (xccomb(33)) then
           scalesr6 = 1.581d0
           scales8  = 0.000d0
c m06
        else if (xccomb(35)) then
           scalesr6 = 1.325d0
           scales8  = 0.000d0
c m062x
        else if (xccomb(36)) then
           scalesr6 = 1.619d0
           scales8  = 0.000d0
c m06hf
        else if (xccomb(34)) then
           scalesr6 = 1.446d0
           scales8  = 0.000d0
c dftb
c          scalesr6 = 1.699d0
c          scales8  = 1.504d0
        endif
      else
        scales6 = 0.0d0
        scales8 = 0.0d0
        scalesr6 = 0.0d0
        scalesr8 = 0.0d0
      endif
      return
      end
c
c     Small function to extract coordination number representation number
c     and atomic number from column two and three of the table
c
      subroutine limit(iat,jat,iadr,jadr)
      implicit none
      integer iat,jat,iadr,jadr
c
      iadr=1
      jadr=1
810   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 810
      endif
820   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 820
      endif
      return
      end
c
      double precision function xc_neff(z)
      implicit none
      integer z
c
      integer nv
c
      if(z.eq.1) then
         xc_neff=0.8d0
         return
      endif
      if(z.le.2) nv=z
      if(z.le.10) nv=z-2
      if(z.gt.10) stop 'neff stop '
      xc_neff=1.17d0+0.33d0*nv
      return
      end
