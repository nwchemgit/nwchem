*
* $Id$
*
      subroutine v_dirac(n2ft3d,ispin,dn,xcp,xce,x)
      implicit double precision(a-h, o-z)
      implicit integer (i-n)
#include "util.fh"

      integer n2ft3d
      integer ispin
      real*8 dn(n2ft3d,2)
      real*8 xcp(n2ft3d,2)
      real*8 xce(n2ft3d,2)
      real*8 x(n2ft3d)

      parameter (one3rd=1.0d0/3.0d0,for3rd=4.0d0/3.0d0)
      parameter (one6th=1.0d0/6.0d0)
      parameter (dncut=1.0d-30)

*---- parameters given by vosko et al -----------------*
      parameter (ap = 3.109070d-02, af = 1.554530d-02)
      parameter (x0p=-1.049800d-01, x0f=-3.250000d-01)
      parameter (bp = 3.727440d+00, bf = 7.060420d+00)
      parameter (cp = 1.293520d+01, cf = 1.805780d+01)
*------------------------------------------------------*

*     constants calculated from vosko's parameters
      parameter (xp  = -4.581653d-01,  xf  = -5.772521d-01)
      parameter (qp  =  6.151991d+00,  qf  =  4.730927d+00)
      parameter (xx0p=  1.255491d+01,  xx0f=  1.586879d+01)
      parameter (cp1 =  3.109070d-02,  cf1 =  1.554530d-02)
      parameter (cp2 =  9.690228d-04,  cf2 =  2.247860d-03)
      parameter (cp3 =  1.049800d-01,  cf3 =  3.250000d-01)
      parameter (cp4 =  3.878329d-02,  cf4 =  5.249122d-02)
      parameter (cp5 =  3.075995d+00,  cf5 =  2.365463d+00)
      parameter (cp6 =  1.863720d+00,  cf6 =  3.530210d+00)
      parameter (dp1 =  6.218140d-02,  df1 =  3.109060d-02)
      parameter (dp2 =  1.938045d-03,  df2 =  4.495720d-03)
      parameter (dp3 =  1.049800d-01,  df3 =  3.250000d-01)
      parameter (dp4 = -3.205972d-02,  df4 = -1.779316d-02)
      parameter (dp5 = -1.192972d-01,  df5 = -1.241661d-01)
      parameter (dp6 =  1.863720d+00,  df6 =  3.530210d+00)
      parameter (dp7 =  9.461748d+00,  df7 =  5.595417d+00)
      parameter (fc  =  1.923661d+00,  fd  =  2.564881d+00)
      parameter (crs =  7.876233d-01)

      call nwpw_timing_start(4)
      pi=4.0d0*datan(1.0d0)

*     square root of wigner radius
!$OMP DO
      do 100 k=1,n2ft3d
         rho=dn(k,1)+dn(k,ispin)+dncut
         x(k)=crs/rho**one6th
  100 end do
!$OMP END DO

*     paramagnetic exchange energy & potential
!$OMP DO
      do 120 k=1,n2ft3d
        xce(k,1)=(xp/x(k)**2)
        xcp(k,1)=for3rd*(xp/x(k)**2)
  120 end do
!$OMP END DO

*     return if spin-restricted dirac
      if(ispin.eq.1) go to 200

*     ferromagnetic exchange-energy & potential
!$OMP DO
      do 140 k=1,n2ft3d
        xce(k,2)=(xf/x(k)**2)
        xcp(k,2)=for3rd*(xf/x(k)**2)
  140 end do
!$OMP END DO

*     spin polarized exchange potential
!$OMP DO
      do 150 k=1,n2ft3d
        rho=dn(k,1)+dn(k,ispin)+dncut
        zup=2.0d0*dn(k,1)/rho
        zdw=2.0d0*dn(k,2)/rho
        f=(zup*(zup**one3rd)+zdw*(zdw**one3rd)-2.0d0)*fc
        xcp(k,1)=(1.0d0-f)*xcp(k,1)+f*xcp(k,2)
        df=(zup**one3rd-zdw**one3rd)*(xce(k,2)-xce(k,1))*fd
        xcp(k,2)=xcp(k,1)-zup*df
        xcp(k,1)=xcp(k,1)+zdw*df
        xce(k,1)=xce(k,1)+f*(xce(k,2)-xce(k,1))
  150 end do
!$OMP END DO

  200 continue
      
      call nwpw_timing_end(4)
      return
      end
