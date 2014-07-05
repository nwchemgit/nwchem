C> \ingroup nwad
C> @{
C>
C> \brief test the implementation of the \f$d(xy)\f$ extractor
C>
C> This is an NWAD unit test. The NWAD module calculates derivatives with
C> respect to the sum of all active variables. As a result 3rd order derivatives
C> with 3 active variables generate a linear combination of partial derivatives.
C> To obtain individual partial derivatives these quantities need to be 
C> extracted from the sum. This program tests the computation of
C> \f$\frac{d^3f}{d(xyz)}\f$.
C> 
      program test_dxyz
      use nwad
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: x, y, z, fx3, fy3, fz3, fx2y, fx2z, fxy2, fxyz
      type(nwad_dble) :: fy2z, fxz2, fyz2
      double precision dfdx,   dfdy,   dfdz
      double precision d2fdx2, d2fdy2, d2fdz2, d2fdxy, d2fdxz, d2fdyz
      double precision d3fdx3, d3fdy3, d3fdz3
      double precision d3fdx2y,d3fdx2z,d3fdxy2,d3fdy2z,d3fdxz2,d3fdyz2
      double precision d3fdxyz
      double precision tmp(3)
      double precision tol
      parameter( tol = 1.0d-10)
      double precision f13, f14, f16, f18, f19, f118, f127, f136, f181
      double precision f23, f29, f227, f527, f536
      double precision f572, f5162
      parameter( f13   = 1.0d0/3.0d0)
      parameter( f14   = 1.0d0/4.0d0)
      parameter( f16   = 1.0d0/6.0d0)
      parameter( f18   = 1.0d0/8.0d0)
      parameter( f19   = 1.0d0/9.0d0)
      parameter( f118  = 1.0d0/18.0d0)
      parameter( f127  = 1.0d0/27.0d0)
      parameter( f136  = 1.0d0/36.0d0)
      parameter( f181  = 1.0d0/81.0d0)
      parameter( f23   = 2.0d0/3.0d0)
      parameter( f29   = 2.0d0/9.0d0)
      parameter( f227  = 2.0d0/27.0d0)
      parameter( f527  = 5.0d0/27.0d0)
      parameter( f536  = 5.0d0/36.0d0)
      parameter( f572  = 5.0d0/72.0d0)
      parameter( f5162 = 5.0d0/162.0d0)
      call random_seed
      do i = 1, npt
        call random_number(tmp)
        tmp = tmp*2.0d0*acos(-1.0d0)
        call submaxima_dxyz(tmp(1),tmp(2),tmp(3),dfdx,dfdy,dfdz,
     +       d2fdx2,d2fdy2,d2fdz2,d2fdxy,d2fdxz,d2fdyz,
     +       d3fdx3,d3fdy3,d3fdz3,
     +       d3fdx2y,d3fdx2z,d3fdxy2,d3fdy2z,d3fdxz2,d3fdyz2,d3fdxyz)
        x = active(tmp(1),3)
        y = inactive(tmp(2))
        z = inactive(tmp(3))
        call subad_dxyz(x,y,z,fx3)
        x = inactive(tmp(1))
        y = active(tmp(2),3)
        z = inactive(tmp(3))
        call subad_dxyz(x,y,z,fy3)
        x = inactive(tmp(1))
        y = inactive(tmp(2))
        z = active(tmp(3),3)
        call subad_dxyz(x,y,z,fz3)
        x = active(tmp(1),2)
        y = active(tmp(2))
        z = inactive(tmp(3))
        call subad_dxyz(x,y,z,fx2y)
        x = active(tmp(1),2)
        y = inactive(tmp(2))
        z = active(tmp(3))
        call subad_dxyz(x,y,z,fx2z)
        x = active(tmp(1))
        y = active(tmp(2),2)
        z = inactive(tmp(3))
        call subad_dxyz(x,y,z,fxy2)
        x = inactive(tmp(1))
        y = active(tmp(2),2)
        z = active(tmp(3))
        call subad_dxyz(x,y,z,fy2z)
        x = active(tmp(1))
        y = inactive(tmp(2))
        z = active(tmp(3),2)
        call subad_dxyz(x,y,z,fxz2)
        x = inactive(tmp(1))
        y = active(tmp(2))
        z = active(tmp(3),2)
        call subad_dxyz(x,y,z,fyz2)
        x = active(tmp(1))
        y = active(tmp(2))
        z = active(tmp(3))
        call subad_dxyz(x,y,z,fxyz)
c       x = active(tmp(1))
c       y = active_neg(tmp(2))
c       call subad_dxy(x,y,fpxmy)
c       x = active_neg(tmp(1))
c       y = active(tmp(2))
c       call subad_dxy(x,y,fmxpy)
c       x = active_neg(tmp(1))
c       y = active_neg(tmp(2))
c       call subad_dxy(x,y,fmxmy)
c       adfxy = extract_dxdy(fx2,fxpy,fy2)
C       if (abs((d2fdxy-adfxy)/(d2fdxy+adfxy)).gt.tol) then
C         write(*,*)"D2FDXY: fail:",i,x%d0,y%d0,d2fdxy
C         write(*,*)"D2FDXY: fail:",i,x%d0,y%d0,adfxy
C         write(*,*)
C       endif
        write(*,*)'dfdx:'
        write(*,*)i,fx3%d1*f13
        write(*,*)i,dfdx
        write(*,*)'dfdy:'
        write(*,*)i,fy3%d1*f13
        write(*,*)i,dfdy
        write(*,*)'dfdz:'
        write(*,*)i,fz3%d1*f13
        write(*,*)i,dfdz
        write(*,*)'d2fdx2:'
        write(*,*)i,fx3%d2*f19
        write(*,*)i,d2fdx2
        write(*,*)'d2fdy2:'
        write(*,*)i,fy3%d2*f19
        write(*,*)i,d2fdy2
        write(*,*)'d2fdz2:'
        write(*,*)i,fz3%d2*f19
        write(*,*)i,d2fdz2
        write(*,*)'d2fdxy:'
        write(*,*)i,f18*(fx2y%d2+fxy2%d2)-f572*(fx3%d2+fy3%d2)
        write(*,*)i,d2fdxy
        write(*,*)'d2fdxz:'
        write(*,*)i,f18*(fx2z%d2+fxz2%d2)-f572*(fx3%d2+fz3%d2)
        write(*,*)i,d2fdxz
        write(*,*)'d2fdyz:'
        write(*,*)i,f18*(fy2z%d2+fyz2%d2)-f572*(fy3%d2+fz3%d2)
        write(*,*)i,d2fdyz
        write(*,*)'d3fdx3:'
        write(*,*)i,fx3%d3*f127
        write(*,*)i,d3fdx3
        write(*,*)'d3fdy3:'
        write(*,*)i,fy3%d3*f127
        write(*,*)i,d3fdy3
        write(*,*)'d3fdz3:'
        write(*,*)i,fz3%d3*f127
        write(*,*)i,d3fdz3
        write(*,*)'d3fdx2y:'
        write(*,*)i,-f5162*fx3%d3+f19*fx2y%d3-f118*fxy2%d3+f181*fy3%d3
        write(*,*)i,d3fdx2y
        write(*,*)'d3fdx2z:'
        write(*,*)i,-f5162*fx3%d3+f19*fx2z%d3-f118*fxz2%d3+f181*fz3%d3
        write(*,*)i,d3fdx2z
        write(*,*)'d3fdxy2:'
        write(*,*)i,-f5162*fy3%d3+f19*fxy2%d3-f118*fx2y%d3+f181*fx3%d3
        write(*,*)i,d3fdxy2
        write(*,*)'d3fdy2z:'
        write(*,*)i,-f5162*fy3%d3+f19*fy2z%d3-f118*fyz2%d3+f181*fz3%d3
        write(*,*)i,d3fdy2z
        write(*,*)'d3fdxz2:'
        write(*,*)i,-f5162*fz3%d3+f19*fxz2%d3-f118*fx2z%d3+f181*fx3%d3
        write(*,*)i,d3fdxz2
        write(*,*)'d3fdyz2:'
        write(*,*)i,-f5162*fz3%d3+f19*fyz2%d3-f118*fy2z%d3+f181*fy3%d3
        write(*,*)i,d3fdyz2
        write(*,*)'d3fdxyz:'
        write(*,*)i,f16*fxyz%d3+f181*(fx3%d3+fy3%d3+fz3%d3)
     +              -f136*(fx2y%d3+fx2z%d3+fxy2%d3+fy2z%d3
     +                     +fxz2%d3+fyz2%d3)
        write(*,*)i,d3fdxyz
        write(*,*)
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_dxyz(x,y,z,f)
      use nwad
      implicit none
      type(nwad_dble) :: x, y, z, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(y)
      f = exp(c * s + cos(x*z) + sin(y*z))
      end
C> @}
