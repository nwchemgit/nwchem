C> \ingroup nwad_tests
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
      use nwad3
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
      double precision adx, ady, adz, adx2, ady2, adz2, adxy, adxz, adyz
      double precision adx3, ady3, adz3, adx2y, adx2z, adxy2, ady2z
      double precision adxz2, adyz2, adxyz
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
        adx   = inter_d3_dx(fx3)
        ady   = inter_d3_dx(fy3)
        adz   = inter_d3_dx(fz3)
        adx2  = inter_d3_dx2(fx3)
        ady2  = inter_d3_dx2(fy3)
        adz2  = inter_d3_dx2(fz3)
        adxy  = inter_d3_dxy(fx3,fx2y,fxy2,fy3)
        adxz  = inter_d3_dxy(fx3,fx2z,fxz2,fz3)
        adyz  = inter_d3_dxy(fy3,fy2z,fyz2,fz3)
        adx3  = inter_d3_dx3(fx3)
        ady3  = inter_d3_dx3(fy3)
        adz3  = inter_d3_dx3(fz3)
        adx2y = inter_d3_dx2y(fx3,fx2y,fxy2,fy3)
        adx2z = inter_d3_dx2y(fx3,fx2z,fxz2,fz3)
        adxy2 = inter_d3_dx2y(fy3,fxy2,fx2y,fx3)
        ady2z = inter_d3_dx2y(fy3,fy2z,fyz2,fz3)
        adxz2 = inter_d3_dx2y(fz3,fxz2,fx2z,fx3)
        adyz2 = inter_d3_dx2y(fz3,fyz2,fy2z,fy3)
        adxyz = inter_d3_dxyz(fx3,fy3,fz3,fx2y,fx2z,fxy2,fy2z,
     +                        fxz2,fyz2,fxyz)
        if (abs((dfdx-adx)/(dfdx+adx)).gt.tol) then
          write(*,*)"DX  : fail:",i,tmp(1),tmp(2),tmp(3),adx,dfdx
          write(*,*)
        endif
        if (abs((dfdy-ady)/(dfdy+ady)).gt.tol) then
          write(*,*)"DY  : fail:",i,tmp(1),tmp(2),tmp(3),ady,dfdy
          write(*,*)
        endif
        if (abs((dfdz-adz)/(dfdz+adz)).gt.tol) then
          write(*,*)"DZ  : fail:",i,tmp(1),tmp(2),tmp(3),adz,dfdz
          write(*,*)
        endif
        if (abs((d2fdx2-adx2)/(d2fdx2+adx2)).gt.tol) then
          write(*,*)"DX2 : fail:",i,tmp(1),tmp(2),tmp(3),adx2,d2fdx2
          write(*,*)
        endif
        if (abs((d2fdy2-ady2)/(d2fdy2+ady2)).gt.tol) then
          write(*,*)"DY2 : fail:",i,tmp(1),tmp(2),tmp(3),ady2,d2fdy2
          write(*,*)
        endif
        if (abs((d2fdz2-adz2)/(d2fdz2+adz2)).gt.tol) then
          write(*,*)"DZ2 : fail:",i,tmp(1),tmp(2),tmp(3),adz2,d2fdz2
          write(*,*)
        endif
        if (abs((d2fdxy-adxy)/(d2fdxy+adxy)).gt.tol) then
          write(*,*)"DXY : fail:",i,tmp(1),tmp(2),tmp(3),adxy,d2fdxy
          write(*,*)
        endif
        if (abs((d2fdxz-adxz)/(d2fdxz+adxz)).gt.tol) then
          write(*,*)"DXZ : fail:",i,tmp(1),tmp(2),tmp(3),adxz,d2fdxz
          write(*,*)
        endif
        if (abs((d2fdyz-adyz)/(d2fdyz+adyz)).gt.tol) then
          write(*,*)"DYZ : fail:",i,tmp(1),tmp(2),tmp(3),adyz,d2fdyz
          write(*,*)
        endif
        if (abs((d3fdx3-adx3)/(d3fdx3+adx3)).gt.tol) then
          write(*,*)"DX3 : fail:",i,tmp(1),tmp(2),tmp(3),adx3,d3fdx3
          write(*,*)
        endif
        if (abs((d3fdy3-ady3)/(d3fdy3+ady3)).gt.tol) then
          write(*,*)"DY3 : fail:",i,tmp(1),tmp(2),tmp(3),ady3,d3fdy3
          write(*,*)
        endif
        if (abs((d3fdz3-adz3)/(d3fdz3+adz3)).gt.tol) then
          write(*,*)"DZ3 : fail:",i,tmp(1),tmp(2),tmp(3),adz3,d3fdz3
          write(*,*)
        endif
        if (abs((d3fdx2y-adx2y)/(d3fdx2y+adx2y)).gt.tol) then
          write(*,*)"DX2Y: fail:",i,tmp(1),tmp(2),tmp(3),adx2y,d3fdx2y
          write(*,*)
        endif
        if (abs((d3fdx2z-adx2z)/(d3fdx2z+adx2z)).gt.tol) then
          write(*,*)"DX2Z: fail:",i,tmp(1),tmp(2),tmp(3),adx2z,d3fdx2z
          write(*,*)
        endif
        if (abs((d3fdxy2-adxy2)/(d3fdxy2+adxy2)).gt.tol) then
          write(*,*)"DXY2: fail:",i,tmp(1),tmp(2),tmp(3),adxy2,d3fdxy2
          write(*,*)
        endif
        if (abs((d3fdy2z-ady2z)/(d3fdy2z+ady2z)).gt.tol) then
          write(*,*)"DY2Z: fail:",i,tmp(1),tmp(2),tmp(3),ady2z,d3fdy2z
          write(*,*)
        endif
        if (abs((d3fdxz2-adxz2)/(d3fdxz2+adxz2)).gt.tol) then
          write(*,*)"DXZ2: fail:",i,tmp(1),tmp(2),tmp(3),adxz2,d3fdxz2
          write(*,*)
        endif
        if (abs((d3fdyz2-adyz2)/(d3fdyz2+adyz2)).gt.tol) then
          write(*,*)"DYZ2: fail:",i,tmp(1),tmp(2),tmp(3),adyz2,d3fdyz2
          write(*,*)
        endif
        if (abs((d3fdxyz-adxyz)/(d3fdxyz+adxyz)).gt.tol) then
          write(*,*)"DXYZ: fail:",i,tmp(1),tmp(2),tmp(3),adxyz,d3fdxyz
          write(*,*)
        endif
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_dxyz(x,y,z,f)
      use nwad3
      implicit none
      type(nwad_dble) :: x, y, z, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(y)
      f = exp(c * s + cos(x*z) + sin(y*z))
      end
C> @}
c $Id$
