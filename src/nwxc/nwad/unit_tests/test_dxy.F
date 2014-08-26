C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the \f$d(xy)\f$ extractor
C>
C> This is an NWAD unit test. The NWAD module calculates derivatives with
C> respect to the sum of all active variables. As a result 2nd order derivatives
C> with 2 active variables generate a linear combination of partial derivatives.
C> To obtain individual partial derivatives these quantities need to be 
C> extracted from the sum. This program tests the computation of
C> \f$\frac{d^2f}{d(xy)}\f$ given \f$\frac{d^2f}{dx^2}\f$,
C> \f$\frac{d^2f}{dy^2}\f$ and \f$\frac{d^2f}{d(x+y)^2}\f$.
C> 
      program test_dxy
      use nwad3
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: x, y, fx2, fy2, fxy
      double precision dfdx, dfdy, d2fdx2, d2fdxy, d2fdy2
      double precision adx, ady, adx2, adxy, ady2
      double precision tmp(2)
      double precision tol
      parameter( tol = 1.0d-10)
      call random_seed
      do i = 1, npt
        call random_number(tmp)
        tmp = tmp*2.0d0*acos(-1.0d0)
        call submaxima_dxy(tmp(1),tmp(2),dfdx,dfdy,d2fdx2,d2fdxy,d2fdy2)
        x = active(tmp(1),2)
        y = inactive(tmp(2))
        call subad_dxy(x,y,fx2)
        x = inactive(tmp(1))
        y = active(tmp(2),2)
        call subad_dxy(x,y,fy2)
        x = active(tmp(1))
        y = active(tmp(2))
        call subad_dxy(x,y,fxy)
        adx  = inter_d2_dx(fx2)
        ady  = inter_d2_dx(fy2)
        adx2 = inter_d2_dx2(fx2)
        ady2 = inter_d2_dx2(fy2)
        adxy = inter_d2_dxy(fx2,fxy,fy2)
        if (abs((dfdx-adx)/(dfdx+adx)).gt.tol) then
          write(*,*)"DX  : fail:",i,tmp(1),tmp(2),adx,dfdx
          write(*,*)
        endif
        if (abs((dfdy-ady)/(dfdy+ady)).gt.tol) then
          write(*,*)"DY  : fail:",i,tmp(1),tmp(2),ady,dfdy
          write(*,*)
        endif
        if (abs((d2fdx2-adx2)/(d2fdx2+adx2)).gt.tol) then
          write(*,*)"DX2 : fail:",i,tmp(1),tmp(2),adx2,d2fdx2
          write(*,*)
        endif
        if (abs((d2fdxy-adxy)/(d2fdxy+adxy)).gt.tol) then
          write(*,*)"DXY : fail:",i,tmp(1),tmp(2),adxy,d2fdxy
          write(*,*)
        endif
        if (abs((d2fdy2-ady2)/(d2fdy2+ady2)).gt.tol) then
          write(*,*)"DY2 : fail:",i,tmp(1),tmp(2),ady2,d2fdy2
          write(*,*)
        endif
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_dxy(x,y,f)
      use nwad3
      implicit none
      type(nwad_dble) :: x, y, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(y)
      f = exp(c * s)
c     f = c * s
      end
C> @}
c $Id$
