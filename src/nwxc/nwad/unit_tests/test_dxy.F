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
      use nwad
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: x, y, fx2, fy2, fpxpy, fpxmy, fmxpy, fmxmy
      double precision dfdx, dfdy, d2fdx2, d2fdxy, d2fdy2, adfxy
      double precision tmp(2)
      double precision tol
      parameter( tol = 1.0d-10)
      call random_seed
      do i = 1, npt
        call random_number(tmp)
        tmp = tmp*2.0d0*acos(-1.0d0)
        call submaxima_dxy(tmp(1),tmp(2),dfdx,dfdy,d2fdx2,d2fdxy,d2fdy2)
        x = active(tmp(1))
        x%d1 = 2
        y = inactive(tmp(2))
        call subad_dxy(x,y,fx2)
        x = inactive(tmp(1))
        y = active(tmp(2))
        y%d1 = 2
        call subad_dxy(x,y,fy2)
        x = active(tmp(1))
        y = active(tmp(2))
        call subad_dxy(x,y,fpxpy)
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
        write(*,*)i,fx2%d1*0.5d0
        write(*,*)i,dfdx
        write(*,*)'dfdy:'
        write(*,*)i,fy2%d1*0.5d0
        write(*,*)i,dfdy
        write(*,*)'d2fdx2:'
        write(*,*)i,fx2%d2*0.25d0
        write(*,*)i,d2fdx2
        write(*,*)'d2fdy2:'
        write(*,*)i,fy2%d2*0.25d0
        write(*,*)i,d2fdy2
        write(*,*)'d2fdxy:'
        write(*,*)i,0.5d0*(fpxpy%d2-0.25d0*(fx2%d2+fy2%d2))
        write(*,*)i,d2fdxy
        write(*,*)
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_dxy(x,y,f)
      use nwad
      implicit none
      type(nwad_dble) :: x, y, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(y)
      f = exp(c * s)
c     f = c * s
      end
C> @}
