C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the \f$d(x)\f$ extractor
C>
C> This is an NWAD unit test. The NWAD module calculates derivatives with
C> respect to the sum of all active variables. 
C> This program tests the computation of
C> \f$\frac{d^2f}{d(x)}\f$.
C> 
      program test_dxy
      use nwad
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: x, y, fx, fy
      double precision dfdx, dfdy
      double precision tmp(2)
      double precision tol
      parameter( tol = 1.0d-10)
      call random_seed
      do i = 1, npt
        call random_number(tmp)
        tmp = tmp*2.0d0*acos(-1.0d0)
        call submaxima_dx(tmp(1),tmp(2),dfdx,dfdy)
        x = active(tmp(1))
        y = inactive(tmp(2))
        call subad_dx(x,y,fx)
        x = inactive(tmp(1))
        y = active(tmp(2))
        call subad_dx(x,y,fy)
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
        write(*,*)i,fx%d1
        write(*,*)i,dfdx
        write(*,*)'dfdy:'
        write(*,*)i,fy%d1
        write(*,*)i,dfdy
        write(*,*)
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_dx(x,y,f)
      use nwad
      implicit none
      type(nwad_dble) :: x, y, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(y)
      f = exp(c * s)
      end
C> @}
