C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the \f$d(x)\f$ extractor
C>
C> This is an NWAD unit test. The NWAD module calculates derivatives with
C> respect to the sum of all active variables. 
C> This program tests the interpolation formula for
C> \f$\frac{df}{dx}\f$.
C> 
      program test_dx
      use nwad3
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: x, y, fx, fy
      double precision adx, ady
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
        adx = inter_d1_dx(fx)
        ady = inter_d1_dx(fy)
        if (abs((dfdx-adx)/(dfdx+adx)).gt.tol) then
          write(*,*)"DX  : fail:",i,tmp(1),tmp(2),adx,dfdx
          write(*,*)
        endif
        if (abs((dfdy-ady)/(dfdy+ady)).gt.tol) then
          write(*,*)"DY  : fail:",i,tmp(1),tmp(2),ady,dfdy
          write(*,*)
        endif
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_dx(x,y,f)
      use nwad3
      implicit none
      type(nwad_dble) :: x, y, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(y)
      f = exp(c * s)
      end
C> @}
c $Id$
