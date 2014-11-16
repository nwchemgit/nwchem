C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of printing the SIN function
C>
C> This is an NWAD unit test for the multivariate implementation.
C> The derivatives of SIN function are compared
C> against analytic derivatives. The function is 3-dimensional as that is
C> sufficient for this test. The input data set is randomly generated.
C> 
      program test_sin_p3
      use nwadp1
      implicit none
      integer :: npt, i, mmvar
      integer :: i1, i2, i3
      integer :: ix, ix2, ix3
      parameter (mmvar = 3)
      parameter (npt = 1)
      integer :: iv(mmvar)
      integer :: numvar
      type(nwad_dble) :: x, y, z, f
      double precision fa,dfa(mmvar),dfa2(mmvar*(mmvar+1)/2)
      double precision dfa3(mmvar*(mmvar+1)*(mmvar+2)/6), tol
      double precision tmp(3), tf, df
      parameter( tol = 1.0d-10)
      call random_seed
      open(UNIT=60,FILE="test_sin_p3.txt")
      call set_unitno(60)
      call set_name(f,"en")
      do i = 1, npt
        call random_number(tmp)
        tmp = tmp*2.0d0*acos(-1.0d0)
        x = set_gamma_aa(tmp(1))
        y = set_gamma_ab(tmp(2))
        z = set_gamma_bb(tmp(3))
        call submaxima_sin_p3(x%d0,y%d0,z%d0,fa,dfa)
        call subad_sin_p3(x,y,z,f)
        if (abs((fa-f%d0)/(fa+f%d0)).gt.tol) then
          write(*,*)"F  : fail:",i,x%d0,y%d0,z%d0,fa
          write(*,*)"F  : fail:",i,x%d0,y%d0,z%d0,f%d0
          write(*,*)
        endif
        numvar = get_nvar(f)
        do i1 = 1, numvar
          call get_d1(f,i1,tf,iv(i1))
          if (abs((dfa(i1)-tf)/(dfa(i1)+tf)).gt.tol) then
            write(*,10)i,i1,x%d0,y%d0,z%d0,dfa(i)
            write(*,10)i,i1,x%d0,y%d0,z%d0,tf
            write(*,*)
          endif
        enddo
      enddo
 10   format("DF : fail:",i4, i2,3e10.3,e18.10)
      end
C>
C> \brief The test routine
C>
      subroutine subad_sin_p3(x,y,z,f)
      use nwadp1
      implicit none
      type(nwad_dble) :: x, y, z, t, f
      t = sin(x) + sin(y) + sin(z)
      f = sin(t)
      end
C> @}
c $Id: test_add.F 26056 2014-08-26 19:03:07Z d3y133 $
