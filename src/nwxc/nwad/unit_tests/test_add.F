C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the ADD operator
C>
C> This is an NWAD unit test. The derivatives of ADD operator are compared
C> against analytic derivatives. The function is 1-dimensional as that is
C> sufficient for this test. The input data set is randomly generated.
C> 
      program test_add
      use nwad3
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: x(npt), f(npt)
      double precision fa(npt),dfa(npt),dfa2(npt),dfa3(npt), tol
      double precision tmp
      parameter( tol = 1.0d-10)
      call random_seed
      do i = 1, npt
        call random_number(tmp)
        tmp = tmp*2.0d0*acos(-1.0d0)
        x(i) = active(tmp)
        call submaxima_add(x(i)%d0,fa(i),dfa(i),dfa2(i),dfa3(i))
        call subad_add(x(i),f(i))
        if (abs((fa(i)-f(i)%d0)/(fa(i)+f(i)%d0)).gt.tol) then
          write(*,*)"F  : fail:",i,x(i)%d0,fa(i)
          write(*,*)"F  : fail:",i,x(i)%d0,f(i)%d0
          write(*,*)
        endif
        if (abs((dfa(i)-f(i)%d1)/(dfa(i)+f(i)%d1)).gt.tol) then
          write(*,*)"DF : fail:",i,x(i)%d0,dfa(i)
          write(*,*)"DF : fail:",i,x(i)%d0,f(i)%d1
          write(*,*)
        endif
        if (abs((dfa2(i)-f(i)%d2)/(dfa2(i)+f(i)%d2)).gt.tol) then
          write(*,*)"DF2: fail:",i,x(i)%d0,dfa2(i)
          write(*,*)"DF2: fail:",i,x(i)%d0,f(i)%d2
          write(*,*)
        endif
        if (abs((dfa3(i)-f(i)%d3)/(dfa3(i)+f(i)%d3)).gt.tol) then
          write(*,*)"DF3: fail:",i,x(i)%d0,dfa3(i)
          write(*,*)"DF3: fail:",i,x(i)%d0,f(i)%d3
          write(*,*)
        endif
      enddo
      end
C>
C> \brief The test routine
C>
      subroutine subad_add(x,f)
      use nwad3
      implicit none
      type(nwad_dble) :: x, f
      type(nwad_dble) :: c, s
      s = sin(x)
      c = cos(x)
      f = c + s
      end
C> @}
c $Id$
