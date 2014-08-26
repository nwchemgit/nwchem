C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the .GE. operator
C>
C> This is an NWAD unit test. The value of the .GE. operator is compared
C> against the result of this operator on traditional double precision
C> values. The input data set is randomly generated.
C> 
      program test_ge
      use nwad3
      implicit none
      integer :: npt, i
      parameter (npt = 100)
      type(nwad_dble) :: nwad_left, nwad_right
      double precision left, right
      double precision tmp(2)
      call random_seed
      do i = 1, npt
        call random_number(tmp)
        left  = tmp(1)
        right = tmp(2)
        nwad_left  = left
        nwad_right = right
        if ((nwad_left.ge.nwad_right).neqv.(left.ge.right)) then
          write(*,*)"nwad.ge.nwad: fail",i,nwad_left.ge.nwad_right,
     +                                     left.ge.right
          write(*,*)
        endif
        if ((nwad_left.ge.right).neqv.(left.ge.right)) then
          write(*,*)"nwad.ge.dble: fail",i,nwad_left.ge.right,
     +                                     left.ge.right
          write(*,*)
        endif
        if ((left.ge.nwad_right).neqv.(left.ge.right)) then
          write(*,*)"dble.ge.nwad: fail",i,left.ge.nwad_right,
     +                                     left.ge.right
          write(*,*)
        endif
      enddo
      end
C> @}
c $Id$
