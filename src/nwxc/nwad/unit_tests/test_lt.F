C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the .LT. operator
C>
C> This is an NWAD unit test. The value of the .LT. operator is compared
C> against the result of this operator on traditional double precision
C> values. The input data set is randomly generated.
C> 
      program test_lt
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
        if ((nwad_left.lt.nwad_right).neqv.(left.lt.right)) then
          write(*,*)"nwad.lt.nwad: fail",i,nwad_left.lt.nwad_right,
     +                                     left.lt.right
          write(*,*)
        endif
        if ((nwad_left.lt.right).neqv.(left.lt.right)) then
          write(*,*)"nwad.lt.dble: fail",i,nwad_left.lt.right,
     +                                     left.lt.right
          write(*,*)
        endif
        if ((left.lt.nwad_right).neqv.(left.lt.right)) then
          write(*,*)"dble.lt.nwad: fail",i,left.lt.nwad_right,
     +                                     left.lt.right
          write(*,*)
        endif
      enddo
      end
C> @}
c $Id$
