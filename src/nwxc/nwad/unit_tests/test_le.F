C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the .LE. operator
C>
C> This is an NWAD unit test. The value of the .LE. operator is compared
C> against the result of this operator on traditional double precision
C> values. The input data set is randomly generated.
C> 
      program test_le
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
        if ((nwad_left.le.nwad_right).neqv.(left.le.right)) then
          write(*,*)"nwad.le.nwad: fail",i,nwad_left.le.nwad_right,
     +                                     left.le.right
          write(*,*)
        endif
        if ((nwad_left.le.right).neqv.(left.le.right)) then
          write(*,*)"nwad.le.dble: fail",i,nwad_left.le.right,
     +                                     left.le.right
          write(*,*)
        endif
        if ((left.le.nwad_right).neqv.(left.le.right)) then
          write(*,*)"dble.le.nwad: fail",i,left.le.nwad_right,
     +                                     left.le.right
          write(*,*)
        endif
      enddo
      end
C> @}
c $Id$
