C> \ingroup nwad_tests
C> @{
C>
C> \brief test the implementation of the .GT. operator
C>
C> This is an NWAD unit test. The value of the .GT. operator is compared
C> against the result of this operator on traditional double precision
C> values. The input data set is randomly generated.
C> 
      program test_gt
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
        if ((nwad_left.gt.nwad_right).neqv.(left.gt.right)) then
          write(*,*)"nwad.gt.nwad: fail",i,nwad_left.gt.nwad_right,
     +                                     left.gt.right
          write(*,*)
        endif
        if ((nwad_left.gt.right).neqv.(left.gt.right)) then
          write(*,*)"nwad.gt.dble: fail",i,nwad_left.gt.right,
     +                                     left.gt.right
          write(*,*)
        endif
        if ((left.gt.nwad_right).neqv.(left.gt.right)) then
          write(*,*)"dble.gt.nwad: fail",i,left.gt.nwad_right,
     +                                     left.gt.right
          write(*,*)
        endif
      enddo
      end
C> @}
c $Id$
