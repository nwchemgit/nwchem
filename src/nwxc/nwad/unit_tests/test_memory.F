C> \ingroup nwad_tests
C> @{
C>
C> \brief test whether we can determine the size of an NWAD_DBLE
C>
C> In NWChem the MA library provides the normal memory allocation mechanisms
C> by partitioning a large memory array. The types it supports are 
C> LOGICALS, INTEGERS, REAL, DOUBLE PRECISION, COMPLEX, and DOUBLE COMPLEX.
C> In order to allocate arrays of type(nwad_dble) we need to play some tricks.
C> 
C> The basic trick is to work out how many double precision words it takes to
C> store one type(nwad_dble) entity, next allocate a memory region that is 
C> large enough to hold the required number of these things, then call a routine
C> passing the memory section and apply a type change through the declarations
C> in the routine called. 
C>
      program test_memory
      use nwad3
      implicit none
      integer         :: npt
      parameter (npt = 100)
      type(nwad_dble) :: x(2)
      integer         :: bytes, doubles, length
      double precision :: array(npt)
C
      bytes = loc(x(2))-loc(x(1))
      doubles = (bytes+8-1)/8
      length = npt/doubles
      write(*,*)'bytes   = ',bytes
      write(*,*)'doubles = ',doubles
      call init_nwad_dble(array,length)
      end
c
      subroutine init_nwad_dble(x,nx)
      use nwad3
      implicit none
      integer                        :: nx
      type(nwad_dble), intent(inout) :: x(nx)
      integer                        :: i
      do i = 1, nx
        x(i) = inactive(0.0d0)
      enddo
      end
C>
C> @}
      
c $Id$
