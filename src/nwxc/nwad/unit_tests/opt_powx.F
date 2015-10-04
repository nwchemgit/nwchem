c
c Program to test the powx function.
c ==================================
c
c Powx can benefit from a few optimizations that reduce the number
c of floating point operations. This code generates pseudo random
c input values for NWAD_DBLE X and DOUBLE PRECISION Y and computes
c result Z. All three quantities are printed for every data point.
c
cdefine NWAD_DEG 0
#include "../nwad.F"
cundef  NWAD_DEG
cdefine NWAD_DEG 1
cinclude "../nwad.F"
cdefine NWAD_PRINT y
cinclude "../nwad.F"
cundef  NWAD_DEG
cundef  NWAD_PRINT
cdefine NWAD_DEG 2
cinclude "../nwad.F"
cundef  NWAD_DEG
cdefine NWAD_DEG 3
cinclude "../nwad.F"
cundef  NWAD_DEG
c
c Dummy routine
c
      subroutine errquit(chr,i,j)
      implicit none
      character*(*) chr
      integer i
      integer j
      write(*,*)chr
      write(*,*)i,j
      stop
      end subroutine errquit
c
c Create NWAD_DBLE value for X
c
      subroutine create_nwad_dble(x)
c
#include "./nwad.fh"
c
      implicit none
c
      type(nwad_dble), intent(out) :: x
      integer                      :: ii
      real                         :: randnum(2)
      real :: r18, r28, r38, r48, r58, r68, r78
      double precision :: pi2
      parameter(r18 = 1.0d0/8.0d0)
      parameter(r28 = 2.0d0/8.0d0)
      parameter(r38 = 3.0d0/8.0d0)
      parameter(r48 = 4.0d0/8.0d0)
      parameter(r58 = 5.0d0/8.0d0)
      parameter(r68 = 6.0d0/8.0d0)
      parameter(r78 = 7.0d0/8.0d0)
      pi2 = 2.0d0*acos(-1.0d0)
c
      x = 0.0d0
      do ii = 1, 8
        call random_number(randnum)
        if (randnum(1).lt.r18) then
          x = sin(x+set_rho_a(pi2*randnum(2)))
        else if (randnum(1).lt.r28) then
          x = sin(x+set_rho_b(pi2*randnum(2)))
        else if (randnum(1).lt.r38) then
          x = sin(x+set_gamma_aa(pi2*randnum(2)))
        else if (randnum(1).lt.r48) then
          x = sin(x+set_gamma_ab(pi2*randnum(2)))
        else if (randnum(1).lt.r58) then
          x = sin(x+set_gamma_bb(pi2*randnum(2)))
        else if (randnum(1).lt.r68) then
          x = sin(x+set_tau_a(pi2*randnum(2)))
        else if (randnum(1).lt.r78) then
          x = sin(x+set_tau_b(pi2*randnum(2)))
        else
        endif
      enddo
c
      end subroutine create_nwad_dble
c
c Create DOUBLE PRECISION value for X
c
      subroutine create_dble(x)
c
      implicit none
c
      double precision, intent(out) :: x
      real                          :: randnum
c
      x = 0.0d0
      call random_number(randnum)
      x = randnum*10.0d0-5.0d0
c
      end subroutine create_dble
c
c Skip a number of random numbers
c
      subroutine skip(nn)
      implicit none
      integer, intent(in) :: nn
      real                :: randnum
      integer             :: ii
      do ii = 1, nn
        call random_number(randnum)
      enddo
      end
c
c Program to test stuff
c
      program opt_powx
#include "./nwad.fh"
      implicit none
      type(nwad_dble)  :: x, z
      double precision :: y
      integer          :: nn_begin
      integer          :: nn_skip
      integer          :: ii
      integer          :: clcnt ! command line argument count
      character*32     :: text
      clcnt = command_argument_count()
      if (clcnt.ne.2) then
        write(*,*)"expecting exactly 2 command line arguments"
        stop 10
      endif
      call get_command_argument(1,text)
      read(text,*)nn_begin
      call get_command_argument(2,text)
      read(text,*)nn_skip
      write(*,*)'skip at the start: ',nn_begin
      write(*,*)'skip in between  : ',nn_skip
      call skip(nn_begin)
      do ii = 1, 10
        call create_nwad_dble(x)
        x = abs(x)
        call create_dble(y)
        z = x**y
        write(*,*)
        write(*,*)"Inputs: X"
        call print(6,x)
        write(*,*)"Inputs: Y"
        write(*,'(1x,e14.6)')y
        write(*,*)"Output: Z"
        call print(6,z)
      enddo
c
      end

