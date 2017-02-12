      subroutine scf_dump_evp_input(rtdb,g_f,g_s)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
!#include "geom.fh"
!#include "util.fh"
!
! call ga_diag( g_fock, g_overlap, g_mos, eigenvalues)
!
      integer, intent(in) :: rtdb
      integer, intent(in) :: g_f ! Fock matrix
      integer, intent(in) :: g_s ! overlap matrix
!
      integer :: style, funit
      integer :: ga_type
      integer :: dimfx, dimfy
      integer :: dimsx, dimsy
      double precision, allocatable :: local(:)
      integer :: alloc_error
      logical :: toobig
!
! style:
!   0 = most verbose format, use maximum memory
!   1 = most verbose format, conserve memory
!
      if (.not.rtdb_get(rtdb,'scf:evp_format', MT_INT, 1, style)) then
          style = 0
      endif
!
! do some basic consistency checks
!
      call ga_inquire(g_f, ga_type, dimfx, dimfy)
      if (ga_type .ne. MT_DBL) then
          call errquit('array type is not double', ga_type, GA_ERR)
      endif
      call ga_inquire(g_s, ga_type, dimsx, dimsy)
      if (ga_type .ne. MT_DBL) then
          call errquit('array type is not double', ga_type, GA_ERR)
      endif
      if ((dimfx .ne. dimsx).or.(dimfy .ne. dimsy)) then
          call errquit('F and S matrices not same size', 0, GA_ERR)
      endif
      if ((dimfx .ne. dimfy).or.(dimsx .ne. dimsy)) then
          call errquit('F and S matrices not square', 0, GA_ERR)
      endif
!
! gather and print
!
      toobig = .false.
      if (style==0) then
          allocate( local(1:dimfx*dimfy), stat=alloc_error)
          if (alloc_error.ne.0) then
              toobig = .true.
          else
              ! open file - make this customizable later
              open(unit=funit,file="eigensystem.dat",
     &             action="write",status="replace")
              ! gather and print F
              call ga_get(g_f, 1, dimfx, 1, dimfy, local, dimfx)
              call scf_dump_local(dimfx, dimfy, local, style, funit)
              ! gather and print S
              call ga_get(g_s, 1, dimsx, 1, dimsy, local, dimsx)
              call scf_dump_local(dimfx, dimfy, local, style, funit)
              ! free memory
              deallocate( local, stat=alloc_error)
              if (alloc_error .ne. 0) then
                  call errquit('deallocate failed', alloc_error, MA_ERR)
              endif
              ! close file
              close(funit)
          endif
      endif
      if (style==1 .or. toobig) then
          call errquit('memory-conserving implementation unfinished',
     &                 0, UNKNOWN_ERR)
      endif
      end subroutine

      subroutine scf_dump_local(dimx, dimy, local, style, fileunit)
      implicit none
      integer, intent(in)  :: dimx, dimy
      double precision, intent(in) :: local(dimx*dimy)
      integer, intent(in) :: style, fileunit
      integer :: i, j
      do i = 1, dimx
          do j = 1, dimy
              write(fileunit, 1) i,',',j,',',local(i*dimy+j)
          enddo
      enddo
    1 format(i6,a1,i6,a1,e30.20e5)
      end subroutine scf_dump_local



