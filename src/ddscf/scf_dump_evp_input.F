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
      integer :: style
      double precision :: thresh
      integer :: funit
      integer :: ga_type
      integer :: dimfx, dimfy
      integer :: dimsx, dimsy
      double precision, allocatable :: local(:)
      integer :: alloc_error
      logical :: toobig
      logical :: nodezero
!
      nodezero = ga_nodeid().eq.0
!
! style:
!   1 = most verbose format, use maximum memory, regular notation
!   2 = most verbose format, use maximum memory, scientific notation
!   3 = most verbose format, conserve memory, regular notation
!   4 = most verbose format, conserve memory, scientific notation
! +10 = add debugging via ga_print
!
      if (.not.rtdb_get(rtdb,'scf:evp_format', MT_INT, 1, style)) then
        style = 1
      endif
      ! make this configurable later
      thresh = 1.e-20
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
!
! these checks matter because I did not validate the code for non-square
! matrices, or even non-symmetric matrices, because F and S are Hermitian.
!
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
      if (style.ge.1) then
        if (nodezero) then
          allocate( local(1:dimfx*dimfy), stat=alloc_error)
          if (alloc_error.ne.0) then
            toobig = .true.
          else
            ! open file - make this customizable later
            open(unit=funit,file="eigensystem.dat",
     &           action="write",status="replace")
            ! gather and print F
            write(funit,*) 'Fock matrix'
            write(funit,*) 'dimfx=',dimfx
            write(funit,*) 'dimfy=',dimfy
            call ga_get(g_f, 1, dimfx, 1, dimfy, local, dimfx)
            call scf_dump_local(dimfx, dimfy, local,
     &                          style, funit, thresh)
            ! gather and print S
            write(funit,*) 'Overlap matrix'
            write(funit,*) 'dimsx=',dimsx
            write(funit,*) 'dimsy=',dimsy
            call ga_get(g_s, 1, dimsx, 1, dimsy, local, dimsx)
            call scf_dump_local(dimfx, dimfy, local,
     &                          style, funit, thresh)
            ! free memory
            deallocate( local, stat=alloc_error)
            if (alloc_error .ne. 0) then
                call errquit('deallocate failed', alloc_error, MA_ERR)
            endif
            ! close file
            close(funit)
          endif
        endif
      endif
      if ((style.ge.3) .or. toobig) then
        call errquit('memory-conserving implementation unfinished',
     &               0, UNKNOWN_ERR)
      endif
      ! debugging
      if (style.ge.11) then
        call ga_print(g_f)
        call ga_print(g_s)
      endif
      end subroutine

      subroutine scf_dump_local(dimx, dimy, local,
     &                          style, fileunit, thresh)
      implicit none
      integer, intent(in)  :: dimx, dimy
      double precision, intent(in) :: local(dimx,dimy)
      integer, intent(in) :: style, fileunit
      double precision, intent(in) :: thresh
      integer :: i, j
      if (style.eq.1) then
        do i = 1, dimx
          do j = 1, dimy
            if ( dabs(local(i,j)) .ge. thresh) then
              write(fileunit, 1) i,' , ',j,' , ',local(i,j)
            endif
          enddo
        enddo
      else if (style.eq.2) then
        do i = 1, dimx
          do j = 1, dimy
            if ( dabs(local(i,j)) .ge. thresh) then
              write(fileunit, 2) i,' , ',j,' , ',local(i,j)
            endif
          enddo
        enddo
      endif
    1 format(i6,a3,i6,a3,f30.15)
    2 format(i6,a3,i6,a3,e30.20e2)
      end subroutine scf_dump_local



