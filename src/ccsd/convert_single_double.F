      subroutine copy_32_to_64(n, a32, a64)
      implicit none
      !
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      !
      integer,       intent(in)  :: n
      real(kind=sp), intent(in)  :: a32(n)
      real(kind=dp), intent(out) :: a64(n)
      !
      integer :: i
#ifdef GNU_GE_4_8
      do concurrent (i=1:n)
#else
      do i=1,n
#endif
          a64(i) = real(a32(i), kind=dp)
      end do
      end subroutine

      subroutine copy_64_to_32(n, a64, a32)
      implicit none
      !
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      !
      integer,       intent(in)  :: n
      real(kind=dp), intent(in)  :: a64(n)
      real(kind=sp), intent(out) :: a32(n)
      !
      integer :: i
#ifdef GNU_GE_4_8
      do concurrent (i=1:n)
#else
      do i=1,n
#endif
          a32(i) = real(a64(i), kind=sp)
      end do
      end subroutine

      subroutine add_32_to_64(n, a32, a64)
      implicit none
      !
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      !
      integer,       intent(in)    :: n
      real(kind=sp), intent(in)    :: a32(n)
      real(kind=dp), intent(inout) :: a64(n)
      !
      integer :: i
#ifdef GNU_GE_4_8
      do concurrent (i=1:n)
#else
      do i=1,n
#endif
          a64(i) = a64(i) + real(a32(i), kind=dp)
      end do
      end subroutine

      subroutine add_64_to_32(n, a64, a32)
      implicit none
      !
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      !
      integer,       intent(in)    :: n
      real(kind=dp), intent(in)    :: a64(n)
      real(kind=sp), intent(inout) :: a32(n)
      !
      integer :: i
#ifdef GNU_GE_4_8
      do concurrent (i=1:n)
#else
      do i=1,n
#endif
          a32(i) = a32(i) + real(a64(i), kind=sp)
      end do
      end subroutine
