! DFTD4 interface for NWChem.
! Developed with assistance from OpenAI Codex
! Implementation and numerical validation
! were reviewed and tested by MPH.
subroutine nw_dftd4_get_params(name, values, status) bind(C)
  use, intrinsic :: iso_c_binding, only : c_char, c_double, c_int, c_null_char
  use dftd4, only : damping_param, rational_damping_param, &
       get_rational_damping
  implicit none
  character(kind=c_char), intent(in) :: name(*)
  real(c_double), intent(out) :: values(6)
  integer(c_int), intent(out) :: status
  character(len=:), allocatable :: method
  class(damping_param), allocatable :: param
  integer :: i, n

  values = 0.0_c_double
  status = 1_c_int
  n = 0
  do i = 1, 256
     if (name(i) == c_null_char) exit
     n = n + 1
  end do
  if (n == 0 .or. n == 256) return

  allocate(character(len=n) :: method)
  do i = 1, n
     method(i:i) = achar(iachar(name(i)))
  end do

  call get_rational_damping(method, param, 1.0_c_double)
  if (.not.allocated(param)) return
  select type (param)
  type is (rational_damping_param)
     values = [param%s6, param%s8, param%s9, param%a1, param%a2, &
          param%alp]
     status = 0_c_int
  end select
end subroutine nw_dftd4_get_params
