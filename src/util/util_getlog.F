! unfinished
      subroutine util_getlog(arg)
          use iso_c_binding, only: c_char, c_null_char, c_size_t
          implicit none
          character*(*), intent(inout) :: arg
          interface
            integer function getlogin_r(buf,bufsize) 
     &                       bind(C, name="getlogin_r")
              use iso_c_binding, only: c_char, c_size_t
              character(kind=c_char), intent(inout) :: buf(*)
              integer(kind=c_size_t), intent(in), value :: bufsize
            end function
          end interface
          integer :: rc
          integer(kind=c_size_t), parameter :: bufsize = 1024
          character(kind=c_char) :: buf(bufsize)
          buf = c_null_char
          rc = getlogin_r(buf,bufsize)
          print*,size(arg)
      end subroutine util_getlog
