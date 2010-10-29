*
*     $Id$
*

*     *********************************
*     *                               *
*     *      nwpw_scratch_write       *
*     *                               *
*     *********************************
      subroutine nwpw_scratch_write(tag,n,A)
      implicit none
      character*(*) tag
      integer n
      real*8 A(*)

*     **** local variables ****
      integer       i
      character*255 full_filename

      call util_file_name_noprefix(tag,
     >                             .true.,
     >                             .true.,
     >                             full_filename)

      open(unit=89,file=full_filename,form='unformatted')
      write(89) (A(i),i=1,n)
      close(89)
      return
      end


*     *********************************
*     *                               *
*     *      nwpw_scratch_read        *
*     *                               *
*     *********************************
      subroutine nwpw_scratch_read(tag,n,A)
      implicit none
      character*(*) tag
      integer n
      real*8 A(*)

*     **** local variables ****
      integer       i
      character*255 full_filename

      call util_file_name_noprefix(tag,
     >                             .true.,
     >                             .true.,
     >                             full_filename)

      open(unit=89,file=full_filename,status='old',form='unformatted')
      read(89) (A(i),i=1,n)
      close(89)
      return
      end

