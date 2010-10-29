c
c $Id$
c

#ifdef USE_SUBGROUPS
      block data util_stdio_data
      implicit none
#include "stdio.fh"
c     
      data  luout /6/

      end
#endif
      block data initial_util_io_data
      implicit none
#include "util_io_data.fh"
c     
      data io_on /.false./
      data io_unit /6/

      end

      function util_get_io_unit(fn)

      implicit none
      integer fn
      logical util_get_io_unit
c 
      integer k
      logical ostatus
c
      do k=80,90
        INQUIRE(UNIT=k,OPENED=ostatus)
        ostatus = .not.ostatus
        if(ostatus) 
     >    INQUIRE(UNIT=k,EXIST=ostatus)
        if(ostatus) then
          fn = k
          util_get_io_unit = .true.
          return
        end if 
      end do
      util_get_io_unit = .false.
      return
      end

      subroutine util_print_xyz(un,nt,c,t)
      implicit none
      integer un
      integer nt
      double precision c(3,nt)
      character*16 t(nt)

      integer i

      write(un,FMT=1) nt

      do i=1,nt
       write(un,FMT=2)
     >            t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00

      end do
1     FORMAT(1X,I9,/, "   ")
2     FORMAT(1X,A5,6X,3(F12.6,4X))

      end

      subroutine util_print_tagged_3d_array(un,nt,scale,c,t)
      implicit none
      integer un
      integer nt
      double precision c(3,nt)
      character*16 t(nt)
      double precision scale

      integer i

      do i=1,nt
       write(un,FMT=2)
     >            i,t(i),c(1,i)*scale,
     >            c(2,i)*scale,
     >            c(3,i)*scale

      end do
2     FORMAT(1X,I5,5X,A16,6X,3(F12.6,4X))

      end

      subroutine util_message(buffer)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "util_io_data.fh"
      character*(*) buffer

      if(.not.io_on) return
c
      write(io_unit,*) buffer,ga_nodeid()
      call util_flush(io_unit) 
      return
      end

      subroutine util_io_enable()
      implicit none
#include "stdio.fh"
#include "util_io_data.fh"
      io_on = .true.
      return
      end

      subroutine util_io_disable()
      implicit none
#include "stdio.fh"
#include "util_io_data.fh"
      io_on = .false.
      return
      end

