c
c $Id: util_io.F,v 1.3 2008-04-14 22:04:21 marat Exp $
c

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
      character*(*) buffer

      return
      write(LuOut,*) buffer,ga_nodeid()
      call util_flush(LuOut) 
      end

