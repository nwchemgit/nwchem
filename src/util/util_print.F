      logical function util_print(name, level)
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
#include "inp.fh"
      character*(*) name
      integer level
      integer ind
c     
      data depth /1/, print_level(1)/print_default/ ! In common
      data n_print(1) /0/, n_no_print(1) /0/
c     
      if (n_print(depth) .gt. 0) then
         if (inp_match(n_print(depth), .false., name, 
     $        print_names(1,depth), ind)) then
            util_print = .true.
*            write(6,*) ' up: explicit on for: ', name
            goto 10
         else if (ind .eq. 0) then
*            write(6,*) ' util_print: ambiguous print option: ', name
            call util_flush(6)
         endif
      endif
c     
      if (n_no_print(depth) .gt. 0) then
         if (inp_match(n_no_print(depth), .false., name, 
     $        no_ print_names(1,depth), ind)) then
            util_print = .false.
*            write(6,*) ' up: explicit off for: ', name
            goto 10
         else if (ind .eq. 0) then
*            write(6,*) ' util_print: ambiguous print option: ', name
            call util_flush(6)
         endif
      endif
c     
*      write(6,*) ' up: general rule for: ', name
      util_print = level .le. print_level(depth)

 10   continue
*      write(6,*) ' up: ', name, ' ', level, print_level(depth), 
*     $     util_print
c     
      end
      subroutine util_print_set_level(level)
      implicit none
#include "printcommon.fh"
      integer level
c
      print_level(depth) = level
c
      end
      subroutine util_print_push
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
c
      if (depth .eq. maxdepth) then
         call errquit('util_print_push: too deep',depth)
      else
         depth = depth + 1
         print_level(depth) = print_default
         n_print(depth) = 0
         n_no_print(depth) = 0
      endif
c
      end
      subroutine util_print_pop
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
c
      if (depth .eq. 1) then
         call errquit('util_print_pop: too high',depth)
      else
         depth = depth - 1
      endif
c
      end
      subroutine util_print_enable(name)
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
      character*(*) name
c
      if (n_print(depth) .eq. maxnames) then
         call errquit('util_print_enable: too many names', 
     $        n_print(depth))
      else
         n_print(depth) = n_print(depth) + 1
         print_names(n_print(depth), depth) = name
*         write(6,*) ' up: enabling: ', name
      endif
c
      end
      subroutine util_print_disable(name)
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
      character*(*) name
c
      if (n_no_print(depth) .eq. maxnames) then
         call errquit('util_print_disable: too many names', 
     $        n_no_print(depth))
      else
         n_no_print(depth) = n_no_print(depth) + 1
         no_print_names(n_no_print(depth), depth) = name
*         write(6,*) ' up: disabling: ', name
      endif
c
      end
      logical function util_print_rtdb_load(rtdb, printentry, 
     $     noprintentry)
      implicit none
#include "printcommon.fh"
#include "printlevels.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb
      character*(*) printentry, noprintentry
c     
      integer maxnames1
      parameter (maxnames1 = maxnames + 1)
      character*32 buf(maxnames1)
      integer i
c     
      call util_print_set_level(print_default)
      do i = 1, maxnames1
         buf(i) = ' '
      enddo
c     
      if (rtdb_cget(rtdb,printentry,maxnames1,buf)) then
         do i = 1, maxnames1
            if (inp_compare(.false.,'none',buf(i))) then
               call util_print_set_level(print_none)
               buf(i) = ' '
               goto 10
            else if (inp_compare(.false.,'low',buf(i))) then
               call util_print_set_level(print_low)
               buf(i) = ' '
               goto 10
            else if (inp_compare(.false.,'medium',buf(i))) then
               call util_print_set_level(print_medium)
               buf(i) = ' '
               goto 10
            else if (inp_compare(.false.,'high',buf(i))) then
               call util_print_set_level(print_high)
               buf(i) = ' '
               goto 10
            else if (inp_compare(.false.,'debug',buf(i))) then
               call util_print_set_level(print_debug)            
               buf(i) = ' '
               goto 10
            else if (inp_compare(.false.,'default',buf(i))) then
               call util_print_set_level(print_default)
               buf(i) = ' '
               goto 10
            endif
         enddo
 10      continue
c     
         do i = 1, maxnames1
            if (buf(i) .ne. ' ') call util_print_enable(buf(i))
         enddo
      endif
c     
      do i = 1, maxnames1
         buf(i) = ' '
      enddo
c     
      if (rtdb_cget(rtdb,noprintentry,maxnames,buf)) then
         do i = 1, maxnames1
            if (buf(i) .ne. ' ') call util_print_disable(buf(i))
         enddo
      endif
c     
      end
