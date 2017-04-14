C$Id$

      block data bd_printcommon
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
      data depth /1/, print_level(1)/print_default/ ! In common
      data n_print(1) /0/, n_no_print(1) /0/
      data oprint_print(1) /.false./
      end

      logical function util_print(name, level)
      implicit none
#include "printlevels.fh"
#include "printcommon.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) name
      external bd_printcommon
      integer level
      integer ind
c     
      if (n_print(depth) .gt. 0) then
         if (inp_match(n_print(depth), .false., name, 
     $        print_names(1,depth), ind)) then
            util_print = .true.
            if (oprint_print(depth)) then
               write(LuOut,*) ' up: explicit on for: ', name
            endif
            goto 10
         else if (ind .eq. 0) then
            write(LuOut,*) ' util_print: ambiguous print option: ', name
            call util_flush(LuOut)
         endif
      endif
c     
      if (n_no_print(depth) .gt. 0) then
         if (inp_match(n_no_print(depth), .false., name, 
     $        no_print_names(1,depth), ind)) then
            util_print = .false.
            if (oprint_print(depth)) then
               write(LuOut,*) ' up: explicit off for: ', name
            endif
            goto 10
         else if (ind .eq. 0) then
            write(LuOut,*) ' util_print: ambiguous print option: ', name
            call util_flush(LuOut)
         endif
      endif
c     
      if (oprint_print(depth)) then
         write(LuOut,*) ' up: general rule for: ', name
      endif
      util_print = level .le. print_level(depth)

 10   continue
      if (oprint_print(depth)) then
         write(LuOut,*) ' up: ', name, ' ', depth,
     &        level, print_level(depth), 
     $        util_print
      endif
c     
      end
      subroutine util_print_get_level(level)
      implicit none
#include "printcommon.fh"
      integer level
c
      level = print_level(depth)
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
#include "errquit.fh"
#include "printlevels.fh"
#include "printcommon.fh"
c
      if (depth .eq. maxdepth) then
         call errquit('util_print_push: too deep',depth, INPUT_ERR)
      else
         depth = depth + 1
*  Next line for no inheritance
*         print_level(depth) = print_default
*  This line for inheritance
         print_level(depth) = print_level(depth-1)
*
         n_print(depth) = 0
         n_no_print(depth) = 0
         oprint_print(depth) = .false.
      endif
c
      end
      subroutine util_print_pop
      implicit none
#include "errquit.fh"
#include "printlevels.fh"
#include "printcommon.fh"
c
      if (depth .eq. 1) then
         call errquit('util_print_pop: too high',depth, INPUT_ERR)
      else
         depth = depth - 1
      endif
c
      end
      subroutine util_print_enable(name)
      implicit none
#include "errquit.fh"
#include "printlevels.fh"
#include "printcommon.fh"
#include "stdio.fh"
      character*(*) name
c
      if (n_print(depth) .eq. maxnames) then
         call errquit('util_print_enable: too many names', 
     $        n_print(depth), INPUT_ERR)
      else
         if (name(1:5) .eq. 'print') then
            oprint_print(depth) = .true.
         else
            n_print(depth) = n_print(depth) + 1
            print_names(n_print(depth), depth) = name
         endif
         if (oprint_print(depth)) then
            write(LuOut,*) ' up: enabling: ', name
         endif
      endif
c
      end
      subroutine util_print_disable(name)
      implicit none
#include "errquit.fh"
#include "printlevels.fh"
#include "printcommon.fh"
#include "stdio.fh"
      character*(*) name
c
      if (n_no_print(depth) .eq. maxnames) then
         call errquit('util_print_disable: too many names', 
     $        n_no_print(depth), INPUT_ERR)
      else
         n_no_print(depth) = n_no_print(depth) + 1
         no_print_names(n_no_print(depth), depth) = name
         if (oprint_print(depth)) then
            write(LuOut,*) ' up: disabling: ', name
         endif
      endif
c
      end
      subroutine util_print_rtdb_load(rtdb, prefix)
      implicit none
#include "printcommon.fh"
#include "printlevels.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb
      character*(*) prefix
c     
      character*80 printentry, noprintentry
      integer maxnames1
      parameter (maxnames1 = maxnames + 1)
      character*32 buf(maxnames1)
      integer i
c
      printentry = ' '
      noprintentry = ' '
      write(printentry,1) prefix(1:inp_strlen(prefix)),
     $     ':print'
      write(noprintentry,1) prefix(1:inp_strlen(prefix)),
     $     ':noprint'
 1    format(a,a)
c     
***      call util_print_set_level(print_default)
      do i = 1, maxnames1
         buf(i) = ' '
      enddo
c
      n_print(depth) = 0
      n_no_print(depth) = 0
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
      subroutine util_print_input(rtdb, prefix)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb
      character*(*) prefix
c
c     Parse input for print and noprint storing result in the database
c     as "<prefix>:print" and "<prefix>:noprint".
c
c     Input is left at end of line.
c
      integer maxnames
      parameter (maxnames = 33)
      character*32 buf(maxnames), test
      character*80 entry
      integer n
c
c     figure out if this is a print/noprint directive and build
c     the name for the database entry
c
      call inp_set_field(0)
      if (.not. inp_a(test))
     $     call errquit('util_print_input: can''t read directive',0,
     &       INPUT_ERR)
c
      entry = ' '
      if (inp_compare(.false.,'print',test)) then
         write(entry,'(a,a)') prefix(1:inp_strlen(prefix)), ':print'
      else if (inp_compare(.false.,'noprint',test)) then
         write(entry,'(a,a)') prefix(1:inp_strlen(prefix)), ':noprint'
      else
         call errquit('util_print_input: not print or noprint?',0,
     &       INPUT_ERR)
      endif
c
      n = 0
 10   if (inp_a(test)) then
         n = n + 1
         if (n .gt. maxnames) call errquit
     $        ('util_print_input: too many names?',n, INPUT_ERR)
         buf(n) = test
         goto 10
      endif
c
      if (.not. rtdb_cput(rtdb,entry,n,buf)) call errquit
     $     ('util_print_input: failed to store in database',rtdb,
     &       RTDB_ERR)
c
      end
      subroutine util_print_common_data
      implicit none
#include "printcommon.fh"
#include "printlevels.fh"
#include "stdio.fh"
      integer id,in
*
      write(luout,*)' ------------ print info -------------- '
      write(luout,*)' depth        = ',depth
      write(luout,*)' print_level  = ',print_level
      write(luout,*)' n_print      = ',n_print
      write(luout,*)' n_no_print   = ',n_no_print
      write(luout,*)' oprint_print = ',oprint_print
      write(luout,*)' '
      do id = 1,depth
        write(luout,*)' at depth ------------------------ ',id
        write(luout,*)'   names  ',
     &        (print_names(in,id),in=1,n_print(id))
        write(luout,*)'   nonames',
     &        (no_print_names(in,id),in=1,n_no_print(id))
        write(luout,*)' '
      enddo
      end
