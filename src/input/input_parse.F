      subroutine input_parse(rtdb)
      implicit none
      integer rtdb              ! [input]
#include "inp.fh"
#include "tcgmsg.fh"
c
c     Read the input looking for high-level directives or names
c     of known modules.  Once all the input has been read attempt
c     to provide sensible defaults for anything that is missing.
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 13)
      character*8 dirs(num_dirs)
      character*255 test
      data dirs /
     $     'memory',   'rtdb',     'geometry', 'basis',    'scf',
     $     'task',     'time',     'set'     , 'enter',    'stop',
     $     'title',    'start',    'restart'/
c
c     Only process 0 is involved here
c
      if (nodeid() .ne. 0) return
c
c     While (input is available)
c
      rewind 5
      call inp_init(5, 6)
 10   if (inp_read()) then
         if (.not. inp_a(test))
     $        call errquit('input_parse: failed to read keyword', 0)
         if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
            write(6,*) ' Skipping unknown directive '
            call inp_outrec()
            goto 10
         endif
c
         call inp_prev_field
         goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     $        1100, 1200, 1300) ind
         call errquit('input_parse: unimplemented directive', ind)
c
c     Memory ... already handled by input_mem_size() but will soon
c                want to do this in a context specific fashion
c
 100     goto 10
c
c     RTDB ... already handled by input_rtdb_name()
c
 200     goto 10
c
c     Geometry
c
 300     call geom_input(rtdb)
         goto 10
c
c     Basis
c
 400     call bas_input(rtdb)
         goto 10
c
c     SCF
c
 500     call scf_input(rtdb)
         goto 10
c
c     Task
c
 600     call task_input(rtdb)
         goto 10
c
c     Time
c
 700     call input_time(rtdb)
         goto 10
c
c     Set
c
 800     call input_set(rtdb)
         goto 10
c
c     Enter ... behave as if we just hit end of input
c
 900     goto 99999
c
c     Stop ... literally just stop
c
 1000    stop
c
c     Title
c
 1100    call input_title(rtdb)
         goto 10
c
c     Start ... already digested by input_start_opt
c
 1200    goto 10
c
c     Restart ... already digested by input_start_opt
c
 1300    goto 10
c
      endif
99999 continue
c
c     All input is exhausted.  Now apply defaults for everything that
c     is missing from the input if is a startup calculation
c
      call input_defaults(rtdb)
c
c     Try and verify anything not checked so far
c
      call input_verify
c
c     Done
c
      end
      subroutine scf_input
      end
      subroutine input_time
      end
      subroutine input_defaults
      end
      subroutine input_verify
      end
