      subroutine input_parse(rtdb)
      implicit none
      integer rtdb              ! [input]
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
c
c     Read the input looking for high-level directives or names
c     of known modules.  Once all the input has been read attempt
c     to provide sensible defaults for anything that is missing.
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      logical status, echo
      parameter (num_dirs = 19)
      character*8 dirs(num_dirs)
      character*255 test
      data dirs /
     $     'memory',   'rtdb',     'geometry', 'basis',    'scf',
     $     'task',     'time',     'set'     , 'enter',    'stop',
     $     'title',    'start',    'restart' , 'stepper',  'echo',
     $     'charge',   'dft',      'md',       'gradient'/
      data echo /.false./
c
c     Only process 0 is involved here
c
      if (ga_nodeid() .ne. 0) return
c
      call util_print_centered(6, 'NWCHEM Input Module', 40, .true.)
      write(6,*)
      write(6,*)
c
c     Do data base operations in sequential mode
c
      status = rtdb_parallel(.false.)
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
     $        1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900) ind
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
c     Stop ... literally just stop ... well allow for parallel world
c
 1000    call errquit('input_parse: STOP', 0)
         goto 10                ! Not executed
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
c     Stepper
c
 1400    call stepper_input(rtdb)
         goto 10
c
c     Echo
c
 1500    echo = .true.
         goto 10
c
c     Charge
c
 1600    call input_charge(rtdb)
         goto 10
c    
c        Dft
c    
 1700    call dft_input(rtdb)
         goto 10
c    
c        MD
c    
 1800    call md_input(rtdb)
         goto 10
c    
c        Gradients
c    
 1900    call gradients_input(rtdb)
         goto 10
c
c
      endif
c
99999 continue
c
c echo input deck if requested
c
      if (echo) call input_echo(5,6)
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
      status = rtdb_parallel(.true.) ! Restore parallel access to db
c
      end
      subroutine input_time(rtdb)
      implicit none
      integer rtdb
      end
      subroutine input_defaults(rtdb)
      implicit none
      integer rtdb
      end
      subroutine input_verify
      implicit none
      end
      subroutine input_charge(rtdb)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      double precision charge
      character*255 test
c
      if (.not. inp_a(test)) call errquit
     &    ('input_charge: failed to read ascii string',911)
      if (test.ne.'charge') call errquit
     &    ('input_charge: failed to read "charge"',911)
c
      if (.not. inp_f(charge)) call errquit
     $     ('input_charge: failed to read real number', 0)
c
      if (.not. rtdb_put(rtdb, 'charge', MT_DBL, 1, charge))
     $     call errquit
     $     ('input_charge: failed to write to rtdb', 0)
c
      end

