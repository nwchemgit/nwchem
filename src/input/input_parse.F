      logical function input_parse(rtdb)
C     $Id: input_parse.F,v 1.25 1996-07-26 16:30:15 d3g681 Exp $
      implicit none
      integer rtdb              ! [input]
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
c     
c     All processes call this routine.  Only process 0 actually reads the 
c     input.  All other processes skip to the end of the routine
c     to wait in the broadcast of the return value.
c
c     Process 0 reads input until a task directive is encountered.
c     After this is read (by task_input) process 0 branches to the 
c     broadcast at the end and all processes return a value of TRUE.  
c     Otherwise FALSE is returned by all processes.
c
      logical status
      integer found_task
      character*255 test
c
      found_task = 0            ! 0=false, 1=true
c
      if (ga_nodeid() .eq. 0) then
c     
         write(6,*)
         write(6,*)
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
 10      if (inp_read()) then
            if (.not. inp_a(test))
     $           call errquit('input_parse: failed to read keyword', 0)
            call inp_prev_field
c     
            if (inp_compare(.false.,test,'task')) then
               call task_input(rtdb)
               found_task = 1
               goto 9999
            else if (inp_compare(.false.,test,'geometry')) then
               call geom_input(rtdb)
            else if (inp_compare(.false.,test,'basis')) then
               call bas_input(rtdb)
            else if (inp_compare(.false.,test,'scf')) then
               call scf_input(rtdb)
            else if (inp_compare(.false.,test,'time')) then
               call input_time(rtdb)
            else if (inp_compare(.false.,test,'set')) then
               call input_set(rtdb)
            else if (inp_compare(.false.,test,'stop')) then
               call errquit('input_parse: STOP', 0)
            else if (inp_compare(.false.,test,'title')) then
               call input_title(rtdb)
            else if (inp_compare(.false.,test,'stepper')) then
               call stepper_input(rtdb)
            else if (inp_compare(.false.,test,'charge')) then
               call input_charge(rtdb)
            else if (inp_compare(.false.,test,'dft')) then
               call dft_input(rtdb)
            else if (inp_compare(.false.,test,'md')) then
               call md_input(rtdb)
            else if (inp_compare(.false.,test,'gradient')) then
               call gradients_input(rtdb)
            else if (inp_compare(.false.,test,'ccsd')) then
               call ccsd_input(rtdb)
            else if (inp_compare(.false.,test,'print')) then
               call util_print_input(rtdb,' ')
            else if (inp_compare(.false.,test,'noprint')) then
               call util_print_input(rtdb,' ')
            else if (inp_compare(.false.,test,'set')) then
               call input_set(rtdb)
            else if (inp_compare(.false.,test,'unset')) then
               call input_unset(rtdb)
            else if (inp_compare(.false.,test,'mcscf')) then
               call mcscf_input(rtdb)
            else if (inp_compare(.false.,test,'plnwv')) then
               call plnwv_input(rtdb)
            else if (inp_compare(.false.,test,'dplot')) then
               call dplot_input(rtdb)
            else if (inp_compare(.false.,test,'property')) then
               call prop_input(rtdb)
            else if (inp_compare(.false.,test,'speech')) then
               call speech_input(rtdb)
c
c     INSERT NEW ELSE IF CLAUSES HERE
c
            else if (inp_compare(.false.,test,'memory')) then
               continue         ! Already digested by input_mem_size
            else if (inp_compare(.false.,test,'start')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'restart')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'echo')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'scratch_dir')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'permanent_dir')) then
               continue         ! Already digested by input_file_info
            else
               write(6,*) ' Skipping unknown directive '
               call inp_outrec()
            endif
c
            goto 10 
         endif                  ! End while
 9999    status = rtdb_parallel(.true.) ! Restore parallel access to db
      endif                     ! End node zero only
c     
c     Global sum of integer variable for simple broadcast
c
      call ga_sync()
      call ga_igop(33,found_task,1,'+')
      input_parse = (found_task .eq. 1)
c     
      end
      subroutine input_time(rtdb)
      implicit none
      integer rtdb
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

