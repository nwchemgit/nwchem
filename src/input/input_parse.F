      logical function input_parse(rtdb)
C     $Id: input_parse.F,v 1.63 2001-08-15 15:57:11 d3e129 Exp $
      implicit none
      integer rtdb              ! [input]
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
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
      integer nlines            ! No. of lines of input read in this call
c
      nlines = 0
      found_task = 0            ! 0=false, 1=true
c
      if (ga_nodeid() .eq. 0) then
c     
         if (util_print('banner',print_low)) then
           write(LuOut,*)
           write(LuOut,*)
           call util_print_centered(LuOut, 
     &          'NWChem Input Module',40,.true.)
           write(LuOut,*)
           write(LuOut,*)
         endif
c     
c     Do data base operations in sequential mode
c     
         status = rtdb_parallel(.false.)
c     
c     While (input is available)
c     
 10      if (inp_read()) then
            nlines = nlines + 1
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
            else if (inp_compare(.false.,test,'python')) then
               call python_input(rtdb)
            else if (inp_compare(.false.,test,'cosmo')) then
               call cosmo_input(rtdb)
            else if (inp_compare(.false.,test,'ecp')) then
               call bas_input(rtdb)
            else if (inp_compare(.false.,test,'so')) then
               call bas_input(rtdb)
            else if (inp_compare(.false.,test,'integral').or.
                     inp_compare(.false.,test,'intgrl')) then
               call intgrl_input(rtdb)
            else if (inp_compare(.false.,test,'scf')) then
               call scf_input(rtdb)
            else if (inp_compare(.false.,test,'mp2')) then
               call mp2_input(rtdb)
            else if (inp_compare(.false.,test,'time')) then
               call input_time(rtdb)
            else if (inp_compare(.false.,test,'set')) then
               call input_set(rtdb)
            else if (inp_compare(.false.,test,'drdy').or.
     &               inp_compare(.false.,test,'dirdyvtst').or.
     &               inp_compare(.false.,test,'drdy_nwchem')) then
               call drdy_input(rtdb)
            else if (inp_compare(.false.,test,'stop')) then
               call errquit('input_parse: STOP', 0)
            else if (inp_compare(.false.,test,'title')) then
               call input_title(rtdb)
            else if (inp_compare(.false.,test,'stepper')) then
               call stepper_input(rtdb)
            else if (inp_compare(.false.,test,'driver')) then
               call driver_input(rtdb)
            else if (inp_compare(.false.,test,'charge')) then
               call input_charge(rtdb)
            else if (inp_compare(.false.,test,'dft')) then
               call dft_input(rtdb)
            else if (inp_compare(.false.,test,'prepare')) then
               call pre_input(rtdb)
            else if (inp_compare(.false.,test,'md')) then
               call md_input(rtdb)
            else if (inp_compare(.false.,test,'esp')) then
               call esp_input(rtdb)
            else if (inp_compare(.false.,test,'analysis')) then
               call ana_input(rtdb)
            else if (inp_compare(.false.,test,'analyze')) then
               call ana_input(rtdb)
            else if (inp_compare(.false.,test,'ideaz')) then
               call ideaz_input(rtdb)
            else if (inp_compare(.false.,test,'gradients')) then
               call gradients_input(rtdb)
            else if (inp_compare(.false.,test,'gradient')) then
               call gradients_input(rtdb)
            else if (inp_compare(.false.,test,'ccsd')) then
               call ccsd_input(rtdb)
            else if (inp_compare(.false.,test,'oniom')) then
               call oniom_input(rtdb)
            else if (inp_compare(.false.,test,'print')) then
               call util_print_input(rtdb,' ')
               call util_print_rtdb_load(rtdb,' ') ! High level print
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
            else if (inp_compare(.false.,test,'nbo')) then
               call nbo_input(rtdb)
            else if (inp_compare(.false.,test,'speech')) then
               call speech_input(rtdb)
            else if (inp_compare(.false.,test,'gapss')) then
               call gap_parse(rtdb)
            else if (inp_compare(.false.,test,'pspw')) then
               call pspw_input(rtdb)
            else if (inp_compare(.false.,test,'pause')) then
	       call util_sleep(60)
            else if (inp_compare(.false.,test,'qmmm')) then
               call qmmm_input(rtdb)
            else if (inp_compare(.false.,test,'relativistic')) then
               call rel_input(rtdb)
            else if (inp_compare(.false.,test,'vscf')) then
               call vscf_input(rtdb)
            else if (inp_compare(.false.,test,'vib')) then
               call freq_vib_input(rtdb)
* this reordered test will match freq and freq*
            else if (inp_compare(.false.,'freq',test)) then  
               call freq_vib_input(rtdb)
            else if (inp_compare(.false.,test,'hessian')) then
               call hess_input(rtdb)
c
c     INSERT NEW ELSE IF CLAUSES ABOVE HERE
c
            else if (inp_compare(.false.,test,'memory')) then
               continue         ! Already digested by input_mem_size
            else if (inp_compare(.false.,test,'start')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'restart')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'continue')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'echo')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'scratch_dir')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'permanent_dir')) then
               continue         ! Already digested by input_file_info
            else if (inp_compare(.false.,test,'ecce_print')) then
               continue         ! Already digested by input_file_info
            else
               call errquit('input_parse: unknown directive',0)
            endif
c
            goto 10 
         endif                  ! End while
 9999    status = rtdb_parallel(.true.) ! Restore parallel access to db
         if (nlines.gt.0 .and. found_task.eq.0 .and. 
     $       util_print('warning',print_low)) then
c
c     Read input but found no task directive.  Not an error
c     but is probably undesirable ... unless doing python loop
c
            write(LuOut,*) 
            write(LuOut,*) ' ! warning: processed input with no task'
            write(LuOut,*)
            call util_flush(LuOut)
         endif
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

