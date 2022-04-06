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
      logical function input_parse(rtdb)

C     $Id$
c
      implicit none
c
#include "errquit.fh"
c
      integer rtdb              ! [input]
c
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
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
     $           call errquit('input_parse: failed to read keyword', 0,
     &       INPUT_ERR)
            call inp_prev_field
c     
            if (inp_compare(.false.,test,'task')) then
               call task_input(rtdb)
               found_task = 1
               goto 9999
            else if (inp_compare(.false.,test,'geometry')) then
               call geom_input(rtdb)
            else if (inp_compare(.false.,test,'bsse')) then
               call bsse_input(rtdb)
            else if (inp_compare(.false.,test,'basis')) then
               call bas_input(rtdb)
            else if (inp_compare(.false.,test,'python')) then
               call python_input(rtdb)
            else if (inp_compare(.false.,test,'cosmo')) then
               call cosmo_input(rtdb)
            else if (inp_compare(.false.,test,'dimqm')) then ! DIM/QM JEM
               call dimqm_input(rtdb)
            else if (inp_compare(.false.,test,'dimpar')) then ! DIM/QM JEM
               call dimpar_input(rtdb)
            else if (inp_compare(.false.,test,'ecp')) then
               call bas_input(rtdb)
            else if (inp_compare(.false.,test,'so')) then
               call bas_input(rtdb)
            else if (inp_compare(.false.,test,'integral').or.
     &               inp_compare(.false.,test,'intgrl')) then
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
               call errquit('input_parse: STOP', 0, INPUT_ERR)
            else if (inp_compare(.false.,test,'title')) then
               call input_title(rtdb)
            else if (inp_compare(.false.,test,'stepper')) then
               call stepper_input(rtdb)
            else if (inp_compare(.false.,test,'mepgs')) then
               call mepgs_input(rtdb)
            else if (inp_compare(.false.,test,'tropt')) then
               call tropt_input(rtdb)
            else if (inp_compare(.false.,test,'driver')) then
               call driver_input(rtdb)
            else if (inp_compare(.false.,test,'string')) then
               call string_input(rtdb)
            else if (inp_compare(.false.,test,'zts')) then
               call string_input(rtdb)
            else if (inp_compare(.false.,test,'qcharge')) then
               call input_qcharge(rtdb)
            else if (inp_compare(.false.,test,'charge')) then
               call input_charge(rtdb)
            else if (inp_compare(.false.,test,'dft')) then
               call dft_input(rtdb)
            else if (inp_compare(.false.,test,'qmd')) then
               call qmd_input(rtdb)
            else if (inp_compare(.false.,test,'occup')) then
               call occup_input(rtdb)
            else if (inp_compare(.false.,test,'prepare')) then
               call pre_input(rtdb)
            else if (inp_compare(.false.,test,'md')) then
               call md_input(rtdb)
            else if (inp_compare(.false.,test,'argos')) then
               call argos_input(rtdb)
            else if (inp_compare(.false.,test,'esp')) then
               call esp_input(rtdb)
            else if (inp_compare(.false.,test,'et')) then
               call et_input(rtdb)
            else if (inp_compare(.false.,test,'analysis')) then
               call ana_input(rtdb)
            else if (inp_compare(.false.,test,'analyze')) then
               call ana_input(rtdb)
            else if (inp_compare(.false.,test,'diana')) then
               call dia_input(rtdb)
            else if (inp_compare(.false.,test,'grad')) then
               call gradients_input(rtdb)
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
c           else if (inp_compare(.false.,test,'plnwv')) then
c              call plnwv_input(rtdb)
            else if (inp_compare(.false.,test,'mm')) then
               call mm_input(rtdb)
            else if (inp_compare(.false.,test,'bq')) then
               call bq_input(rtdb)
            else if (inp_compare(.false.,test,'constraints')) then
               call cons_input(rtdb)
            else if (inp_compare(.false.,test,'dplot')) then
               call dplot_input(rtdb)
            else if (inp_compare(.false.,test,'property')) then
               call prop_input(rtdb)
            else if (inp_compare(.false.,test,'speech')) then
               call speech_input(rtdb)
            else if (inp_compare(.false.,test,'neb')) then
               call neb_input(rtdb)
            else if (inp_compare(.false.,test,'pspw')) then
               call nwpw_input(rtdb)
            else if (inp_compare(.false.,test,'nwpw')) then
               call nwpw_input(rtdb)
            else if (inp_compare(.false.,test,'smd')) then
               call smd_input(rtdb)
            else if (inp_compare(.false.,test,'pause')) then
               call util_sleep(60)
            else if (inp_compare(.false.,test,'rism')) then
               call rism_input(rtdb)
            else if (inp_compare(.false.,test,'qmmm')) then
               call qmmm_input(rtdb)
            else if (inp_compare(.false.,test,'ccca')) then
               call ccca_input(rtdb)
c            else if (inp_compare(.false.,test,'embed')) then
c               call embed_input(rtdb)
            else if (inp_compare(.false.,test,'relativistic')) then
               call rel_input(rtdb)
            else if (inp_compare(.false.,test,'nbo')) then
               call nbo_input(rtdb)
            else if (inp_compare(.false.,test,'vscf')) then
               call vscf_input(rtdb)
            else if (inp_compare(.false.,test,'raman')) then
               call raman_input(rtdb)
            else if (inp_compare(.false.,test,'dntmc')) then
               call dntmc_input(rtdb)
            else if (inp_compare(.false.,test,'vib')) then
               call freq_vib_input(rtdb)
* this reordered test will match freq and freq*
            else if (inp_compare(.false.,'freq',test)) then  
               call freq_vib_input(rtdb)
            else if (inp_compare(.false.,test,'hessian')) then
               call hess_input(rtdb)
            else if (inp_compare(.false.,test,'tddft')) then
               call tddft_input(rtdb)
            else if (inp_compare(.false.,test,'rt_tddft')) then
               call rt_tddft_input(rtdb)
            else if (inp_compare(.false.,test,'fcidump')) then
               call fcidump_input(rtdb)
            else if (inp_compare(.false.,test,'mymd')) then
                call mymd_input(rtdb)
            else if (inp_compare(.false.,test,'mymc')) then
                call mymc_input(rtdb)
            else if (inp_compare(.false.,test,'tce')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'tce'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ucc')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'ucc'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'uci')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'uci'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ump')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'umbpt')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'uccd')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'uccd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ulccd')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'ulccd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'uccsd')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'uccsd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ulccsd')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'ulccsd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'uqcisd')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'uqcisd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'uccsdt')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'uccsdt'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'uccsdtq')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'uccsdtq'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ucisd')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'ucisd'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ucisdt')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'ucisdt'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ucisdtq')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'ucisdtq'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ump2')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt2'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ump3')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt3'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'ump4')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt4'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'umbpt2')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt2'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'umbpt3')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt3'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)
            else if (inp_compare(.false.,test,'umbpt4')) then
               if (.not.rtdb_cput(rtdb,'tce:module',1,'umbpt4'))
     1         call errquit('input_parse: failed writing to rtdb',0,
     &       RTDB_ERR)
               call tce_input(rtdb)

#ifdef MRCC_METHODS
c#warning "MRCC methods are currently under development"
ckbn mrcc-r-1 -3
            else if (inp_compare(.false.,test,'mrccdata')) then
               call tce_mrcc_input(rtdb)
#endif
            else if (inp_compare(.false.,test,'gw')) then
               call gw_input(rtdb)
            else if (inp_compare(.false.,test,'xtb')) then
               call xtb_input(rtdb)
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
               call errquit('input_parse: unknown directive '//test,0,
     &       INPUT_ERR)
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
c            write(LuOut,*) 
c            write(LuOut,*) ' ! warning: processed input with no task'
c            write(LuOut,*)
c            call util_flush(LuOut)
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
      subroutine input_qcharge(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      double precision qcharge
      character*255 test
c
      if (.not. inp_a(test)) call errquit
     &    ('input_qcharge: failed to read ascii string',911, INPUT_ERR)
      if (.not.inp_compare(.false.,test,'qcharge')) call errquit
     &    ('input_qcharge: failed to read "qcharge"',911, INPUT_ERR)
c
      if (.not. inp_f(qcharge)) call errquit
     $     ('input_qcharge: failed to read real number', 0, INPUT_ERR)
c
      if (.not. rtdb_put(rtdb, 'qcharge', MT_DBL, 1, qcharge))
     $     call errquit
     $     ('input_qcharge: failed to write to rtdb', 0, RTDB_ERR)
c
      end
      subroutine input_charge(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      double precision charge
      character*255 test
c
      if (.not. inp_a(test)) call errquit
     &    ('input_charge: failed to read ascii string',911, INPUT_ERR)
      if (.not.inp_compare(.false.,test,'charge')) call errquit
     &    ('input_charge: failed to read "charge"',911, INPUT_ERR)
c
      if (.not. inp_f(charge)) call errquit
     $     ('input_charge: failed to read real number', 0, INPUT_ERR)
c
      if (.not. rtdb_put(rtdb, 'charge', MT_DBL, 1, charge))
     $     call errquit
     $     ('input_charge: failed to write to rtdb', 0, RTDB_ERR)
c
      end

