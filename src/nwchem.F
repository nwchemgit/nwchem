      program nwchem
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "pstat.fh"
#include "util.fh"
#include "inp.fh"
      integer rtdb
      integer stack
      integer heap
      integer global
      logical status
      logical overify
#ifdef PARALLEL_DIAG
      integer junk(4)
#endif
c
c     $Id: nwchem.F,v 1.77 1996-05-13 16:21:53 gg502 Exp $
*------------------------------------------------*
* COPYRIGHT (C) 1994, 1995, 1996,                *
*         Pacific Northwest National Laboratory, * 
*         Battelle Memorial Institute.           *
*------------------------------------------------*
*------------> All Rights Reserved <-------------*
*------------------------------------------------*
c
*--------------------------------------------------------------------*
*                             DISCLAIMER                             *
*--------------------------------------------------------------------*
* This material was prepared as an account of work sponsored by an   *
* agency of the United States Government.  Neither the United States *
* Government nor the United States Department of Energy, nor         *
* Battelle, nor any of their employees, MAKES ANY WARRANTY, EXPRESS  *
* OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR   *
* THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION,      *
* APPARATUS, PRODUCT, SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS  *
* THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.            *
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
*                            LIMITED USE                             *
*--------------------------------------------------------------------*
* This software (including any documentation) is being made          *
* available to you for your internal use only, solely for use in     *
* performance of work directly for the U.S. Federal Government or    *
* work under contracts with the U.S. Department of Energy or other   *
* U.S. Federal Government agencies.  This software is a version      *
* which has not yet been evaluated and cleared for commercialization.*
* Adherence to this notice may be necessary for the author, Battelle *
* Memorial Institute, to successfully assert copyright in and        *
* commercialize this software.  This software is not intended for    * 
* duplication or distribution to third parties without the           *
* permission of the Manager of Software Products at Pacific          *
* Northwest National Laboratory, Richland, Washington, 99352.        *
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
*                            ACKNOWLEDGMENT                          *
*--------------------------------------------------------------------*
* This software and its documentation were produced with Government  *
* support under Contract Number DE-AC06-76RLO-1830 awarded by the    *
* United States Department of Energy.  The Government retains a      *
* paid-up non-exclusive, irrevocable worldwide license to reproduce, *
* prepare derivative works, perform publicly and display publicly by *
* or for the Government, including the right to distribute to other  *
* Government contractors.                                            *
*--------------------------------------------------------------------*
c     
c
      character*255 input_filename
      double precision total_wall, total_cpu
c
c     Create parallel processes and initialize IPC layer
c
      call pbeginf()
#ifdef PARALLEL_DIAG
      call mxpara(junk)
#endif
c
c     Initialize timers so they are relative to job start
c
**      call ieeetrap()
      total_wall = -util_wallsec()
      total_cpu  = -util_cpusec()
c
c     Only process 0 opens the input file 
c     (note that ga_nodeid() does not work yet!)
c
      if (nodeid() .eq. 0) then   
c
c     Get the name of the input file from the command line
c
         call get_input_filename(input_filename)
c
         open(5, file=input_filename, status='old', form='formatted',
     $        err=10000)
         goto 10001
10000    call errquit('nwchem: failed to open the input file', 0)
10001    continue
c
      else
#ifndef KSR
         close(5)
#endif
      endif
c
c     Look for memory directive in the input ... must eventually
c     open the rtdb first so that can get memory directive out of that
c     if it is not in the input
c
C     The user input model has well-defined categories of memory,
C     each of which has a specific size.  How we use these limits
C     depends on the platform we are running on.
C
      call input_mem_size(stack, heap, global, overify)
c
c     Initialize local memory allocator & global array tools
C
C     GA allocations come out of MA space, so lump them together
C     and let MA impose the limits on GA sizes instead of actually
C     using the global limit.
C
      if ( ga_uses_ma() ) then
         if (.not. ma_init(mt_dbl, stack, heap+global))
     $      call errquit('nwchem: ma_init failed', -1)
         call ga_initialize
C
C     GA allocations are separate from MA, so the separate limit
C     must be enforced.  Note GA only understands bytes.
C
      else
         if (.not. ma_init(mt_dbl, stack, heap))
     $      call errquit('nwchem: ma_init failed', -1)
         call ga_initialize_ltd(ma_sizeof(mt_dbl, global, mt_byte) )
      endif
C
      status = ma_set_auto_verify(overify)
      status = ma_set_hard_fail(.true.)
c
c     Initialize pstat
c
      if (.not. pstat_init(20, 1, ' ')) call errquit
     $     ('nwchem: pstat_init failed', 0)
c
c     Now are ready to summarize the environment
c
      call nwchem_banner(input_filename)
c
c     Read the input and open the RTDB
c
      call start(rtdb)
c
      if (util_print('debug tcgmsg',print_never)) 
     $     call setdbg(1)
c
c     Call driver routine to do something useful
c
      call control(rtdb)
c
c     Close the RTDB
c
      if (util_print('rtdbvalues', print_debug)) then
         if (.not. rtdb_print(rtdb, .true.))
     $     call errquit('control: rtdb_print failed', 0)
      else if (util_print('rtdb', print_high)) then
         if (.not. rtdb_print(rtdb, .false.))
     $     call errquit('control: rtdb_print failed', 0)
      endif
c
      if (.not. rtdb_close(rtdb, 'keep'))
     $     call errquit('nwchem: rtdb_close failed', rtdb)
c
      if (util_print('rtdb', print_high) .or.
     $    util_print('rtdbvalues', print_high)) then
         call rtdb_print_usage  ! Called after closing so memory leaks apparent
      endif
c
c     Tidy up pstat
c
      if (.not. pstat_terminate()) call errquit
     $     ('nwchem: pstat_terminate failed', 0)
c
c     Print memory and other info
c
      call ga_sync()
      if (ga_nodeid() .eq. 0) then
         if (util_print('ga summary', print_default))
     $        call ga_summarize(0)
         if (util_print('ma summary', print_default))
     $        call ma_summarize_allocated_blocks()
         if (util_print('ma stats', print_high))
     $        call ma_print_stats()
         if (util_print('version', print_debug))
     $        call util_version
         if (util_print('citation', print_low))
     &       call util_cite()
         if (util_print('total time', print_low)) then
            total_wall = total_wall + util_wallsec()
            total_cpu  = total_cpu  + util_cpusec()
            write(6,1) total_cpu, total_wall
 1          format(/' Total times  cpu:',f11.1,'s     wall:',f11.1,'s')
            call util_flush(6)
         endif
      endif
      call ga_sync()
c
c     Tidy up the global memory
c
      call ga_terminate()
c
c     Tidy up the parallel environment
c
      call pend()
c
c     Done
c
      end
      subroutine start(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
      integer rtdb              ! [output]
      character*255 rtdb_name
      character*255 file_prefix
      logical startup
c
c     Figure out if this is a start/restart and what the database
c     and file prefix are
c
      call input_start_opt(rtdb_name, file_prefix, startup)
c
c     Actually open the database and insert the file_prefix
c
c     Note that only process 0 has the database name ... that is OK.
c
      if (startup) then
         if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $        call errquit('start: rtdb_open empty failed ', 0)
      else
         if (.not. rtdb_open(rtdb_name, 'old', rtdb))
     $        call errquit('start: rtdb_open old failed ', 0)
      endif
      call util_file_prefix_set(file_prefix)
      call util_file_prefix_rtdb_store(rtdb) ! Not ideal location for this
c
      if (ga_nodeid() .eq. 0) then
         write(6,1) file_prefix(1:inp_strlen(file_prefix))
 1       format(' file prefix = ', a//)
         call util_flush(6)
      endif
c
c     Parse input data and shove into the database
c
      call input_parse(rtdb)
c
c     Extract any high level print options
c
      call util_print_rtdb_load(rtdb,' ')
c
c     Done
c
      end
      subroutine control(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb              ! [input] data base handle
c     
c     Drive the calculation.  Eventually this routine will hook up with the
c     yet-to-be-written task module.  Right now task_input makes an rtdb
c     entry task that comprises just the name of the routine to be called.
c     If no task is present then perform some default action (exit).
c     
      character*255 task_name
      integer ntasks, ind
      logical status 
      logical  scf, task_check_sum, selci, eopt, rjhtest
      external scf, task_check_sum, selci, eopt, rjhtest
      logical rimp2, oimp2
      external rimp2, oimp2
      logical nwdft
      external nwdft
      logical nwargos
      external nwargos
      logical mp2_grad
      external mp2_grad
      logical cphf
      external cphf
      parameter (ntasks = 26)
      character*10 task_list(ntasks)
      data task_list /'rjhtest', ! 1
     $     'raktest',           ! 2
     $     'rtdbprint',         ! 3
     $     'scf',               ! 4
     $     'oimp2',             ! 5
     $     'gradients',         ! 6
     $     'optscf',            ! 7
     $     'atwtest',           ! 8
     $     'rimp2',             ! 9
     $     'checksum',          ! 10
     $     'dft',               ! 11
     $     'md',                ! 12
     $     'dftgrad',           ! 13
     $     'optdft',            ! 14
     $     'direct_mp2',        ! 15
     $     'cphf',           ! 16
     $     'ccsd',              ! 17
     $     'freqscf',           ! 18
     $     'freqdft',           ! 19
     $     'mcscf',             ! 20
     $     'nwargos',           ! 21
     $     'fourindex',         ! 22
     $     'selci',             ! 23
     $     'dplot',             ! 24
     $     'mp2_grad',          ! 25
     $     'optimize'           ! 26
     $     /
c     
      if (.not. rtdb_cget(rtdb, 'task', 1, task_name)) then
c     
c     Default action ... do nothing
c     
         return
      endif
c     
      if (inp_match(ntasks, .false., task_name, task_list, ind)) then
         goto (10,  20,  30,  40,  50,  60,  70,  80,  90, 100, 
     $        110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
     $        210, 220, 230, 240, 250, 260) ind
c     
         write(6,*) ' Unimplemented task ', 
     $        task_name(1:inp_strlen(task_name))
         goto 1000
c     
c     rjhtest
c     
 10      status = rjhtest(rtdb)
         goto 1000
c     
c     raktest
c     
 20      call raktest(rtdb)
         goto 1000
c     
c     rtdbprint
c     
 30      if (.not. rtdb_print(rtdb, .true.))
     $        call errquit('control: rtdb_print failed', 0)
         goto 1000
c     
c     scf ... eventually there will be no such task name
c     
 40      status = scf(rtdb)
         goto 1000
C     
C     OI-MP2 = Orbital Invariant MP2 (different from RI-MP2)
C     
 50      status = oimp2(rtdb)
         goto 1000
c     
c     gradients (scf has to be done first)
c     
 60      call scf_gradient(rtdb)
         goto 1000
c     
c     optscf ... scf optimization
c     
 70      call scf_optimize(rtdb)
         goto 1000
c     
c     atwtest (4-index)
c     
 80      call atwtest(rtdb)
         goto 1000
c     
C     RI-MP2 energy (SCF must be done separately first)
c     
 90      status = rimp2(rtdb)
         goto 1000
c     
c     integral check sum
c     
 100     status = task_check_sum(rtdb)
         goto 1000
c     
c     dft
c     
 110     continue 
         status = nwdft(rtdb)
         goto 1000
c     
c     MD
c     
 120     continue
         call mdmain(rtdb)
         goto 1000
c     
c     dft gradient
c     
 130     continue
C     write(6,*) 'NWdft gradient not linked into this executable'
         call dft_grad_top(rtdb)
         goto 1000
c     
c     optdft ... dft optimization
c     
 140     call dft_optimize(rtdb)
         goto 1000
c     
c     conventional direct MP2 (requires moints module)
c     
 150     call direct_mp2(rtdb)
         goto 1000
c     
c     Miles test 
c     
 160     status=cphf(rtdb)
         goto 1000
c     
c     CCSD
c     
 170     call ccsd(rtdb)
         goto 1000
c     
c     SCF frequencies by finite difference (in stepper)
c     
 180     call scf_freq(rtdb)
         goto 1000
c     
c     DFT frequencies by finite difference (in stepper)
c     
 190     call dft_freq(rtdb)
         goto 1000
c     
c     MCSCF 
c     
 200     call mcscf(rtdb)
         goto 1000
c     
c     nw_Argos
c     
  210 continue
      status=nwargos(rtdb)
      goto 1000
c     
c     Four index
c     
 220     call moints_full(rtdb)
         goto 1000
c
c     Selected CI
c
 230     status = selci(rtdb)
         goto 1000
c     
c     Plot of the density matrix
c     
 240     call DPlot(rtdb)
         goto 1000
c
c     semi-direct mp2 energy and gradient module
c     
 250     status=mp2_grad(rtdb)
         goto 1000
c
c     soon to be generic optimization interface
c     currently just a hook for the energy optimization
c
 260     status = eopt(rtdb)
         goto 1000
c
      else
c     
         write(6,*) ' task_name = ', 
     $        task_name(1:inp_strlen(task_name))
         call errquit('control: unknown task name', 0)
      endif
c     
c     Eventually this will be where the next task is identified
c     
 1000 continue
c     
c     Finished ... eventually need to track status of tasks so that
c     restarts happen correctly.  For now assume that the task completed
c     and delete task from the database so that a restart will have no
c     effect.
c     
      status = rtdb_delete(rtdb, 'task') ! quietly discard the status
c     
      end
      subroutine nwchem_banner(input_filename)
      implicit none
#include "global.fh"
#include "inp.fh"
#include "mafdecls.fh"
c
c     Print a banner at the top of the output to identify the origin
c     of the code along with hostname, no. of processors, date, ...
c
      character*(*) input_filename
      character*80 hostname, compiled, srcdir
      character*256 executable
      character*26 date
#if defined(CRAY)
      integer ilen,ierror
#endif
      integer heap, stack, global, global_b, total, nproc
      logical status
      character*3 cstatus
      Character*40 GStatus
c
      if (ga_nodeid() .ne. 0) return
c
      call util_print_centered(6,
     $     '                                                          ',
     $     40, .true.)
      write(6,*)
      write(6,*)
      call util_print_centered(6,
     $     'Northwest Computational Chemisty Package (NWChem)',
     $     40, .true.)
      write(6,*)
      write(6,*)
      call util_print_centered(6,
     $      'Environmental Molecular Sciences Laboratory',
     $     40, .false.)
      call util_print_centered(6,
     $     'Pacific Northwest National Laboratory',
     $     40, .false.)
      call util_print_centered(6,'Richland, WA 99352',
     $     40, .false.)
      write(6,*)
      call util_print_centered(6,
     $     '                                                          ',
     $     40, .true.)
      write(6,*)
      write(6,*)
      call util_legal()
c
#if defined(CRAY)
      call pxfgetarg(0,executable,ilen,ierror)
#else
      call getarg(0,executable)
#endif
      call util_hostname(hostname)
      call util_date(date)
      nproc = ga_nnodes()
c
      compiled = COMPILATION_DATE
      srcdir   = 
     $ COMPILATION_DIR
      call util_print_centered(6, 'Job information', 20, .true.)
      write(6,2) 
     $     hostname(1:inp_strlen(hostname)),
     $     executable(1:inp_strlen(executable)),
     $     date(1:inp_strlen(date)),
     $     compiled(1:inp_strlen(compiled)),
     $     srcdir(1:inp_strlen(srcdir)),
     $     input_filename(1:inp_strlen(input_filename)),
     $     nproc
 2    format(/
     $     '    hostname = ', a/
     $     '    program  = ', a/
     $     '    date     = ', a/
     $     '    compiled = ', a/
     $     '    source   = ', a/
     $     '    input    = ', a/
     $     '    nproc    = ', i4//)
c
c     Determine, without altering, setting of memory verification
c
      status = ma_set_auto_verify(.true.)
      if (status) then
         cstatus = 'yes'
      else
         cstatus = 'no '
      endif
      status = ma_set_auto_verify(status)
C
C     Memory size information.  Note:  Subtleties of MA & friends
C     will likely result in actual heap and stack sizes _slightly_ 
C     larger than requested at initialization, and other slight
C     "inconsistencies".  It has to do with the overhead associated
C     with an allocation.
C
      heap  = ma_inquire_heap(MT_Dbl)
      stack = ma_inquire_stack(MT_Dbl)
C
C     If GA cannot determine a limit on the available memory, it
C     will return a negative number, which MA_SizeOf does not
C     appreciate.  This should not happen in NWChem, but we might
C     as well play it safe.
C
      global_b = ga_memory_avail()
      If ( Global_B .ge. 0) then
         global = MA_SizeOf(MT_Byte, global_b, MT_Dbl)
      Else
         global = 0
      EndIf
C
      Total = Heap + Stack
      If ( .NOT. GA_Uses_MA() ) Total = Total + Global
C
      If ( GA_Uses_MA() ) then
         GStatus = '(allocated from within heap+stack)'
      Else
         GStatus = '(allocated separately from heap & stack)'
      EndIf
c
      call util_print_centered(6, 'Memory information', 20, .true.)
c
      write(6,3) heap, stack, global, GStatus, total, cstatus
 3    format(/
     $     '    heap   = ', i10,' doubles'/
     $     '    stack  = ', i10,' doubles'/
     $     '    global = ', i10,' doubles ', A/
     $     '    total  = ', i10,' doubles'/
     $     '    verify = ',a3//)
c
      call util_flush(6)
c
      end
