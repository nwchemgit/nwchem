      program nwchem
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#ifdef USE_TCGMSG
#include "tcgmsg.fh"
#else
      integer nodeid
      external nodeid
#endif
#include "pstat.fh"
#include "util.fh"
#include "inp.fh"
#include "bgj_common.fh"
#include "stdio.fh"
      integer rtdb
      integer stack
      integer heap
      integer global
      logical status
      logical overify, ohardfail
#ifdef PSCALE
      integer io_code
#else
      integer*4 io_code
#endif
c
c  $Id$
c
c  ======================================================================================================
C> \mainpage Northwest Computational Chemistry Package (NWChem) 7.2.0
C>
C> NWChem is an open-source computational chemistry package distributed under the terms of
C> the Educational Community License (ECL) 2.0
C>
C> This software and its documentation were developed at the EMSL at Pacific Northwest National Laboratory,
C> a multiprogram national laboratory, operated for the U.S. Department of Energy by Battelle under
C> Contract Number DE-AC05-76RL01830. Support for this work was provided by the Department of Energy
C> Office of Biological and Environmental Research, Office of Basic Energy Science, and the Office of
C> Advanced Scientific Computing.
C>
C> Licensed under the Educational Community License, Version 2.0 (the "License"); you may
C> not use this file except in compliance with the License. You may obtain a copy of the
C> License at <a href="https://opensource.org/licenses/ECL-2.0">https://opensource.org/licenses/ECL-2.0</a>.
C>
C> Unless required by applicable law or agreed to in writing, software distributed under the
C> License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
C> either express or implied. See the License for the specific language governing
C> permissions and limitations under the License.
C>
C> Further information, including user documentation and forums, may be found at
C> <a href="https://nwchemgit.github.io">https:/nwchemgit.github.io</a>. Alternatively,
C> the paper
C>
C> * M. Valiev, E.J. Bylaska, N. Govind, K. Kowalski, T.P. Straatsma, H.J.J. Van Dam,
C>   D. Wang, J. Nieplocha, E. Apra, T.L. Windus, W.A. de Jong (2010)<br>
C>   "NWChem: A comprehensive and scalable open-source solution for large scale molecular simulations"<br>
C>   <I>Computer Physics Communications</I>, <b>181</b>, 1477–1489, DOI: <a href="https://doi.org/10.1016/j.cpc.2010.04.018">10.1016/j.cpc.2010.04.018</a>
C>
C> provides details on the codes capabilities.
C>
C> Copyright (c) 1994-2022 Pacific Northwest National Laboratory, Battelle Memorial Institute
C>
C> Environmental Molecular Sciences Laboratory (EMSL)<br>
C> Pacific Northwest National Laboratory<br>
C> Richland, WA 99352
c
c  ======================================================================================================
c
      character*(nw_max_path_len) input_filename, rtdb_name
      double precision total_wall, total_cpu
#ifdef USE_OFFLOAD
      integer ppnout
      logical offload_enabled
      external offload_enabled
      integer offload_device
      external offload_device
#endif

      logical ostartup, ocontinue, orestart
      logical  input_parse
      external input_parse
#if defined(USE_OPENMP)
      integer omp_get_max_threads
      external omp_get_max_threads
#endif
      logical do_printaff
      integer linux_printaff
c
c     Create parallel processes and initialize IPC layer
c
      call pbeginf()

c
c     MXINIT is needed by PeIGS and PFFT to initialize
c     the communication fabric they use.
c
      call mxinit
c
c     Initialize timers so they are relative to job start
c
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
         open(LuIn, file=input_filename, status='old', form='formatted',
     $        err=10000)
         goto 10001
10000    call errquit('nwchem: failed to open the input file', 0,
     &       INPUT_ERR)
10001    continue
c
      else
#if !(defined(KSR) || defined(IBM) || defined(FUJITSU_SOLARIS) )
         close(LuIn)
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
      call input_mem_size(stack, heap, global, overify, ohardfail)
c
c     Initialize local memory allocator & global array tools
C
      call ga_initialize_ltd(ma_sizeof(mt_dbl,global,mt_byte))
      ! this must happen after GA and before MA!
      call util_setup_gpu_affinity()
      if ( ga_uses_ma() ) then
        if (.not.ma_init(mt_dbl, stack, heap+global))
     &      call errquit('nwchem.F: ma_init failed (ga_uses_ma=T)',911,
     &      MA_ERR)
      else
        if (.not.ma_init(mt_dbl,stack,heap))
     &      call errquit('nwchem.F: ma_init failed (ga_uses_ma=F)',911,
     &      MA_ERR)
      endif
!
!     Touch OpenMP here so that any runtime initialization happens up-front.
!     This ensures that any printout that the OpenMP runtime generates,
!     such as affinity information, appears at the top of the output file.
!     Otherwise, it might not appear until e.g. the CCSD module, at which
!     point it will pollute the output file in an undesirable way.
!
!     Do not move this in front of GA/MPI/TCGMSG initialization, since the
!     OpenMP runtime may inherit affinity information from MPI that is only
!     determined during MPI initialization.
!
!     Format definition is outside of preprocessor protection to ensure the
!     label is not accidentally reused, since that will not be caught by
!     testing that does not enable OpenMP.
!
 99   format(2x,'NWChem w/ OpenMP: maximum threads = ',i4)
      do_printaff=.false.
#if defined(USE_OPENMP)
      !$omp parallel
      !$omp master
      if (ga_nodeid().eq.0) write(luout,99) omp_get_max_threads()
      !$omp end master
      !$omp end parallel
#endif
c
c     set no. threads for threaded BLAS to 1
      call util_blas_set_num_threads(1)
!
      call rtdb_init()
c
*
*     More for amusement then efficiency force all MA allocated entities 
*     to be aligned at the beginning of a 128 byte cache line
*
!      if (.not. ma_set_numalign(7))
!     $     call errquit('nwchem.F: ma_set_numalign failed',911, MA_ERR)
c aligned to 64byte record
      if (.not. ma_set_numalign(6))
     $     call errquit('nwchem.F: ma_set_numalign failed',911, MA_ERR)
*
*
*old:------------------------------------------------------- START ---------
*old:C     GA allocations come out of MA space, so lump them together
*old:C     and let MA impose the limits on GA sizes instead of actually
*old:C     using the global limit.
*old:C
*old:      if ( ga_uses_ma() ) then
*old:         if (.not. ma_init(mt_dbl, stack, heap+global))
*old:     $      call errquit('nwchem: ma_init failed', -1)
*old:         call ga_initialize
*old:C
*old:C     GA allocations are separate from MA, so the separate limit
*old:C     must be enforced.  Note GA only understands bytes.
*old:C
*old:      else
*old:         if (.not. ma_init(mt_dbl, stack, heap))
*old:     $      call errquit('nwchem: ma_init failed', -1)
*old:         call ga_initialize_ltd(ma_sizeof(mt_dbl, global, mt_byte) )
*old:      endif
*old:------------------------------------------------------- END -----------
c
*****      call nxtval_ga_initialize()
c
c     Trap SIGFPE after GA to override handler
c
*****      call ieeetrap()
#if defined(LINUXALPHA)
      call dec_fpe  ! To avoid underflow problems on Alpha in Texas
#endif
#ifdef LINUX
c      uncommenting this line turns on sigfpe trapping under linux
!      call linux_trapfpe()
#endif
#if defined(MACX)
c      uncommenting this line turns on sigfpe trapping under Mac OSX
!      call macx_trapfpe()
#endif
c
c     Hard fail is good for development but means that we cannot
c     respond to allocation problems. Disable by default.
c
      status = ma_set_auto_verify(overify)
      status = ma_set_hard_fail(ohardfail)
      status = ma_set_error_print(ohardfail)
c
c     Initialize pstat
c
      if (.not. pstat_init(20, 1, ' ')) call errquit
     $     ('nwchem: pstat_init failed', 0, UNKNOWN_ERR)
c
      call input_file_info(input_filename, rtdb_name,
     &    ostartup, ocontinue)
c
c     Now are ready to summarize the environment
c
      call nwchem_banner(input_filename, rtdb_name,
     &    ostartup, ocontinue)
c
c     Actually open the database and store the file prefix 
c
c     Note that only process 0 has the database name ... that is OK.
c
      if (ostartup) then
         if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $        call errquit('start: rtdb_open empty failed ', 0,
     &       RTDB_ERR)
      else
         if (.not. rtdb_open(rtdb_name, 'old', rtdb))
     $        call errquit('start: rtdb_open old failed ', 0, RTDB_ERR)
      endif
C
C     initialize nxtask
C
      call nxtask_init(rtdb)
c!!! BGJ
      bgj_rtdb = rtdb
c!!! BGJ
c
c
      
      if (ostartup.or.ocontinue) then
        orestart = .false.
      else
        orestart = .true.
      endif
c
      call util_set_rtdb_state(rtdb,ostartup,ocontinue,orestart)
c
      call util_file_info_rtdb(rtdb) ! Save file info for restart
c
      call movecs_ecce_print_on()
c
      call geom_hnd_parallel(.true.)
c
      call perfm_start()

#ifdef USE_OFFLOAD
      call util_getppn(ppnout)
      if(ppnout.eq.0) call errquit('util_getppn failed ',0,UERR)
      if(ga_nodeid().eq.0)write(luout,*)  ga_nodeid(),' ppn ',ppnout
      if(offload_enabled()) then
      if(ga_nodeid().lt.ppnout)
     w        write(luout,'(I8,A,I2)')  ga_nodeid(),
     w              ' offload enabled, GPU: ', offload_device()
      endif
      call ga_sync()
#endif

      if (orestart .or. ocontinue) call nw_print_restart_info(rtdb)
c
c if continue then go right to task stored on rtdb do not further parse
c input.  If input is required then user should have used restart
c
      if (ocontinue) call task(rtdb) 
c
c
c     Parse input data, shove into the database and execute tasks
c
 10   if (input_parse(rtdb)) then ! while (tasks to do)
         call util_print_rtdb_load(rtdb,' ') ! High level print
         if (util_print('tcgmsg',print_never)) then 
            call setdbg(1)
         else
            call setdbg(0)
         endif
#ifdef CATAMOUNT
         call util_allocga()
#endif
      if (.not.rtdb_get(rtdb,'do_printaff',mt_log,1,do_printaff))
     &     do_printaff = .false.
      if(do_printaff) then
         if(ga_nodeid().eq.0)
     W        write(luout,*) '%%%%%% CPU Affinity %%%%%%%'
         call util_flush(luout)
         call ga_sync()
         if(linux_printaff().ne.0) call errquit(
     A   ' call to linux_printaff failed ',0,CAPMIS_ERR)
      endif
         call task(rtdb)
         goto 10                !    end while
      endif
c
c     Close the RTDB
c
      call util_print_rtdb_load(rtdb,' ') ! High level print
      if (util_print('rtdbvalues', print_debug)) then
         if (.not. rtdb_print(rtdb, .true.))
     $     call errquit('control: rtdb_print failed', 0, RTDB_ERR)
      else if (util_print('rtdb', print_high)) then
         if (.not. rtdb_print(rtdb, .false.))
     $     call errquit('control: rtdb_print failed', 0, RTDB_ERR)
      endif
c
      if (.not. rtdb_close(rtdb, 'keep'))
     $     call errquit('nwchem: rtdb_close failed', rtdb, RTDB_ERR)
c
      if (util_print('rtdb', print_high) .or.
     $    util_print('rtdbvalues', print_high)) then
         call rtdb_print_usage  ! Called after closing so memory leaks apparent
      endif
c
c     Tidy up pstat
c
      if (.not. pstat_terminate()) call errquit
     $     ('nwchem: pstat_terminate failed', 0, UNKNOWN_ERR)
c
****      call nxtval_ga_terminate()
c
c     Print memory and other info
c
      call ga_sync()
      if (ga_nodeid() .eq. 0) then
         if (util_print('ga summary', print_default))
     $        call ga_summarize(0)
         if (util_print('ga stats', print_high)) then
            call ga_print_stats()
            write(LuOut,*)
         endif
         if (util_print('ma summary', print_default))
     $        call ma_summarize_allocated_blocks()
         if (util_print('ma stats', print_high)) then
            call ma_print_stats(.true.)
         else if (util_print('ma stats', print_default)) then
            call ma_print_stats(.false.)
         endif
         if (util_print('version', print_debug))
     $        call util_version
         if (util_print('citation', print_none))
     &       call util_cite()
      endif
      call perfm_end()
      total_wall = total_wall + util_wallsec()
      total_cpu  = total_cpu  + util_cpusec()
      if (ga_nodeid() .eq. 0) then
        if (util_print('total time', print_none)) then
          write(LuOut,1,iostat=io_code) total_cpu, total_wall
 1        format(/' Total times  cpu:',f11.1,'s     wall:',f11.1,'s')
          if(io_code.ne.0.and.
     A         util_print('total time', print_high))  then
             write(luout,*) ' iostat ',io_code
           endif
           call util_flush(LuOut)
        endif
      endif
      call ecce_print1('all tasks cpu time',mt_dbl,total_cpu,1)
      call ecce_print1('all tasks wall time',mt_dbl,total_wall,1)
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
c     Close the ECCE' output file if there is one
c
      call ecce_print_file_close()
c
c     Done
c
      end
      subroutine nwchem_banner(input_filename, rtdb_name, ostartup,
     &    ocontinue)
      implicit none
#include "global.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
c
c     Print a banner at the top of the output to identify the origin
c     of the code along with hostname, no. of processors, date, ...
c
      character*(*) input_filename, rtdb_name
      logical ostartup
      logical ocontinue
c
      double precision dbletomb
      character*80 hostname, compiled, batchid
      character*(nw_max_path_len) executable, srcdir, file_prefix
      character*(nw_max_path_len) branch, nwchem_rev, ga_rev
      character*26 date
#if defined(CRAY)
      integer ilen,ierror
#endif
      integer heap, stack, global, global_b, total, nproc
      logical status
#ifdef CATAMOUNT
      integer istatus,setvbuf3f
      external setvbuf3f
#endif
      logical util_scalapack_info
      external util_scalapack_info
      character*3 cstatus, fstatus
      character*10 cstart
      Character*40 GStatus
c
#ifdef CATAMOUNT
      istatus=setvbuf3f(luout, 0, 4096)
      istatus=setvbuf3f(0, 0, 4096)
#endif
      if (ga_nodeid() .ne.0) goto 10101
c
      call util_print_centered(LuOut,
     $     '                                                          ',
     $     40, .true.)
      write(LuOut,*)
      write(LuOut,*)
      call util_print_centered(LuOut,
     $     'Northwest Computational Chemistry Package (NWChem) 7.2.0',
     $     40, .true.)
      write(LuOut,*)
      write(LuOut,*)
      call util_print_centered(LuOut,
     $      'Environmental Molecular Sciences Laboratory',
     $     40, .false.)
      call util_print_centered(LuOut,
     $     'Pacific Northwest National Laboratory',
     $     40, .false.)
      call util_print_centered(LuOut,'Richland, WA 99352',
     $     40, .false.)
      write(LuOut,*)
      call util_print_centered(LuOut,'Copyright (c) 1994-2022',
     $     40, .false.)
      call util_print_centered(LuOut,
     $     'Pacific Northwest National Laboratory',
     $     40, .false.)
      call util_print_centered(LuOut,'Battelle Memorial Institute',
     $     40, .false.)
      write(LuOut,*)
      call util_print_centered(LuOut,
     $'NWChem is an open-source computational chemistry package'
     $ , 40, .false.)
      call util_print_centered(LuOut,
     $'distributed under the terms of the'
     $ , 40, .false.)
      call util_print_centered(LuOut,
     $'Educational Community License (ECL) 2.0'
     $ , 40, .false.)
      call util_print_centered(LuOut,
     $'A copy of the license is included with this distribution' 
     $ , 40, .false.)
      call util_print_centered(LuOut,
     $'in the LICENSE.TXT file'
     $ , 40, .false.)
      write(LuOut,*)
      call util_legal()
c
      call util_getarg(0,executable)
      call util_hostname(hostname)
      call util_date(date)
      nproc = ga_nnodes()
      if (ostartup) then
         cstart = 'startup'
      else
        if (ocontinue) then
          cstart = 'continue'
        else
          cstart = 'restart'
        endif
      endif
      call util_file_prefix(' ',file_prefix)
c
      compiled =
     & COMPILATION_DATE
      call util_nwchem_srcdir(srcdir)
      srcdir = srcdir(1:max(1,inp_strlen(srcdir)-4))
      branch =
     & NWCHEM_BRANCH   
      call util_nwchem_version(nwchem_rev)
      call util_ga_version(ga_rev)
C
C     See if we can get the batch system id.  We understand NQS and
C     LoadLeveler at present
C
      BatchID = ' '
      Call Util_GetEnv('QSUB_REQID', batchid)
      If ( Inp_StrLen(batchid) .eq. 0 ) then
         Call Util_GetEnv('LOADL_JOB_NAME', batchid)
        If ( inp_strlen(batchid) .eq. 0 ) then
            call util_getenv('LSB_JOBID', batchid)
         endif
      EndIf
C
      call util_print_centered(LuOut, 'Job information', 17, .true.)
      write(LuOut,2) 
     $     hostname(1:inp_strlen(hostname)),
     $     executable(1:inp_strlen(executable)),
     $     date(1:inp_strlen(date)),
     $     compiled(1:inp_strlen(compiled)),
     $     srcdir(1:inp_strlen(srcdir)),
     &     branch(1:inp_strlen(branch)),
     &     nwchem_rev(1:inp_strlen(nwchem_rev)),
     &     ga_rev(1:inp_strlen(ga_rev)),
     S     util_scalapack_info(),
     $     input_filename(1:inp_strlen(input_filename)),
     $     file_prefix(1:inp_strlen(file_prefix)),
     $     rtdb_name(1:inp_strlen(rtdb_name)),
     $     cstart(1:inp_strlen(cstart)),
     $     nproc,
     $     util_batch_job_time_remaining()
C
C  Print version information to ecce.out file
C
      compiled(inp_strlen(compiled)+1:inp_strlen(compiled)+25) = 
     $         ' Version 7.2.0'
      call ecce_print_version(compiled(1:inp_strlen(compiled)))
C
      If ( Inp_StrLen(BatchID) .gt. 0) then
         Write(LuOut,22) BatchID(:Inp_StrLen(BatchID))
      Else
         Write(LuOut,'(//)') 
      EndIf
C
 2    format(/
     $     '    hostname        = ', a/
     $     '    program         = ', a/
     $     '    date            = ', a/
     $     '    compiled        = ', a/
     $     '    source          = ', a/
     $     '    nwchem branch   = ', a/
     $     '    nwchem revision = ', a/
     $     '    ga revision     = ', a/
     $     '    use scalapack   = ', l1/
     $     '    input           = ', a/
     $     '    prefix          = ', a/
     $     '    data base       = ', a/
     $     '    status          = ', a/
     $     '    nproc           = ', i8/,
     $     '    time left       = ', i6,'s')
 22   format(
     $     '    batch job id    = ', a/)
c
c     Determine, without altering, setting of memory verification
c     and hardfail.
c
      status = ma_set_auto_verify(.true.)
      if (status) then
         cstatus = 'yes'
      else
         cstatus = 'no '
      endif
      status = ma_set_auto_verify(status)
c
      status = ma_set_hard_fail(.true.)
      if (status) then
         fstatus = 'yes'
      else
         fstatus = 'no '
      endif
      status = ma_set_hard_fail(status)
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
         GStatus = '(within heap+stack)'
      Else
         GStatus = '(distinct from heap & stack)'
      EndIf
c
      call util_print_centered(LuOut, 'Memory information', 19, .true.)
c
*old:      write(LuOut,3) heap, stack, global, GStatus, total, cstatus,
*old:     $     fstatus
*old: 3    format(/
*old:     $     '    heap      = ', i10,' doubles'/
*old:     $     '    stack     = ', i10,' doubles'/
*old:     $     '    global    = ', i10,' doubles ', A/
*old:     $     '    total     = ', i10,' doubles'/
*old:     $     '    verify    = ', a3/
*old:     $     '    hardfail  = ', a3//)
      dbletomb = 8.0d00/(1024.0d00*1024.0d00)
      write(LuOut,3)
     &    heap,((dble(heap)*dbletomb)),
     &    stack,((dble(stack)*dbletomb)),
     &    global,((dble(global)*dbletomb)),
     &    GStatus(1:inp_strlen(gstatus)),
     &    total,((dble(total)*dbletomb)),
     &    cstatus,
     &    fstatus
 3    format(/
     $     '    heap     = ', i10,' doubles',' = ',f8.1,' Mbytes',/,
     $     '    stack    = ', i10,' doubles',' = ',f8.1,' Mbytes',/,
     $     '    global   = ', i10,' doubles',' = ',f8.1,' Mbytes',
     &     ' ',A,/,
     $     '    total    = ', i10,' doubles',' = ',f8.1,' Mbytes',/,
     $     '    verify   = ', a3,/,
     $     '    hardfail = ', a3,//)
c
      call util_print_centered(LuOut, 
     &     'Directory information', 20, .true.)
      write(LuOut,*)
c
c     Below here back to parallel execution
c
10101 call util_file_print_dirs()
      if(ga_nodeid().eq.0) write(LuOut,*)
c
      end
c
      subroutine nw_print_restart_info(rtdb)
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb
c
c     Summarize the status of a restart calculation
c
      logical mode, task_qmmm, ignore, lstatus
      character*30 operation, status
      character*32 theory
c
      if (ga_nodeid() .eq. 0) then
         mode = rtdb_parallel(.false.)
         if (.not.rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     &        task_qmmm = .false.
c      
         if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $        operation = ' '
c
         if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     $        theory = ' '
c
         if (.not. rtdb_get(rtdb, 'task:ignore', mt_log, 1, ignore))
     $        ignore = .false.
c
         if (rtdb_get(rtdb, 'task:status', mt_log, 1, lstatus)) then
            if (lstatus) then
               status = 'ok'
            else
               status = 'fail'
            endif
         else
            status = 'unknown'
         endif
c
         if (theory.ne.' ' .or. operation.ne.' ' .or. task_qmmm) then
            call util_print_centered(LuOut,'Previous task information',
     $           22, .true.)
            write(LuOut,1) theory, operation, status, task_qmmm, ignore
 1          format(/,
     $           '    Theory    = ', a,/,
     $           '    Operation = ', a,/,
     $           '    Status    = ', a,/,
     $           '    Qmmm      = ', l1,/,
     $           '    Ignore    = ', l1,/)
         endif
         call geom_print_known_geoms(rtdb)
         call bas_print_known_bases(rtdb)
         call util_flush(LuOut)
         mode = rtdb_parallel(mode)
      endif
c
      call ga_sync()
c
      end
      
