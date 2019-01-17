c
c     $Id$
c
      subroutine perfm_start()
      implicit none
#include "perfm.fh"
#ifdef USE_PAPI
      integer ii
      integer*4 arg1
#endif
      double precision util_wallsec,util_cpusec
      external util_wallsec,util_cpusec
#include "yflop.fh"
#ifdef FDIST
#include "dist.fh"
      call initialize_dist()
#endif
#ifdef USE_PAPI
      lfailed=.false.
      retval = PAPI_VER_CURRENT
      call PAPIf_library_init(retval)
      if ( retval.NE.PAPI_VER_CURRENT) then
        write(luout,*) ga_nodeid(),' *** FAILED papi_init '
       call perfm_error(retval,.false.,
     $ ' FAILED PAPI_INIT ')
        lfailed=.true.
        return
      endif
#ifdef PAPI_FPINS
      arg1=PAPI_FP_INS
      call PAPIf_query_event(arg1, retval)
      if (retval .NE. PAPI_OK) call 
     $ perfm_error(retval,.true.,
     $ ' perfm_start: could not query papi FP_INS event ')

      call PAPIf_set_debug(PAPI_VERB_ECONT, retval)
      if (retval .NE. PAPI_OK) call 
     $ perfm_error(retval,.true.,
     $ ' perfm_start: could not query papi setdebug ')
      uso=util_wallsec()
      proc0=util_cpusec()
      call PAPIf_flips( real_time, proc_time, flpins, mflops, retval )
      if (retval .NE. PAPI_OK) call 
     $ perfm_error(retval,.true.,
     $ ' perfm_start: could not init papi flpins ')
#else
      call PAPIf_create_eventset(es, retval)
      if (retval .NE. PAPI_OK) call errquit
     $ (' perfm_start: could not papi eventset ',retval,0)
!      eventname(1)='FP_COMP_OPS_EXE'
      nevents=2
      eventname(1)='PAPI_FP_INS'
      eventname(2)='PAPI_FP_INS'
      do ii=1,nevents
         write(0,*) ' adding event no. ',ii
      call PAPIf_event_name_to_code(eventname(ii),
     $     nativees, retval )
      if (retval .NE. PAPI_OK) call
     $ perfm_error(retval,.true.,
     $     'perfm_start: could not papi event2name ')
      call PAPIf_add_event( es, nativees, retval )
      if (retval .NE. PAPI_OK) call 
     $ perfm_error(retval,.true.,
     $     'perfm_start: could not papi addevent ')
      enddo
c
      call PAPIf_start(es, retval)
      uso=util_wallsec()
      proc0=util_cpusec()
      if (retval .NE. PAPI_OK) call errquit
     $ (' perfm_start: could not papi start ',retval,0)
#endif
#elif USE_YFLOP
c     initialize
      flops_ycount=0
#endif
      return
      end
#ifdef USE_PAPI
c
c
c
      subroutine perfm_error(retval,lstop,where)
      implicit none
#include "perfm.fh"
      character*(*) where
      logical lstop
         call PAPIF_perror(retval,papi_errstr,iret)
         write(*,*) 'PAPI error: ',
     ,        papi_errstr(1:inp_strlen(papi_errstr))
         if(lstop)
     A        call errquit(where,retval,0)
         return
         end
#endif
