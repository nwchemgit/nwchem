c
c     $Id: perfm_start.F,v 1.1 2004-10-11 19:39:05 edo Exp $
c
      subroutine perfm_start()
      implicit none
#include "perfm.fh"
#ifdef USE_PAPI
      integer ii
      retval = PAPI_VER_CURRENT
      call PAPIf_library_init(retval)
      if ( retval.NE.PAPI_VER_CURRENT) call errquit
     $ (' perfm_start: wrong papi header from ver ',retval,0)
#ifdef PAPI_FPINS
      call PAPIf_query_event(PAPI_FP_INS, retval)
      if (retval .NE. PAPI_OK) call errquit
     $ (' perfm_start: could not query papi FP_INS event ',retval,0)

      call PAPIf_set_debug(PAPI_VERB_ECONT, retval)
      if (retval .NE. PAPI_OK) call errquit
     $ (' perfm_start: could not query papi setdebug ',retval,0)

      call PAPIf_flips( real_time, proc_time, flpins, mflops, retval )
      if (retval .NE. PAPI_OK) call errquit
     $ (' perfm_start: could not init papi flpins ',retval,0)
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
     $ perfm_error(retval,
     $     'perfm_start: could not papi event2name ')
      call PAPIf_add_event( es, nativees, retval )
      if (retval .NE. PAPI_OK) call 
     $ perfm_error(retval,
     $     'perfm_start: could not papi addevent ')
      enddo
c
      call PAPIf_get_real_usec(uso)
      call PAPIf_start(es, retval)
      if (retval .NE. PAPI_OK) call errquit
     $ (' perfm_start: could not papi start ',retval,0)
#endif
#endif
      return
      end
#ifdef USE_PAPI
c
c
c
      subroutine perfm_error(retval,where)
      implicit none
#include "perfm.fh"
      character*(*) where
         call PAPIF_perror(retval,papi_errstr,iret)
         write(*,*) 'PAPI error: ',
     ,        papi_errstr(1:inp_strlen(papi_errstr))
         call errquit
     $ (where,retval,0)
         return
         end
#endif
