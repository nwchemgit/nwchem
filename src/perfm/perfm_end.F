c
c     $Id: perfm_end.F,v 1.1 2004-10-11 19:39:05 edo Exp $
c
      subroutine perfm_end()
      implicit none
#include "perfm.fh"
      double precision conv2m
      parameter(conv2m=1d-6)
#ifdef USE_PAPI
      integer*8 flpsum,flpmax,flpmin
      integer retvals
#ifdef PAPI_FPINS
      call ga_sync()
      call PAPIf_flips( real_time, proc_time, flpins, mflops, retval)
      call ga_sync()
      retvals=retval
      if (retval .NE. PAPI_OK) call perfm_error(retval,
     $     ' perfm_end: could not get papi flpins ')
c
c     sum,max and min dgops
c     
      flpsum=flpins
      flpmax=flpins
      flpmin=flpins
      call ga_dgop(msg_perfm_max,flpmax,1, 'max')
      call ga_dgop(msg_perfm_min,flpmin,1, 'min')
      call ga_dgop(msg_perfm_sum,flpsum,1, '+')
c
      if(ga_nodeid().eq.0) then
         write(luout,*) 'Real_time: ', real_time
         write(luout,*) ' Proc_time: ', proc_time
         write(luout,*) ' Total flpins: ', flpins
         write(luout,*) ' MFLOPS: ', mflops
         write(luout,12) ' Aggregate MFLOP : ', 
     ,        int(flpsum*conv2m)
         write(luout,12) ' Aggregate MFLOPS (Real_time): ', 
     ,        int(flpsum/real_time*conv2m)
         write(luout,12) ' Aggregate MFLOPS (Proc_time): ', 
     ,        int(flpsum/proc_time*conv2m)
         write(luout,12) ' Average   MFLOPS: ', 
     ,        int(conv2m*flpsum/real_time/ga_nnodes())
         write(luout,12) ' Max       MFLOPS: ', 
     ,        int(flpmax/real_time*conv2m)
         write(luout,12) ' Min       MFLOPS: ', 
     ,        int(flpmin/real_time*conv2m)
 12      format(A32,I10) 
      endif
#else
      call PAPIf_stop(es, values(1), retval)
      if (retval .NE. PAPI_OK) call
     $ perfm_error(retval,
     $     'perfm_start: could not stop papi')
      call PAPIf_get_real_usec(usn)
      usn = usn - uso
      write(luout,*) 'Real_time: ', usn, 'secs'
      ll=inp_strlen(eventname(1))
      write(luout,*) eventname(1)(1:ll),'      : ',values(1)
      write(luout ,*) eventname(1)(1:ll),'/sec : ',values(1)/usn
#endif
#endif
      return
      end
