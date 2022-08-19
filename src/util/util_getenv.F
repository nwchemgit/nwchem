      subroutine util_getenv(env_name, env_value)
* $Id$
      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      character*(*) env_name  ! [input] the name of the environment variable to check for
      character*(*) env_value ! [output] the value of the named environement variable
#if defined(CRAY)
      integer*4 lname, lvalue, lerror
#endif
*
#if defined(CRAY)
      lname  = 0
      lvalue = 0
      lerror = 0
      call pxfgetenv(env_name,lname,env_value,lvalue,lerror)
      if (lerror.eq.0) then  ! 0 ok found env_name
        return
      else if (lerror.eq.22) then  ! 22 env_name is not in environment
        env_value = ' '
        return
      else
        call errquit('util_getenv:pxfgetenv fatal error',
     *        lerror, UNKNOWN_ERR)
      endif
#else
      call getenv(env_name, env_value)
#endif
*
      end

#ifdef USE_OFFLOAD
      integer function offload_span()
      implicit none
      character(64) env_string
      integer span
      integer stat
      span = -1
      call util_getenv('NWC_OFFLOAD_SPAN', env_string)
      if ((len(env_string).ne.0).or.(env_string.ne.' ')) then
        read (env_string,*,iostat=stat) span
      endif
      if (stat.ne.0) span = 1
      if (span.lt.1) span = 1
      offload_span = span
      end
#endif


#ifdef USE_OFFLOAD
c     Thus function is a place holder for a function that
c     will assign offload capabilities to the calling GA rank.
      logical function offload_enabled()
      implicit none
#include "global.fh"
      integer omp_get_num_devices
      external omp_get_num_devices
      character(64) env_string
      integer span
      external offload_span
      integer offload_span
      call util_getenv('NWC_NO_OFFLOAD', env_string)
      if (env_string.eq.'T'.or.env_string.eq.'t'.or.
     c    env_string.eq.'1') then
        offload_enabled=.false.
      else if (omp_get_num_devices().eq.0) then
        offload_enabled=.false.
      else
#ifdef DEBUG         
         write(6,*) '$$ omp_get_num_devices=',
     C        omp_get_num_devices()
#endif         
        span = offload_span()
        if (mod(ga_nodeid(), span).eq.0) then
          offload_enabled=.true.
        else
          offload_enabled=.false.
       endif
#ifdef DEBUG       
       write(6,*) '$$ offl_span ',span,
     C      ' offload_enabled=',offload_enabled
#endif       
      endif
      end
#endif


c     This function is a place holder for a function that
c     will assign actual offload devices to the GA rank w/
c     offloading enabled.
#ifdef USE_OFFLOAD
      integer function offload_device()
      implicit none
#include "global.fh"
#if USE_OPENMP
      external offload_span
      integer offload_span
      integer omp_get_num_devices
      external omp_get_num_devices
      integer omp_get_default_device
      external omp_get_default_device
      external omp_set_default_device
#endif
      character(64) env_string
      integer ndev
      integer span

      span = offload_span()
#if USE_OPENMP
      ndev = omp_get_num_devices()
      write(6,*) ' omp_get_num_devices() ',ndev
      offload_device = mod(ga_nodeid() / span, ndev)
      call omp_set_default_device(offload_device)
#else
      offload_device = -1
#endif
      end
#endif
