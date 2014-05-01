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
      integer function util_nwcnumthreads(rtdb)
      implicit none
#include "rtdb.fh"      
      integer rtdb
      character*20 env_value 
      integer nnn,num_devices
c
#if   defined(USE_OPENMP) && defined(USE_OFFLOAD)
#define RANKS_PER_DEVICE 2
      external util_mic_get_num_devices
      integer*4 util_mic_get_num_devices
      num_devices=util_mic_get_num_devices()
#elif defined(TCE_CUDA)
#define RANKS_PER_DEVICE 1
      if(.not.rtdb_get(rtdb,'tce:cuda',mt_int,1,num_devices))
     C     num_devices=0
#else
!#warning ndev=0
      num_devices=0
      call errquit(' util_nwcnumthreads: used for coprocessors',0,0)
#endif
#ifndef RANKS_PER_DEVICE 
#define RANKS_PER_DEVICE 1
#endif
      call util_getppn(nnn)
      util_nwcnumthreads=nnn/RANKS_PER_DEVICE/
     D     num_devices
      return
      end
      logical function offload_master()
      implicit none
#include "global.fh"
      external util_nwcnumthreads
      integer util_nwcnumthreads
c                                                                               
      integer nn
c                                                                               
      offload_master=.false.
      nn=util_nwcnumthreads()
      if(nn.eq.0) return
      offload_master=mod(ga_nodeid(),util_nwcnumthreads()).eq.0
      return
      end
      integer function util_micdev()
#if defined(USE_OPENMP) && defined(USE_OFFLOAD)
      use mic_lib
      implicit none
#include "global.fh"
#include "errquit.fh"
      external offload_master
      logical offload_master
c                                                                               
      integer numdevs
      integer*4 util_mic_get_device
      external util_mic_get_device
c                                                                               
      if(.not.offload_master()) call errquit(
     A     ' need to be offload master',0,0)
      util_micdev=util_mic_get_device()
#else
      util_micdev=-1
      call errquit(' micdev: wrong installation',0,0)
#endif
      return
      end
