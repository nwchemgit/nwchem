      integer function nxtasksub(nproc,ichunk,p_handle)
      implicit none
#include "errquit.fh"
#include "global.fh"      
      logical odlb
      common /dlbinfo/ odlb
c#ifdef GANXTVAL
      logical f1stcall
      save f1stcall
      data f1stcall /.true./
      external util_tcesublock
      integer util_tcesublock
c#endif
c
      integer nproc, ichunk     ! [input]
      integer p_handle
c
      integer icount, nleft, nxtval, junk
      external nxtask
      save icount, nleft
      data nleft, icount /0, 0/
c#ifdef GANXTVAL
      if(f1stcall) then
         junk=util_tcesublock(0,p_handle)
         f1stcall=.false.
      endif
c#endif
c     
c     wrapper round nxtval() to increase granularity
c     and thus reduce no. of requests to shared counter
c
c     !! ichunk MUST not be changed during a parallel loop
c     !! AND all processors MUST agree on the value
c     
*      if(nproc.gt.0) then

      odlb = .true.

      if (odlb) then
      if(nproc.eq.ga_pgroup_nnodes(p_handle)) then
         if(nleft.eq.0) then
            icount = util_tcesublock(nproc,p_handle) * ichunk
            nleft = ichunk
         endif
         nxtasksub = icount
         icount = icount + 1
         nleft = nleft -1
*      else if (nproc.lt.0) then
      else if (nproc.eq.-ga_pgroup_nnodes(p_handle)) then
         nleft = 0
         nxtasksub = -1
         junk = util_tcesublock(nproc,p_handle) 
         f1stcall = .true.
      else 
         call errquit('nxtask: invalid nproc ', nproc, UNKNOWN_ERR)
      endif
c     
c     following does dumb static load balancing
c     
      else

         call errquit('nxtask: STATIC DIST !!! ', nproc, UNKNOWN_ERR)

      if(nproc.gt.0) then
         if (nleft .eq. 0) then
            icount = ga_nodeid()
            nleft = 1
         endif
         nxtasksub = icount
         icount = icount + ga_nnodes()
      else
         nleft = 0
         nxtasksub = 0
      endif

      endif
c     
c      if(nxtasksub.eq.-1) call ga_pgroup_sync(p_handle)

c      if(p_handle.eq.1) then
c        write(6,*) ga_nodeid(), ' got task ', nxtasksub
c      endif
c      call forflush(6)
c     
      end


c $Id$
