C$Id: nxtask.F,v 1.2 2004-10-21 01:36:52 edo Exp $
      subroutine nxtask_init(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
C  This is needed so that subgroups work, GA has MPI_COMM_WORLD in it
C  But, it does not work on all platforms
C#define GANXTVAL 1
      integer rtdb
c
c     load dlb info
c     odlb = if false use static load balance
c
      logical odlb
      common /dlbinfo/ odlb

      if (.not.rtdb_get(rtdb,'loadbalancing:dynamic',MT_LOG,1,odlb))
     $ odlb = .true.
      end

      integer function nxtask(nproc, ichunk)
      implicit none
#include "errquit.fh"
#include "global.fh"      
      logical odlb
      common /dlbinfo/ odlb
#ifdef GANXTVAL
      logical f1stcall
      save f1stcall
      data f1stcall /.true./
      external util_gnxtval
      integer util_gnxtval
#endif
c
      integer nproc, ichunk     ! [input]
c
      integer icount, nleft, nxtval, junk
      external nxtval
      save icount, nleft
      data nleft, icount /0, 0/
#ifdef GANXTVAL
      if(f1stcall) then
         junk=util_gnxtval(0)
         f1stcall=.false.
      endif
#endif
c     
c     wrapper round nxtval() to increase granularity
c     and thus reduce no. of requests to shared counter
c
c     !! ichunk MUST not be changed during a parallel loop
c     !! AND all processors MUST agree on the value
c     
*      if(nproc.gt.0) then

      if (odlb) then
      if(nproc.eq.ga_nnodes()) then
         if(nleft.eq.0) then
#ifdef GANXTVAL
            icount = util_gnxtval(nproc) * ichunk
#else
            icount = nxtval(nproc) * ichunk
#endif
            nleft = ichunk
         endif
         nxtask = icount
         icount = icount + 1
         nleft = nleft -1
*      else if (nproc.lt.0) then
      else if (nproc.eq.-ga_nnodes()) then
         nleft = 0
         nxtask = -1
#ifdef GANXTVAL
         junk = util_gnxtval(nproc)
         f1stcall = .true.
#else
         junk = nxtval(nproc)
#endif
      else 
         call errquit('nxtask: invalid nproc ', nproc, UNKNOWN_ERR)
      endif
c     
c     following does dumb static load balancing
c     
      else

      if(nproc.gt.0) then
         if (nleft .eq. 0) then
            icount = ga_nodeid()
            nleft = 1
         endif
         nxtask = icount
         icount = icount + ga_nnodes()
      else
         nleft = 0
         nxtask = 0
      endif

      endif
c     
*      write(6,*) ga_nodeid(), ' got task ', nxtask
*      call forflush(6)
c     
      end
      logical function util_statldb(lstatic,rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      logical lstatic
c
c     set static load balancing if lstatic eq t
c     set dynamic load balancing if lstatic eq f
c     return value = true if old oldb is static
c
      logical odlb
      common /dlbinfo/ odlb
      integer rtdb
      
      util_statldb=.not.odlb
      if(lstatic) then
         odlb=.false.
      else
         odlb=.true.
      endif
      if (.not.rtdb_put(rtdb,'loadbalancing:dynamic',MT_LOG,1,odlb))
     &     call errquit('util_statldb: put failed',0,RTDB_ERR)

      return
      end
