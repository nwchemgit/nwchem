C$Id$
      subroutine nxtask_init(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
C  This is needed so that subgroups work, GA has MPI_COMM_WORLD in it
C  But, it does not work on all platforms
#define GANXTVAL 1
#define SKIPFIRST_ 1 
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
      integer nproc  ! [in]
      integer ichunk ! [in]
c
      external nxtask0
      integer nxtask0
c
      nxtask=nxtask0(nproc, ichunk,0)
      return
      end
      integer function nxtask0(nproc, ichunk, skipfirst)
      implicit none
#include "errquit.fh"
#include "global.fh"      
      logical odlb
      common /dlbinfo/ odlb
      logical f1stcall
      save f1stcall
      data f1stcall /.true./
      external util_shcount
      integer util_shcount
c
      integer nproc, ichunk     ! [input]
      integer skipfirst ! [in] .eq.0 then noskipfirst
c
      integer icount, nleft, nxtval, junk
      external nxtval
      save icount, nleft
      data nleft, icount /0, 0/
c     we do not like negative values for skipfirst
      if(skipfirst.lt.0) call errquit(
     A     ' nxtask0: skipfirst has negative value',skipfirst,0)
      if(skipfirst.gt.1.and.ichunk.gt.1) call errquit(
     A     ' nxtask0: case with skip1st.gt.1.and.ichunk.gt.1 not ready',
     S     skipfirst,0)
      if (odlb) then
        if(f1stcall) then
c     skip read inc for first call                                              
           junk=util_shcount(0)
           if(skipfirst.eq.0) f1stcall=.false.
        endif
      endif
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
            nleft = ichunk
c     skip read inc for first call   
            if(f1stcall) then
               nleft=nleft*skipfirst
               f1stcall=.false.
               icount=(ga_nodeid()) * ichunk * skipfirst
            else
               icount = (util_shcount(nproc) +
     C              skipfirst*ga_nnodes()) * ichunk
            endif
         endif
         nxtask0 = icount
         icount = icount + 1
         nleft = nleft -1
*      else if (nproc.lt.0) then
      else if (nproc.eq.-ga_nnodes()) then
         nleft = 0
         nxtask0 = -1
         junk = util_shcount(nproc)
         f1stcall = .true.
      else 
         call errquit('nxtask0: invalid nproc ', nproc, UNKNOWN_ERR)
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
         nxtask0 = icount
         icount = icount + ga_nnodes()
      else
         nleft = 0
         nxtask0 = 0
      endif

      endif
c     
*      write(6,*) ga_nodeid(), ' got task ', nxtask0
*      call util_flush(6)
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
      integer function util_shcount(nproc)
      implicit none
      integer nproc
#ifdef GANXTVAL
      integer util_gnxtval
      external util_gnxtval
      util_shcount=util_gnxtval(nproc)
#else
      integer nxtval
      external nxtval
c     no init with tcg nxtval??
      if(nproc.eq.0) return
      util_shcount = nxtval(nproc)
#endif
       return
       end
