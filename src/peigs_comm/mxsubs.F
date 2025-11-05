*
* $Id$
*
*======================================================================
*
* DISCLAIMER
*
* This material was prepared as an account of work sponsored by an
* agency of the United States Government.  Neither the United States
* Government nor the United States Department of Energy, nor Battelle,
* nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
* ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
* COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
* SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
* INFRINGE PRIVATELY OWNED RIGHTS.
*
* ACKNOWLEDGMENT
*
* This software and its documentation were produced with Government
* support under Contract Number DE-AC06-76RLO-1830 awarded by the United
* States Department of Energy.  The Government retains a paid-up
* non-exclusive, irrevocable worldwide license to reproduce, prepare
* derivative works, perform publicly and display publicly by or for the
* Government, including the right to distribute to other Government
* contractors.
*
*======================================================================
*
*  -- PEIGS  routine (version 2.1) --
*     Pacific Northwest Laboratory
*     July 28, 1995
*
*======================================================================
c	mxsubs.cpp,v 1.7 1995/04/05 21:53:49 d3c423 Exp
c
c defining UNSAFE causes mxread to ignore the count and source node,
c taking it on faith that what was received is what was expected.
#define UNSAFE

C     Portable interface to communications routines, following EXPRESS
C     semantics.  These are (supposed to be) ALL the parallel 
C     control/comm routines that are used by the programs in this
C     directory.  Note that subroutines ksingl, kasync, and kmulti
C     are related to node I/O and are inherently non-portable.
C     When not running in an EXPRESS environment, those routines
C     expand to nothing.
C
C     Rik Littlefield, 1/16/91.
C
#ifdef ALOG
      include 'alogevents.h'
#endif
      subroutine mxtime( job, xtime )
      double precision xtime
      integer          job
c
c     mxtime -- job = 0, then initialize total mxread + mxwrit time
c                        to 0.0 for this processor
c
c                   <>0, returns total (mxread time) + (mxwrit time)
c                        since last call to mxtime with job = 0
c                        for this processor

      double precision t1, t2, time
      common /timemx/  t1, t2, time

      double precision mxclock
      external         mxclock

      if( job .eq. 0 ) then

c       make sure mxclock called once before start timing
        t1 = mxclock()

        time = 0.d0
      endif

      xtime = time

      end
      subroutine mxinit
C     mxinit -- initialization.  Note that this implementation assumes
c               that EXPRESS fills in the 'symbolic constants' in
c               /xpress/.  These constants just need to be recognizably
c               different from valid node numbers, numbers of nodes,
c               and processor counts.  (At present, the only one
c            used by the diagonalizer is 'ialnod', meaning "all nodes".)
c

      common /xpress/ nocare, norder, nonode, ihost, ialnod, ialprc
#include "mxsubs.h"
      integer ienv(4)
#ifdef ALOG
#include "alogf.h"
#endif
      integer ifirst
      integer nocare, norder, ihost, ialnod, ialprc, nonode
	double precision mxclock
#ifdef KSR
	real tarray(2)
	real t1
#endif
#ifdef MPI
      LOGICAL FLAG
      INTEGER IERROR
#endif
      save ifirst
      data ifirst /1/

c        Do not reinitialize the underlying comm package.

      if (ifirst.eq.0) return
      ifirst = 0

#ifdef EXPRESS
      call kxinit
#else
      nocare = -1
      norder = -2
      nonode = -3
      ihost = -4
      ialnod = -5
      ialprc = -6
#endif
#ifdef MPI
      CALL MPI_INITIALIZED( FLAG, IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXINIT: ERROR IN MPI_INITIALIZED'
        STOP
      ENDIF

      IF( .NOT. FLAG) THEN

        CALL MPI_INIT( IERROR )

        IF( IERROR .NE. 0 ) THEN
          WRITE(*,*) ' ROUTINE MXINIT: ERROR IN MPI_INIT'
          STOP
        ENDIF

      ENDIF
#endif

#ifdef TCGMSG
c
      t1 = mxclock()
      call mxpara (ienv)
      me = ienv(1)
      nproc = ienv(2)
#endif
#ifdef KSR
      t1 = mxclock()
      call mxpara (ienv)
      me = ienv(1)
      nproc = ienv(2)
#endif
#ifdef iPSC_NATIVE
#endif
#ifdef PICL
      call tracenode (2000000,1,1)
      call tracelevel (4,4,0)
      call open0 (nproc,me,host)
#endif
#ifdef UNIPROC
#endif
      call mxpara (ienv)
      me = ienv(1)
      nproc = ienv(2)
#ifdef ALOG
      ipid = me
      call alogfs (ipid,0)
      call alogfg (ipid, AE_OPEN, 0, sdata, 0)
#endif
      return
      end

      subroutine mxpara (ienv)
c     mxpara -- return environment.
c                 ienv(1) = me (this processor number)
c                 ienv(2) = nprocs (total available node processors)
c                 ienv(3..4) reserved
c

      integer ienv(4)
      integer nodeid, nnodes, mynodes, numnodes
#include "mxsubs.h"
#ifdef MPI
	include 'mpif.h'
      INTEGER IERROR
#endif
#ifdef iPSC_NATIVE
	external  mynodes, numnodes
#endif
#ifdef TCGMSG
 	external nodeid, nnodes
#endif
#ifdef UNIPROC
      ienv(1) = 0
      ienv(2) = 1
#endif
#ifdef EXPRESS
      call kxpara (ienv)
#endif
#ifdef MPI
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, IENV(1), IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXPARA: ERROR IN MPI_COMM_RANK'
        STOP
      ENDIF

      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, IENV(2), IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXPARA: ERROR IN MPI_COMM_SIZE'
        STOP
      ENDIF
#endif
#ifdef TCGMSG
      ienv(1) = nodeid()
      ienv(2) = nnodes()
#endif
#ifdef iPSC_NATIVE
      ienv(1) = mynode ()
      ienv(2) = numnodes ()
#endif
#ifdef PICL
      ienv(1) = me
      ienv(2) = nproc
#endif
      me     = ienv(1)
      nproc  = ienv(2)
      host   = 0
c     
      return
      end

      integer function mxtick ()
c     mxtick -- return cpu time in arbitrary increments.
c
#include "mxsubs.h"
      
#ifdef MPI
      double precision mxclock
      external mxclock
      mxtick = mxclock()*1000
#endif

#ifdef TPUTER_NODE
      mxtick = kxtime ()
#endif
#ifdef EXPRESS
      mxtick = kxtick ()
#endif
#ifdef TCGMSG
      mxtick = mtime()
#endif
#ifdef iPSC_NATIVE
      mxtick = mclock ()
#endif
#ifdef PICL
      double precision clock0
c     .      write (*,*) me,' calling clock0'
      mxtick = 1000 * clock0 ()
c     .      write (*,*) me,' return from clock0'
c..   mxtick = mclock ()
#endif
#ifdef UNIPROC
      real et(2)
c     call etime (et)
      mxtick = et(1)*1000
#endif
      return
      end
      
	integer function mxread (data,isize,iproc,itype)
c     mxread -- read data from another processor.
c               As used by the Jacobi diagonalizer, iproc and itype
c               are always specified exactly.  This routine does not
c               need to handle promiscuous types, iproc=ialnod, or
c               other such complexities.
c
	integer data(*), isize, iproc, itype
	integer icopy
#ifdef TIMING
      double precision t1, t2, time
      common /timemx/  t1, t2, time
      double precision mxclock
      external         mxclock
#endif
#ifdef TCGMSG
	external rcv
#endif
#ifdef ALOG
#include "alogf.h"
#endif
#ifdef MPI
	include 'mpif.h'
      INTEGER IERROR, STATUS(MPI_STATUS_SIZE)
#endif
#ifdef TIMING
      t1 = mxclock()
#endif
#ifdef ALOG
      call alogfg (ipid, AE_RCV_ENTER, itype, sdata, 0)
#endif
#ifdef UNIPROC
c      mxread = -999
c      write (*,*) '!!! mxread called in uniprocessor mode'
c
c      stop 999
#endif
#ifdef EXPRESS
      mxread = kxread (data,isize,iproc,itype)
#endif
#ifdef TCGMSG
c     icopy = mxmynd() + 100*iproc
c     icopy = mod(icopy, 32767) + 1
c     
	icopy = mod(itype, 500) + 25001
c     
c     write(*,*) ' me ', mxmynd(), ' rcv ', iproc, ' isize ',
c     $     isize, ' itype ', itype
c     call ffflush(6)
c     
      call rcv(icopy, data, isize, len, iproc, ifrom, 1)
      mxread = len
#endif
#ifdef iPSC_NATIVE
	me = mxmynd()
c
c     -- Note corruption of message type to circumvent iPSC idiocy
c     -- of not allowing crecv to specify node number
c	      write(*,*) 'node ', me, ' rcving type ', 
c     $	me +10000*iproc, 'from ', iproc              		
c
c     call crecv (mod(itype, 512)+10000*iproc,data,isize)
c     
	call crecv (10000*iproc+me,data,isize)
#endif
#ifdef MPI
      icopy = mod(itype, 32767) + 1
      CALL MPI_RECV( DATA, ISIZE, MPI_BYTE, IPROC, ICOPY,
     $               MPI_COMM_WORLD, STATUS, IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXREAD: ERROR IN MPI_RECV'
        STOP
      ENDIF
#endif
#ifdef PICL
c       -- Note corruption of message type to circumvent iPSC idiocy
c       -- of not allowing crecv to specify node number
      call recv0 (data, isize, itype+10000*iproc)
      call recvinfo0 (nbytes,jtype,node)
      mxread = nbytes
      node = infonode ()
      if (node .ne. iproc) then
         write (*,*) '!!!!! mxread: unexpected node -- ',node,' !=',
     &                iproc
         stop 999
      endif
#endif
#ifdef ALOG
      sdata(1) = itype
      sdata(2) = isize
      call alogfg (ipid, AE_RCV_EXIT, iproc, sdata, 2)
#endif
#ifdef TIMING
      t2 = mxclock()
      time = time + (t2 - t1)
#endif
	mxread = isize
      return
      end

      integer function mxwrit (data,isize,iproc,itype)
c     mxwrit -- write data to another processor.
c               As used by the Jacobi diagonalizer, iproc and itype
c               are always specified exactly.  This routine does not
c               need to handle promiscuous types, iproc=ialnod, or
c               other such complexities.
c
      integer data(*), isize, iproc, itype
#ifdef TIMING
      double precision t1, t2, time
      common /timemx/  t1, t2, time
      double precision mxclock
      external         mxclock
#endif
#include "mxsubs.h"
#ifdef ALOG
#include "alogf.h"
#endif
#ifdef MPI
	include 'mpif.h'
      INTEGER IERROR
#endif
#ifdef TIMING
      t1 = mxclock()
#endif
#ifdef UNIPROC
c      mxwrit = -999
c      write (*,*) '!!! mxwrit called in uniprocessor mode'
c      stop 999
#endif
#ifdef EXPRESS
      mxwrit = kxwrit (data,isize,iproc,itype)
#endif
#ifdef TCGMSG
c      icopy = iproc + 100*mxmynd()
c
c     write(*,*) ' me ', mxmynd(), 'snd ', iproc,  ' isize ',
c
c     icopy = iproc + 100*mxmynd()
c     
c     write(*,*) ' me ', mxmynd(), 'snd ', iproc,  ' isize ',
c     $ 	    isize, ' itype ', itype
c     
c     icopy = mod(icopy, 32767) + 1
      icopy = mod(itype, 500) + 25001
c     
c
c to get within range of tcgmsg message types
c
      call snd(icopy, data, isize, iproc, 1)
      mxwrit = isize
#endif
#ifdef MPI
      icopy = mod(itype, 32767) + 1
      CALL MPI_SEND( DATA, ISIZE, MPI_BYTE, IPROC, ICOPY,
     $               MPI_COMM_WORLD, IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXWRIT: ERROR IN MPI_SEND'
        STOP
      ENDIF
#endif
#ifdef iPSC_NATIVE
c        -- Note corruption of type.  See mxread for rationale.
c      write(*,*) 'node ', me, 'sending type', iproc+10000*me, 'to ', iproc
c	
c      call csend (mod(itype, 512)+10000*me , data, isize, iproc, 0)
c
      call csend(iproc+10000*me, data, isize, iproc, 0)
      mxwrit = isize
#endif
#ifdef PICL
c        -- Note corruption of type.  See mxread for rationale.
      call send0 (data, isize, itype+10000*me, iproc)
      mxwrit = isize
#endif
#ifdef ALOG
      sdata(1) = itype
      sdata(2) = isize
      call alogfg (ipid, AE_SEND, iproc, sdata, 2)
#endif
#ifdef TIMING
      t2 = mxclock()
      time = time + (t2 - t1)
#endif
	mxwrit = isize
      return
      end

      integer function mxsync ()
c     mxsync -- synchronize all processors.
c
#include "mxsubs.h"
#ifdef ALOG
#include "alogf.h"
#endif
#ifdef MPI
	include 'mpif.h'
      INTEGER IERROR
#endif
#ifdef UNIPROC
      mxsync = 0
#endif
#ifdef EXPRESS
      mxsync = kxsync ()
#endif
#ifdef MPI
      CALL MPI_BARRIER( MPI_COMM_WORLD, IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXSYNC: ERROR IN MPI_BARRIER'
        STOP
      ENDIF
#endif

#ifdef iPSC_NATIVE
#ifndef DELTA
      call gsync ()
      mxsync = 0
#else
c
c     We should be able to just call gsync here, but on early release
c     Delta (X011), having many processors idle at gsync makes other
c     communications go too slow due a many-to-one implementation that
c     causes long queues on node 0.  So here, we substitute a
c     binary-tree based version that runs in log(P) time and queues
c     no more than two messages at any node.  (Replaces old mxbrod/mxcomb
c     hack.)
c
#include "xpress.com"
      integer parent, child1, child2, itype

      parent = (me-1)/2
      child1 = 2*me + 1
      child2 = child1 + 1
      itype = 9876

      if (child1 .lt. nproc) call mxread (junk,4,child1,itype)
      if (child2 .lt. nproc) call mxread (junk,4,child2,itype)
      if (me .gt. 0) then
        call mxwrit (junk,4,parent,itype)
        call mxread (junk,4,parent,itype+1)
      endif
      if (child1 .lt. nproc) call mxwrit (junk,4,child1,itype+1)
      if (child2 .lt. nproc) call mxwrit (junk,4,child2,itype+1)

      mxsync = 0
#endif
#endif

#ifdef PICL
c       the following does not work correctly and would not have the
c      	proper timing effect anyway.  We fake it.
c.      write (*	,*) me,' in mxsync calling sync0'
c.      call sync0 ()
c.      write (*,*) me,' in mxsync returned from sync0'
      integer nulfnc, mxsync, ijunk
      external nulfnc
      ijunk = 0
      call mxbrod (ijunk,0,4,ialnod,0,7652)
      call mxcomb (ijunk,nulfnc,4,1,ialnod,0,7653)
      mxsync = 0
#endif
#ifdef ALOG
      call alogfg (ipid, AE_SYNC, 0, sdata, 0)
#endif
#ifdef TCGMSG
	ijunk = 1000*mxnprc() + 1
	call synch(ijunk)
        mxsync = 0
#endif
	mxsync = 0
      return
      end


      integer function nulfnc(a,b,isize)
      nulfnc = 0
      return
      end

      integer function mxbcst (data,isize,iproc,itype)
c     mxbcst -- broadcast (from iproc to all others).
c               This implementation uses the slow
c               but reliable approach -- sending processor
c               explicitly sends to each of the receiving ones.
c
      integer data(*)
#include "mxsubs.h"
c.	write (*,*) 'mxbcst: me,iproc,isize,itype=',me,iproc,isize,itype
      if (me.eq.iproc) then
         do 10 ito = 0,nproc-1
            if (ito.ne.me) then
               ireslt = mxwrit (data,isize,ito,itype)
c.               write (*,*) 'send: me,ito,ireslt = ',me,ito,ireslt
            endif
10       continue
      else
         ireslt = mxread (data,isize,iproc,itype)
c.         write (*,*) 'receive: me,iproc,ireslt = ',me,iproc,ireslt
      endif
c.      call mxsync ()
      mxbcst = ireslt
      return
      end

      subroutine mxpend
#ifdef ALOG
#include "alogf.h"
#endif
#ifdef MPI
	include 'mpif.h'
      INTEGER IERROR
#endif
#ifdef PICL
      call close0 ()
#endif
#ifdef MPI
      CALL MPI_FINALIZE( IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXPEND: ERROR IN MPI_FINALIZE'
        STOP
      ENDIF
#endif
#ifdef TCGMSG
        call ga_error(' peigs error: mxpend',0)
#endif
#ifdef ALOG
      call alogfg (ipid, AE_CLOSE, 0, sdata, 0)
      call alogfo ()
#endif
      return
      end

#ifndef TCGMSG
      subroutine llog ()
      return
      end
      subroutine parrerr (n)
	integer n
      write (*,*) 'parrerr: error ',n
      stop
      end
#endif
#ifndef EXPRESS
	subroutine kflush(i)
	integer i
	end
	subroutine kasync(i)
	integer i
	end
	subroutine kmulti(i)
	integer i
	end
	subroutine ksingl(i)
	integer i
	end
#endif
	integer function mxmynd ()
c
c	The obvious "mynode" function, for those tired of using the
c	'env' vector all the time.
c
#include "mxsubs.h"
        mxmynd = me
        return
        end

	integer function mxnprc ()
c
c	The obvious "numprocs" function, for those tired of using the
c	'env' vector all the time.
c
#include "mxsubs.h"
        mxnprc = nproc
        return
        end
c
	double precision function mxclock ()
c		
c	
#ifdef TCGMSG
	double precision tcgtime
	external tcgtime
c	
	mxclock = tcgtime()
#endif
#ifdef MPI
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL         MPI_WTIME
      MXCLOCK = MPI_WTIME()
#endif
#ifdef iPSC_NATIVE
	double precision dclock
	external dclock
	mxclock = dclock()
#endif	
#ifdef UNIPROC
c	integer clock
c	external clock
	mxclock = 0.0d0
#endif
	return
	end

	subroutine mxend
#ifdef MPI
	include 'mpif.h'
      INTEGER IERROR
      CALL MPI_FINALIZE( IERROR )

      IF( IERROR .NE. 0 ) THEN
        WRITE(*,*) ' ROUTINE MXEND: ERROR IN MPI_FINALIZE'
        STOP
      ENDIF
#endif
#ifdef PICL
	call close0(1)
#endif
#ifdef TCGMSG
cedo	call pend()
      call ga_error(' peigs error: mxend',0)
#endif
	return
	end	

