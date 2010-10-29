c********************************************
      program test
*
* $Id$
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "pstat.fh"
#include "util.fh"
#include "inp.fh"

#include "myprocs.fh"
#include "parms.fh"
      integer idummy(64),p
      integer id1,id2,id3, stack,heap,global
c---------------------

      if (proctype.EQ.mpi) then
         call mpiInit
      else if (proctype.EQ.tcg) then
        call pbeginf()
      else 
        print*,'FATAL: invalid proctype = ',proctype
        stop
      endif

      me=nodeid()
      nproc=nnodes()

c	print*,me,' lives with proctype = ',proctype

      heap = nd1*nd2*nd3*16
      stack = 1024*128
      if (.NOT. MA_init(mt_dbl, stack, heap)) then
         print*,' MA_init failed'
      endif

      if (me.EQ.0) print*,' *** FFT Type: ',fft_type,
     .             ' with ',procs_type,
     .             ' and nproc = ',nproc,' ***'

c     do p = 0, nproc-1
c        if (me.EQ.p) print*,me,' test fft3 lives'
c        call mysync()
c     enddo
c     call print2D(idummy,1,4,'after mpi init')

        id3 = MAX(md3,nproc*2)

        if (id3.GT.nd3) then
          if (me.EQ.0) then
            print*,' too many processors for matrix size',nproc,id3
          endif
c         call mpiStop
          call pend
        endif

        do while (id3.LE.nd3)
          call fft3(id3,id3,id3)
         call mysync
          id3 = id3*2 
        enddo

c     call mpiStop
      call pend

      end

c******************************************************************
      subroutine fft3(id1,id2,id3)
      implicit none
      INTEGER id1,id2,id3
#include "myprocs.fh"
#include "parms.fh"
      double precision x1(nd1*nd2*nd3*2)
      double precision x2(nd1*nd2*nd3*2)
      double precision x3(nd1*nd2*nd3*2)
c --
      integer blksz, islab, isign
      integer j,k,owner,mapindx
      double precision t0(2),t1(2),myclock, inorm, infnorm
      integer map1(nd3), map2(nd3), doMap
      integer d1, d2, d3
      double precision fftime(2), util_wallsec
      integer idum(1)
      double precision rdum(1)
c------------------------------
c     print*,me,' enters fft3 with: ',id1,id2,id3
 
      call mysync

c     *initialize scratch arrays y1 and y2
      blksz = id3/nproc

      owner = 0
      k = 1
      do d3 = 1, id3
        if (d3.GT.blksz*(owner+1)) then
          owner = owner+1
        endif
        if (me.EQ.owner) then
          do d2 = 1, id2
            do d1 = 1, id1
              x1(k) = dble(d1-1)
              x1(k+1) = 0.D0
              k = k+2
            enddo
          enddo
        endif
      enddo

#ifdef FFTTIME
      call mysync
      fftime(1) = util_wallsec()
#endif


c ** forward fft
      doMap = 1
      isign = 1
      call nwfft3d(id1,id2,id3,isign,x1,x2,map1,map2,doMap,
     .             0,idum,0,rdum,istatus)


c ** reverse fft
      doMap = 0
      isign = -1
      call nwfft3d(id1,id2,id3,isign,x2,x3,map1,map2,doMap,
     .             0,idum,0,rdum,istatus)

#ifdef FFTTIME
      fftime(2) = util_wallsec()
      if (me.EQ.0) write(6,11) id1,id2,id3,fftime(2)-fftime(1)
11    format(3(I3,1X),'fft wallclock time (sec): ',F10.5)
#endif


c     check infinity norm
      infnorm = inorm(x1,x3,id1*id2*blksz*2)
      if (infnorm.GT.0.000001) then
         write(6,*) me,id1,id2,id3,' infinity norm=',infnorm
      endif


      return
      end

c******************************************************************


      subroutine print2D(x,id1,id2,blksz,s)
      implicit none
      double precision x(*)
      character*(*) s
      integer id1,id2, d1,d2,d3, i, p, ioff, blksz
#include "myprocs.fh"

      call mysync
      if (me.EQ.0) then
         write(6,*) 
         write(6,*) s
      endif
      do p = 0, nproc-1
        if (p.EQ.me) then
          do d3 = 1, blksz
            do d2 = 1, id2
              ioff = id1*id2*(d3-1)*2 + id1*(d2-1)*2 
c              ioff = id1*(d2-1) + d1
              write(6,10) me,d2,d3,(x(i),i=1+ioff,id1*2+ioff) 
            enddo
          enddo
        endif
        call mysync
      enddo
10    format(I2,'(y,z)',2(1X,I2),' : ',16F4.1)

      return
      end

c************************************************************************
      subroutine xerbla
      print*,' xerbla stub'
      return
      end
c************************************************************************

      subroutine mysync
#include "msgtypesf.h"
      call synch(msgint)
      return
      end

c************************************************************************
