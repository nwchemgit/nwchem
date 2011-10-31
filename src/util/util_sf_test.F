c
c FNAME - filename for test program
c
#define BASE_NAME 'shared.file'
#ifdef  HPIODIR
#  define FNAME   HPIODIR//BASE_NAME
#else
#  define FNAME   BASE_NAME
#endif
#define SYNCH ga_sync

      subroutine util_sf_test
c
      implicit none
#include "sf.fh"
c
      integer size
      integer dimnsn
      integer maxid
      parameter(maxid=5)           ! max number of outstanding I/O requests
      parameter(dimnsn=8192)
      real*8 buffer(dimnsn,maxid)  ! need buffering for all maxid requests
      double precision tt0,ttw, ttr, tcgtime
c     integer stack, heap
      integer idlist(maxid), curid
      character*80 errmsg
      logical status
c
      integer nodeid, nnodes
      integer me,nproc,handle,i,start,end,j, rc, chunk
c
c     call pbeginf()
c     call ga_initialize()
c
      curid = 0
      me = nodeid()
      nproc = nnodes()
      size  = maxid*dimnsn*nproc
c
      status = .true.
      rc = sf_create(FNAME,
     $     dble(16*size),dble(8*size),dble(8*dimnsn),handle)
      if (rc.ne.0) status = .false.
      if (me.eq.0) then
         if (status) then
         write(6,'(A)') ' sf_create ........................ OK'
         else
         write(6,'(A)') ' sf_create ........................ Failed'
         endif
         call ffflush(6)
      endif
c
      call SYNCH
      chunk = (size+nproc-1)/nproc
      start = me*chunk+1
      end = min((start+chunk-1),size)
      tt0 = tcgtime()
      status = .true.
c
c     print *,'me=',me,'writing:', start, end
c     everybody writes chunk of data
      if(start.le.end) then
         do i = start, end,dimnsn
            do j = 1, min(dimnsn,(end-i+1))
               buffer(j,curid+1) = dble(i+j-1)
            enddo

            if(curid .eq. maxid)then
               rc=sf_waitall(idlist,maxid)
               curid = 0
            endif
            curid = curid+1
            rc = sf_write(handle,  8*dble(i-1), 
     &           8*dble(min(dimnsn,(end-i+1))), buffer(1,curid), 
     &           idlist(curid))
            if (rc.ne.0) status = .false.

         enddo
      endif
c
      rc=sf_waitall(idlist,curid)
      if(rc.ne.0) status = .false.
      curid = 0
      if (me.eq.0) then
         if (status) then
         write(6,'(A)') ' sf_write/sf_waitall .............. OK'
         else
         write(6,'(A)') ' sf_write/sf_waitall .............. Failed'
         endif
         call ffflush(6)
      endif

      ttw = tcgtime() -tt0
c
      call ga_dgop(88,ttw,1,'max')
      call SYNCH
c
c
c     everybody reads different chunk of data
      start = (nproc-me-1)*chunk+1
      end = min((start+chunk-1),size)
c     print *,'me=',me,'reading:', start, end
      tt0 = tcgtime()
      status = .true.
      do i = start,end,dimnsn

c           read and test data chunk by chunk
            rc = sf_read(handle, 8*dble(i-1),
     .           8*dble(min(dimnsn,(end-i+1))), buffer, idlist(1))
            if (rc.ne.0)then
               status = .false.
            endif
            rc=sf_wait(idlist(1))
            if (rc.ne.0) status = .false.
c
            do j = 1,min(dimnsn,(end-i+1))
               if(buffer(j,1).ne.dble(i+j-1)) then
                  status = .false.
               endif
            enddo
      enddo
      ttr = tcgtime() -tt0
c
      call ga_dgop(88,ttr,1,'max')
      call SYNCH
      if (me.eq.0) then
         if (status) then
         write(6,'(A)') ' sf_read/sf_wait .................. OK'
         else
         write(6,'(A)') ' sf_read/sf_wait .................. Failed'
         endif
         call ffflush(6)
      endif
c
      status = .true.
      rc = sf_destroy(handle)
      if(rc.ne.0) status = .false.
      if (me.eq.0) then
         if (status) then
         write(6,'(A)') ' sf_destroy ....................... OK'
         else
         write(6,'(A)') ' sf_destroy ....................... Failed'
         endif
         call ffflush(6)
      endif
      if(me.eq.0)then
c       write(*,*)' '
c       write(*,*)'test passed ', 8*maxid*dimnsn,' bytes'
c       write(*,*) 8.0e-6*maxid*dimnsn/ttw,' MB/s write rate'
c       write(*,*) 8.0e-6*maxid*dimnsn/ttr,' MB/s read rate'
c       write(*,*)' '
        write(6,'(A,F7.1,A)') ' sf_write rate .................... ',
     1  8.0e-6*maxid*dimnsn/ttw,' MB/s'
        write(6,'(A,F7.1,A)') ' sf_read rate ..................... ',
     1  8.0e-6*maxid*dimnsn/ttr,' MB/s'
        call ffflush(6)
      endif
c
      call SYNCH
c     call ga_terminate
c     call pend()
      end
c $Id$
