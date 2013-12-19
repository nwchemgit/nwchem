      subroutine util_mygabcast(g_a,m,n,a,ld)
      implicit none
#include "global.fh"
#include <mpif.h>
      integer g_a
      integer m,n,ld
      double precision a(*)
c
      integer i
      integer*4 ierr,len,ierr2
      integer*8 bigint,len8,nsteps,istart,end8
cold      parameter (bigint=1.999999*2.**30)
      parameter (bigint=2.**29/8)
      character*(MPI_MAX_ERROR_STRING) errstr
c
      if(ga_nodeid().eq.0)
     G     call ga_get(g_a,1,m,1,n,a,ld)
      len8=m*n
      if(len8.lt.bigint) then
         len=m*n
         call MPI_Bcast(a, len, MPI_DOUBLE_PRECISION, 0, 
     G        MPI_COMM_WORLD, ierr)
         else
            nsteps=ceiling(dble(len8)/dble(bigint))
c            write(0,*) ' nsteps ',nsteps,
c     F           dble(len8)/dble(bigint)
            istart=1
            do i=1,nsteps
                      len=bigint
               end8=istart+bigint-1
               if(end8.gt.len8) then
                  end8=len8-istart+1
                  len=end8
               endif
#ifdef DBG
               if(ga_nodeid().eq.0) then
                  write(6,'(A,I10,A,I10,A,I10,A,I10)') 
     P                 ' istart ',istart,'len ',len,
     I                 ' is+len-1 ',istart+len-1, 
     L                 ' len8 ',len8
               endif
#endif
               call MPI_Bcast(a(istart), len, MPI_DOUBLE_PRECISION, 0, 
     G              MPI_COMM_WORLD, ierr)
#ifdef DBG
               call ga_sync()
#endif
               istart=istart+len
            enddo
        endif
      call MPI_ERROR_STRING(ierr, errstr, len, ierr2)
      if(ierr.ne.0) then
         if(ga_nodeid().eq.0.)
     W        write(0,*) ' mpibcast ',ierr,errstr(1:len)
         call errquit(' utilmygabcast failed ',0,0)
      endif
      return
      end
c $Id$
