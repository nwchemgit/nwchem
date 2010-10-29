      program dumpdb
C$Id$
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      character*128 name
      integer rtdb
      integer crap
c
      name = 'h2o.db'
c
      call pbeginf
      if (.not. ma_init(MT_DBL, -1,-1))
     $     call errquit(' ma_init failed ', 0, MA_ERR)
      if (.not.rtdb_open( name, 'empty', rtdb))
     $     call errquit(' open 1 failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'Test ', MT_INT, 1, 1))
     $     call errquit(' put failed ', 1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'Test ', MT_INT, 1, crap))
     $     call errquit(' get failed ', 1, RTDB_ERR)
      write(6,*) ' crap = ', crap
      call ffflush(6)
      crap = 0
      if (.not. rtdb_print( rtdb, .true.))
     $     call errquit(' print failed ', 1, RTDB_ERR)
      if (.not. rtdb_close(rtdb, 'keep'))
     $     call errquit(' close failed ', 1, RTDB_ERR)
c
      if (.not. rtdb_open( name, 'old', rtdb))
     $     call errquit(' open 2 failed', 1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'Test ', MT_INT, 1, crap))
     $     call errquit(' get 2 failed ', 1, RTDB_ERR)
      write(6,*) ' crap = ', crap
      call ffflush(6)
      if (.not. rtdb_put(rtdb, 'Crap ', MT_INT, 1, crap))
     $     call errquit(' put 2 failed ', 1, RTDB_ERR)
      if (.not. rtdb_print( rtdb, .true.))
     $     call errquit(' print failed ', 1, RTDB_ERR)
      if (.not. rtdb_close(rtdb, 'keep'))
     $     call errquit(' close failed ', 1, RTDB_ERR)
      call pend
c
      end
      subroutine ffflush(i)
      call flush(i)
      end
