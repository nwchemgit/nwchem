      program test
      implicit none
c
#include 'rtdb.fh'
#include 'mafdecls.fh'
c
      integer rtdb
      integer itest(3), ibuf(3)
      real ftest(4), fbuf(4)
      double precision dtest(5), dbuf(5)
      character *20 cbuf(4), ccbuf(4)
      character *20 name
      character *26 date
      integer status, type, nelem
c
      data itest/1, 2, 3/
      data ftest/1.0e0, 2.0e0, 3.0e0, 4.0e0/
      data dtest/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0/
c
      cbuf(1) = 'Have'
      cbuf(2) = 'a'
      cbuf(3) = 'nice'
      cbuf(4) = 'day, Robert!'
c
      call pbeginf()
      if (.not. ma_init(MT_DBL, -1, -1)) stop 99
      if (.not. rtdb_open('fortran.db', 'empty', rtdb)) stop 1
c
      if (.not. rtdb_put(rtdb, 'Integer array', MT_INT, 3, itest))
     $     stop 2
      if (.not. rtdb_put(rtdb, 'Real array', MT_REAL, 4, ftest))
     $     stop 3
      if (.not. rtdb_put(rtdb, 'Double array', MT_DBL, 5, dtest))
     $     stop 4
      if (.not. rtdb_cput(rtdb, 'Character variable', 4, cbuf))
     $     stop 44
c
      if (.not. rtdb_print(rtdb, .true.)) stop 55
c
      if (.not. rtdb_get(rtdb, 'Integer array', MT_INT, 3, ibuf))
     $     stop 5
      write(6,*) ' Ibuf ', ibuf
      if (.not. rtdb_get(rtdb, 'Real array', MT_REAL, 4, fbuf))
     $     stop 6
      write(6,*) ' Fbuf ', fbuf
      if (.not. rtdb_get(rtdb, 'Double array', MT_DBL, 5, dbuf))
     $     stop 7
      write(6,*) ' Dbuf ', dbuf
      if (.not. rtdb_cget(rtdb, 'Character variable', 4, ccbuf))
     $     stop 8
      write(6,*) ' Cbuf 1 ', ccbuf(1)
      write(6,*) ' Cbuf 2 ', ccbuf(2)
      write(6,*) ' Cbuf 3 ', ccbuf(3)
      write(6,*) ' Cbuf 4 ', ccbuf(4)
c
      status = rtdb_first(rtdb, name)
 10   if (status) then
         if (.not. rtdb_get_info(rtdb, name, type, nelem, date)) stop 9
         write(6,*) name, type, nelem, date
         status = rtdb_next(rtdb, name)
         goto 10
      endif
c
      if (.not. rtdb_close(rtdb, 'delete')) stop 10
      call pend()
c
      end
