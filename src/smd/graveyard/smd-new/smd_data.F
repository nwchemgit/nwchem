*
* $Id$
*
      function smd_data_create(srtdb,dname,ndim,dtype)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "srtdb.fh"
      integer srtdb
      character*(*) dname
      integer ndim
      integer dtype
      logical smd_data_create
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,handle,id,datum(4)

      pname = "smd_data_create"


      smd_data_create = .false.
      buffer = dname
      if(.not.ma_alloc_get(dtype,
     +   ndim,
     +   dname,
     +   handle,
     +   id)) then
         call smd_message_post(pname//'memory allocation '//buffer)
         call smd_message_post(pname//'memory allocation1 '//buffer)
         call smd_message_post(pname//'memory allocation2 '//buffer)
         call smd_message_post(pname//'memory allocation3 '//buffer)
         return
      end if

      if(dtype.eq.MT_DBL) then
        call dfill(ndim,
     +             0.d0,
     +             dbl_mb(id),
     +             1)
      end if

      if(dtype.eq.MT_INT) then
        call ifill(ndim,
     +             0,
     +             int_mb(id),
     +             1)
      end if
      if(dtype.eq.MT_LOG) then
        do i=1,ndim
          log_mb(id+i-1) = .false.
        end do
      end if
c     should probably fill logical and char arrays

      datum(1) = dtype
      datum(2) = ndim
      datum(3) = id
      datum(4) = handle
      if (.not.srtdb_put(srtdb,dname,MT_INT,4,datum)) then
        call smd_message_post(pname//'srtdb_put '//buffer)
        return
      end if

      smd_data_create = .true.
      return
      end

      function smd_data_create_get(srtdb,dname,ndim,dtype,id)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "srtdb.fh"
      integer srtdb
      character*(*) dname
      integer ndim
      integer dtype
      integer id
      logical smd_data_create_get
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,handle,datum(4)

      pname = "smd_data_create_get"


      smd_data_create_get = .false.
      buffer = dname
      if(.not.ma_alloc_get(dtype,
     +   ndim,
     +   dname,
     +   handle,
     +   id)) then
         call smd_message_post(pname//'memory allocation '//buffer)
         call smd_message_post(pname//'memory allocation1 '//buffer)
         call smd_message_post(pname//'memory allocation2 '//buffer)
         call smd_message_post(pname//'memory allocation3 '//buffer)
         return
      end if

      if(dtype.eq.MT_DBL) then
        call dfill(ndim,
     +             0.d0,
     +             dbl_mb(id),
     +             1)
      end if

      if(dtype.eq.MT_INT) then
        call ifill(ndim,
     +             0,
     +             int_mb(id),
     +             1)
      end if
      if(dtype.eq.MT_LOG) then
        do i=1,ndim
          log_mb(id+i-1) = .false.
        end do
      end if
c     should probably fill logical and char arrays

      datum(1) = dtype
      datum(2) = ndim
      datum(3) = id
      datum(4) = handle
      if (.not.srtdb_put(srtdb,dname,MT_INT,4,datum)) then
        call smd_message_post(pname//'srtdb_put '//buffer)
        return
      end if

      smd_data_create_get = .true.
      return
      end

      function smd_data_get_index(srtdb,dname,descriptor)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "srtdb.fh"
      integer srtdb
      character*(*) dname
      integer descriptor
      logical smd_data_get_index
c
      logical result
      character*30 pname
      character*72 buffer
      integer datum(3)

      pname = "smd_data_get_index"

      smd_data_get_index = .false.

      if (.not.srtdb_get(srtdb,dname,MT_INT,4,datum)) then
         call smd_message_post(pname//" srtdb_get")
         return
      end if

      descriptor = datum(3)
      smd_data_get_index = .true.
      return
      end

      subroutine smd_data_print(srtdb,dname,un)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "srtdb.fh"
      integer srtdb
      character*(*) dname
      integer un
c
      logical result
      character*30 pname
      character*72 buffer
      integer i,j,datum(3)

      pname = "smd_data_print"

      if (.not.srtdb_get(srtdb,dname,MT_INT,4,datum)) then
         call errquit(pname//" srtdb_get",0,0)
         return
      end if

      write(un,*) "SMD DATA ARRAY"
      write(un,*) "name:",dname
      
      if(datum(1).eq.MT_INT) then
        write(un,*) "size:",datum(2)
        do i=1,datum(2)
          write(un,*) int_mb(datum(3)+i-1)
        end do
      else if (datum(1).eq.MT_DBL) then
        write(un,*) "size:",datum(2)
        do i=1,datum(2)
          write(un,*) dbl_mb(datum(3)+i-1)
        end do
      else if (datum(1).eq.MT_LOG) then
        write(un,*) "size:",datum(2)
        do i=1,datum(2)
          write(un,'(A16)') log_mb(datum(3)+i-1)
        end do
      else if (datum(1).eq.MT_BYTE) then
        write(un,*) "size:",datum(2)/16
c        call smd_ma_chararray_print(un,datum(2)/16,byte_mb(datum(3)))
        do i=1,datum(2)/16
          write(un,'(16A1)') (byte_mb(datum(3)+16*(i-1)+j-1),j=1,16)
        end do
      end if

      return
      end

      subroutine smd_ma_chararray_print(un,n,a)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "srtdb.fh"
      integer n
      integer un
      character*16 a(n)
c
      integer i
      do i=1,n
        write(un,*) a(i)
      end do

      end subroutine

      function smd_data_destroy(srtdb,dname)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "srtdb.fh"
      integer srtdb
      character*(*) dname
      logical smd_data_destroy
c
      logical result
      character*30 pname
      character*72 buffer
      integer datum(4)

      pname = "smd_data_destroy"

      smd_data_destroy = .false.

      if (.not.srtdb_get(srtdb,dname,MT_INT,4,datum)) then
         call smd_message_post(pname//" srtdb_get")
         return
      end if

      if(.not.ma_free_heap(datum(4))) then
         call smd_message_post(pname//" ma_free_heap")
         return
      end if

      smd_data_destroy = .true.
      return
      end
