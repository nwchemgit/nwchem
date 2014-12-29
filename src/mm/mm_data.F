c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      block data initial_mm_data
      implicit none
#include "mm_data.fh"
c     
      data nvdw /0/
      data nr /0/

      end

      function  mm_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      logical mm_init
      integer rtdb
c
      character*180 filename
      integer       nmax,n
c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
c
      pname = "mm_init"
c
      if (nvdw .eq. 0) then
          call mm_data_load_vdw(rtdb)
      else
          goto 911
      end if
      if (nr .eq. 0) then
          call mm_data_load_coords(rtdb)
      else
          goto 911
      end if
      mm_init = .true.
      return
911   mm_init = .false.
      return
      end

      function  mm_end(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      logical mm_end
      integer rtdb
c
      character*30 pname
      pname = "mm_end"
c
      if (nvdw .gt. 0) then
        if (.not.ma_free_heap(h_bvdw)) goto 911
        if (.not.ma_free_heap(h_avdw)) goto 911
        if (.not.ma_free_heap(h_jvdw)) goto 911
        if (.not.ma_free_heap(h_ivdw)) goto 911
      end if
      mm_end = .true.
      call util_flush(6)
      return
911   mm_end = .false.
      return
      end

      subroutine mm_test()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"

      write(*,*) "handles",nvdw,h_ivdw,h_jvdw,h_avdw,h_bvdw
      end

      subroutine mm_data_load_coords(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      integer rtdb
c
      character*180 filename
      integer       nmax,n
c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
c
      pname = "mm_input_load_coords"
c
      tag="mm:coords:load"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      call mm_open_file(filename,fn)
c     ----------------------
c     count numbers of lines
c     ----------------------
      call util_xyz_nlines(nr,fn)
      write(*,*) "nr=",nr
Cc
C      if ( .not. ma_alloc_get( MT_INT, nvdw, 
C     &      'vdw scratch i',
C     &      h_ivdw, i_ivdw) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      nvdw, MA_ERR)
C
C      if ( .not. ma_alloc_get( MT_INT, nvdw, 
C     &      'vdw scratch j',
C     &      h_jvdw, i_jvdw) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      nvdw, MA_ERR)
C
C      if ( .not. ma_alloc_get( MT_DBL, nvdw, 
C     &      'vdw scratch A',
C     &      h_avdw, i_avdw) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      nvdw, MA_ERR)
C
C      if ( .not. ma_alloc_get( MT_DBL, nvdw, 
C     &      'vdw scratch B',
C     &      h_bvdw, i_bvdw) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      nvdw, MA_ERR)
C
C      call inp_init(fn,LuOut)
C      n = 0
C10    continue
C      if (.not.inp_read()) goto 30
C      message = "checking for number of fields"
C      if(inp_n_field().ne.4) goto 911
C      message = "checking for bounds"
C      if(n.gt.nvdw) goto 911
C      message = "reading vdw params"
C      n = n+1
C      write(*,*) "n=",n
C      if (.not.inp_i(int_mb(i_ivdw+n-1))) goto 911
C      if (.not.inp_i(int_mb(i_jvdw+n-1))) goto 911
C      if (.not.inp_f(dbl_mb(i_avdw+n-1))) goto 911
C      if (.not.inp_f(dbl_mb(i_bvdw+n-1))) goto 911
C      goto 10
C30    continue
C      close(fn)
C      call inp_init(LuIn,LuOut)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_data_load_vdw(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      integer rtdb
c
      character*180 filename
      integer       nmax,n
c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
c
      pname = "mm_input_load_vdw"
c
      tag="mm:vdw:load"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      call mm_open_file(filename,fn)
c     ------------------------- 
c     count numbers of vdw pairs
c     ------------------------- 
      call inp_init(fn,LuOut)
      n = 0
      do
        if (.not.inp_read()) exit
        n=n+1
      end do
      nvdw=n
      rewind(fn)
      call inp_init(fn,LuOut)
c
      max_vdw = 2500
      call util_xyz_nlines(max_vdw,fn)
      write(*,*) "nvdw,max_vdw",nvdw,max_vdw
c
      if ( .not. ma_alloc_get( MT_INT, nvdw, 
     &      'vdw scratch i',
     &      h_ivdw, i_ivdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      if ( .not. ma_alloc_get( MT_INT, nvdw, 
     &      'vdw scratch j',
     &      h_jvdw, i_jvdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      if ( .not. ma_alloc_get( MT_DBL, nvdw, 
     &      'vdw scratch A',
     &      h_avdw, i_avdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      if ( .not. ma_alloc_get( MT_DBL, nvdw, 
     &      'vdw scratch B',
     &      h_bvdw, i_bvdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      call inp_init(fn,LuOut)
      n = 0
10    continue
      if (.not.inp_read()) goto 30
      message = "checking for number of fields"
      if(inp_n_field().ne.4) goto 911
      message = "checking for bounds"
      if(n.gt.nvdw) goto 911
      message = "reading vdw params"
      n = n+1
      write(*,*) "n=",n
      if (.not.inp_i(int_mb(i_ivdw+n-1))) goto 911
      if (.not.inp_i(int_mb(i_jvdw+n-1))) goto 911
      if (.not.inp_f(dbl_mb(i_avdw+n-1))) goto 911
      if (.not.inp_f(dbl_mb(i_bvdw+n-1))) goto 911
      goto 10
30    continue
      close(fn)
      call inp_init(LuIn,LuOut)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_open_file(filename,fn)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      integer       fn
c
      character*180 buffer
      character*180 message
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
c
      pname = "mm_open_file"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
c     first try to open file in the run directory
      buffer = filename
      message = "loading vdw information from "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=10)
      goto 800
10    continue
c     now try perm directory
      call util_file_name_resolve(buffer, .false.)
      message = "loading vdw information from "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=911)
800   continue
      write(luout,*) "Successfully "//trim(message)
      write(luout,*)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end


