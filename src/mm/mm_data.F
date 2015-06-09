c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      function  mmi_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      logical mmi_init
      integer rtdb
c
      logical ignore

      ignore = MA_set_hard_fail(.true.)
      ignore=MA_set_auto_verify(.true.)
      call mm_coords_init(rtdb)
      call mm_vdw_init(rtdb)

      mmi_init = .true.
      return
911   mmi_init = .false.
      return
      end

      function  mmi_end()
      implicit none
      logical mmi_end
c
      character*30 pname
      pname = "mm_end"
      call mm_coords_end()
      call mm_vdw_end()
      mmi_end = .true.
      return
911   mmi_end = .false.
      return
      end

      subroutine mm_test(n,t,c)
      implicit none

      integer n
      character*(16) t(n)
      double precision c(3,n)
      integer i

      do i=1,n
         write(6,*) t(i),c(1,i),c(2,i),c(3,i)
      end do
      
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
      message = "opening file "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=10)
      goto 800
10    continue
c     now try perm directory
      call util_file_name_resolve(buffer, .false.)
      message = "opening file "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=911)
800   continue
      write(luout,*) "Successfully "//trim(message)
      write(luout,*)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"

      integer rtdb
      double precision e
c
      logical rtdb_mode
      integer master
c
      master=0
      rtdb_mode = rtdb_parallel(.false.)
      if(ga_nodeid().eq.master) then
        call mm_vdw_add_energy(rtdb,e)
      end if
      call ga_sync()
      call ga_brdcst(msg_smd,e,
     >     ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_sync()
      rtdb_mode = rtdb_parallel(rtdb_mode)
      end

      subroutine mm_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"

      integer rtdb
      double precision e
      integer n
      double precision g(3,n)
c
      integer i,j
      character*30 pname
      integer master
      logical rtdb_mode

      master=0
      rtdb_mode = rtdb_parallel(.false.)
      if(ga_nodeid().eq.master) then
        call mm_vdw_add_egrad(rtdb,e,n,g)
      end if
      call ga_sync()
      call ga_brdcst(msg_smd,e,
     >     ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_brdcst(msg_smd,g,
     >     3*n*ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_sync()
      rtdb_mode = rtdb_parallel(rtdb_mode)

      return
911   call errquit("error "//trim(pname),0,-1)
      end

