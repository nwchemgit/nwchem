c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      block data initial_mm_data
      implicit none
#include "mm_data.fh"
c     
      data ntot /0/
      data nvdw /0/
      data nmm /0/
      data nqm /0/

      end

      function  mmi_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
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
#include "mm_data.fh"
#include "mafdecls.fh"

      integer rtdb
      double precision e
      call mm_vdw_add_energy(rtdb,e)
      end

      subroutine mm_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mm_data.fh"
#include "mafdecls.fh"

      integer rtdb
      double precision e
      integer n
      double precision g(3,n)
c
      integer i,j
      character*30 pname

      call mm_vdw_add_egrad(rtdb,e,n,g)

      write(6,1000) "QM + MM/VDW",
     $     'x','y','z','x','y','z'
      do i = 1, nqm
         write(6,2000) i, byte_mb(i_t+16*(i-1)),
     $        (dbl_mb(i_r+3*(i-1)+j),j=0,2),
     $        g(1,i),g(2,i),g(3,i)
      enddo
      write(6,*)



 1000 format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $     'coordinates',
     $     24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000 format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
      write(6,*)
      call util_flush(6)

      return
911   call errquit("error "//trim(pname),0,-1)
      end

