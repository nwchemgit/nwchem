C $Id$
      subroutine drdy_nwc_get_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
*::passed::
      integer rtdb
      double precision energy
*::local::
      logical task_status
*
      energy = 0.0d00
*      
      if (.not.rtdb_get(rtdb,'task:status', mt_log, 1, task_status))
     &      call errquit(
     &      'drdy_nwc_get_energy: rtdb read status failed',911,
     &       RTDB_ERR)
      if (.not.task_status) then
        write(luout,*) ' task status is failed '
        call errquit('drdy_nwc_get_energy: fatal error',911, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy)) then
        write(luout,*)' failed to read task:energy '
        call errquit('drdy_nwc_get_energy: fatal error',911, RTDB_ERR)
      endif
      end
      subroutine drdy_nwc_get_coords(rtdb,xyz)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
*::passed::
      integer rtdb
      double precision xyz(3,*)
*::local::
      integer geom
      logical task_status
*
      if (.not.rtdb_get(rtdb,'task:status', mt_log, 1, task_status))
     &      call errquit(
     &      'drdy_nwc_get_coords: rtdb read status failed',911,
     &       RTDB_ERR)
      if (.not.task_status) then
        write(luout,*) ' task status is failed '
        call errquit('drdy_nwc_get_coords: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_create(geom,'geometry')) then
        write(luout,*)' failed to create geometry'
        call errquit('drdy_nwc_get_coords: fatal error',911, GEOM_ERR)
      endif
      if (.not.geom_rtdb_load(rtdb, geom,'geometry')) then
        write(luout,*)' failed to read geometry'
        call errquit('drdy_nwc_get_coords: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_cart_coords_get(geom,xyz)) then
        write(luout,*)' failed to get geometry coords'
        call errquit('drdy_nwc_get_coords: fatal error',911,
     &       GEOM_ERR)
      endif
      if (.not.geom_destroy(geom)) then
        write(luout,*)' failed to destroy geometry'
        call errquit('drdy_nwc_get_coords: fatal error',911,
     &       GEOM_ERR)
      endif
      end
      subroutine drdy_nwc_get_gradient(rtdb,grad)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb
      double precision grad(*)
*
      logical task_status
      integer geom
      integer mycent
*
      if (.not.rtdb_get(rtdb,'task:status', mt_log, 1, task_status))
     &      call errquit(
     &      'drdy_nwc_get_gradient: rtdb read status failed',911,
     &       RTDB_ERR)
      if (.not.task_status) then
        write(luout,*) ' task status is failed '
        call errquit('drdy_nwc_get_gradient: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_create(geom,'geometry')) then
        write(luout,*)' failed to create geometry'
        call errquit('drdy_nwc_get_gradient: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_rtdb_load(rtdb, geom,'geometry')) then
        write(luout,*)' failed to read geometry'
        call errquit('drdy_nwc_get_gradient: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_ncent(geom,mycent)) then
        write(luout,*)' failed to get mycent geometry'
        call errquit('drdy_nwc_get_gradient: fatal error',911,
     &       GEOM_ERR)
      endif
      if (.not.rtdb_get(rtdb,'task:gradient',mt_dbl,3*mycent,grad)) then
        write(luout,*)' failed to get task gradient'
        call errquit('drdy_nwc_get_gradient: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_destroy(geom)) then
        write(luout,*)' failed to destroy geometry'
        call errquit('drdy_nwc_get_gradient: fatal error',911,
     &       RTDB_ERR)
      endif
*
      end
      subroutine drdy_nwc_get_hessian(rtdb,hess)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb
      double precision hess(*)
*
      integer geom, mycent, nhesst, iii
      logical task_status
      character*255 hess_file
      double precision dbl_tmp
*
      if (.not.rtdb_get(rtdb,'task:status', mt_log, 1, task_status))
     &      call errquit(
     &      'drdy_nwc_get_hessian: rtdb read status failed',911,
     &       RTDB_ERR)
      if (.not.task_status) then
        write(luout,*) ' task status is failed '
        call errquit('drdy_nwc_get_hessian: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_create(geom,'geometry')) then
        write(luout,*)' failed to create geometry'
        call errquit('drdy_nwc_get_hessian: fatal error',911,
     &       GEOM_ERR)
      endif
      if (.not.geom_rtdb_load(rtdb, geom,'geometry')) then
        write(luout,*)' failed to read geometry'
        call errquit('drdy_nwc_get_hessian: fatal error',911,
     &       RTDB_ERR)
      endif
      if (.not.geom_ncent(geom,mycent)) then
        write(luout,*)' failed to get mycent geometry'
        call errquit('drdy_nwc_get_hessian: fatal error',911, GEOM_ERR)
      endif
      nhesst = 3*mycent
      nhesst = (nhesst*(nhesst+1))/2
      call util_file_name('hess',  .false., .false.,hess_file)
      open(unit=69,file=hess_file,form='formatted',status='old',
     &      err=99900,access='sequential')
      do iii = 1,nhesst
        read(69,*,err=99901,end=99902)dbl_tmp
        hess(iii) = dbl_tmp
      enddo
      close(unit=69,status='keep')
*
      if (.not.geom_destroy(geom)) then
        write(luout,*)' failed to destroy geometry'
        call errquit('drdy_nwc_get_hessian: fatal error',911, GEOM_ERR)
      endif
*
      return
*
99900 continue
      write(luout,*)'hess_file => ',hess_file
      write(luout,*)
     &      'drdy_nwc_get_hessian: error opening file: "hess_file"'
      call errquit('drdy_nwc_get_hessian: fatal error',911,
     &       UNKNOWN_ERR)
99901 continue
      write(luout,*)'hess_file => ',hess_file
      write(luout,*)
     &      'drdy_nwc_get_hessian: error reading file: "hess_file"'
      call errquit('drdy_nwc_get_hessian: fatal error',911,
     &       UNKNOWN_ERR)
99902 continue
      write(luout,*)'hess_file => ',hess_file
      write(luout,*)
     &      'drdy_nwc_get_hessian: unexpected EOF when reading file:',
     &      ' "hess_file"'
      call errquit('drdy_nwc_get_hessian: fatal error',911,
     &       UNKNOWN_ERR)
      end
