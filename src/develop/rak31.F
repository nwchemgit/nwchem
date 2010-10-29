*
* $Id$
*
      logical function rak_justrunvib(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "geom.fh"
      integer rtdb
c
      character*255 filehess
      logical ignore
      integer geom, nat
      ignore = rtdb_parallel(.false.)
      if ((ga_nodeid()).eq.0) then
         if (.not. rtdb_cget(rtdb, 'task:hessian file name', 1,
     $        filehess)) then
            call util_file_name('hess',  .false., .false.,filehess)
         endif
c
c     create/load reference geometry
c
        if (.not.geom_create(geom,'geometry')) call errquit
     $      ('task_freq:geom_create failed?',1, GEOM_ERR)
        if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $      call errquit
     $      ('task_freq:geom_rtdb_load failed?',2, RTDB_ERR)
        if (.not. geom_ncent(geom,nat)) call errquit
     $      ('task_freq:geom_ncent failed?',3, GEOM_ERR)
        if (.not. geom_destroy(geom)) call errquit
     $      ('task_freq:geom_destroy failed?',911, GEOM_ERR)
        call vib_vib(rtdb,filehess,.true.,
     $      0,.false.,0,.false.,nat)
      endif
      call ga_sync()
      ignore = rtdb_parallel(.true.)
      rak_justrunvib = .true.
      end
