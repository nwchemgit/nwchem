       subroutine mjdtest(rtdb)
*
* $Id$
*
       implicit none
 #include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       logical status
       character*255 movecs_in, movecs_out
c       double precision charge,energy
       double precision energy
c
       if (.not. context_push('mjdtest'))
     $      call errquit('mjdtest: context_push failed',0, UNKNOWN_ERR)
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
      call util_print_centered(6, 'NWCHEM CPHF Module', 40, .true.)
      write(6,*)
      write(6,*)
c
c     Extract high level info from the data-base setting defaults
c
c     
      call mjd_get_db_handles( rtdb, geom, basis, movecs_in,
     $     movecs_out )
c      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
c     $     charge = 0.0d0
c
      call cphf_poliz(rtdb, geom, basis, movecs_in )
c
       if (.not. context_pop('mjdtest'))
     $     call errquit('context_pop failed',0, UNKNOWN_ERR)
       end

       subroutine mjd_get_db_handles( rtdb, geom, basis, movecs_in,
     $                                movecs_out )
       implicit none
 #include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
       integer rtdb
       integer geom, basis
       character*(*) movecs_in, movecs_out

       if (.not. geom_create(geom, 'geometry'))
     $     call errquit('get_db_handles: geom_create?', 0, GEOM_ERR)
       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('get_db_handles: no geometry ', 0, RTDB_ERR)
       if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('get_db_handles: bas_create?', 0, BASIS_ERR)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('get_db_handles: no ao basis ', 0, RTDB_ERR)
       if (.not. context_rtdb_match(rtdb, 'input mo vector file',
     $     movecs_in)) movecs_in = 'movecs'
c
       if (.not. geom_print(geom)) call errquit('cphf', 0, GEOM_ERR)
       if (.not. bas_print(basis)) call errquit('cphf', 0, BASIS_ERR)
c
       return
       end

