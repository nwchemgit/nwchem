       subroutine mjdtest(rtdb)
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
c
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       logical status
       character*255 movecs_in, movecs_out
       double precision charge,energy
c
c
       if (.not. context_push('mjdtest'))
     $      call errquit('mjdtest: context_push failed',0)
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c     Extract high level info from the data-base setting defaults
c
c     
       call mjd_get_db_handles( rtdb, geom, basis, movecs_in,
     $                          movecs_out )
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
c
c
c
      call cphf_poliz( rtdb, geom, basis, movecs_in )
c
c
c
       if (.not. context_pop('mjdtest'))
     $     call errquit('context_pop failed',0)
       end










       subroutine mjd_get_db_handles( rtdb, geom, basis, movecs_in,
     $                                movecs_out )
       implicit none
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
c
c
       integer rtdb
       integer geom, basis
       character*(*) movecs_in, movecs_out

       if (.not. geom_create(geom, 'geometry'))
     $     call errquit('get_db_handles: geom_create?', 0)
       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('get_db_handles: no geometry ', 0)
       if (.not. bas_create(basis, 'mo basis'))
     $     call errquit('get_db_handles: bas_create?', 0)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $     call errquit('get_db_handles: no mo basis ', 0)
       if (.not. context_rtdb_match(rtdb, 'input mo vector file',
     $     movecs_in)) movecs_in = 'movecs'
       if (.not. context_rtdb_match(rtdb, 'output mo vector file',
     $     movecs_out)) movecs_out = 'movecs'


       return
       end

