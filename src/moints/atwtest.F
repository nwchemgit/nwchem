       logical function atwtest(rtdb)
C$Id$
       implicit none
#include "errquit.fh"
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
       character*255 title
       character*255 movecs_in, movecs
c
c
       if (.not. context_push('atwtest'))
     $      call errquit('atwtest: context_push failed',0,
     &       UNKNOWN_ERR)
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf: no geometry ', 0, GEOM_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('scf: no ao basis ', 0, RTDB_ERR)
      movecs = 'movecs'
      if (.not. context_rtdb_match(rtdb, 'input mo vector file',
     $     movecs_in)) movecs_in = movecs
c
c
c
      call moint_driver( rtdb, geom, basis, movecs_in )
c
c
c
      if (.not.geom_destroy(geom))
     $  call errquit('cannot destroy geom handle?',0, GEOM_ERR)
      if (.not.bas_destroy(basis))
     $  call errquit('cannot destroy basis handle?',0, BASIS_ERR)
      if (.not. context_pop('atwtest'))
     $   call errquit('context_pop failed',0, MEM_ERR)

      atwtest = .true.
      return
      end




