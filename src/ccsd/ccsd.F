       Logical Function ccsd(rtdb)
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
c
c
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       logical status
       character*255 title
       character*255 movecs_in, movecs
       external moint_driver
       logical moint_driver,ccsd_driver
       logical do_old
c
       Call Util_Print_Push
       Call Util_Print_RTDB_Load(RTDB, 'ccsd')
       call ecce_print_module_entry('ccsd')
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c     Extract high level info from the data-base setting defaults
c
C ADRIAN
C
c ------OLD CCSD CODE
c
c$$$      do_old=.false.
c$$$      if (do_old)then
c$$$       if (.not. rtdb_cget(rtdb, 'title', 1, title))
c$$$     $      title = ' '
c$$$       if (.not. geom_create(geom, 'geometry'))
c$$$     $      call errquit('scf: geom_create?', 0)
c$$$       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
c$$$     $      call errquit('scf: no geometry ', 0)
c$$$       if (.not. bas_create(basis, 'ao basis'))
c$$$     $      call errquit('scf: bas_create?', 0)
c$$$       if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
c$$$     $      call errquit('scf: no ao basis ', 0)
c$$$       movecs = 'movecs'
c$$$       if (.not. context_rtdb_match(rtdb, 'input mo vector file',
c$$$     $      movecs_in)) movecs_in = movecs
c$$$       status = moint_driver( rtdb, geom, basis, movecs_in )
c$$$      endif
c
c ------NEW AO CCSD CODE
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('ccsd: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('ccsd: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('ccsd: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('ccsd: no ao basis ', 0)
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in)) 
     $     call util_file_name('movecs', .false.,.false.,movecs_in)
c
      status = ccsd_driver( rtdb, geom, basis, movecs_in )
C
      If ( .NOT. Geom_Destroy(Geom) )
     $   Call ErrQuit('CCSD: Failed to destroy geometry', Geom)
      If ( .NOT. Bas_Destroy(Basis) )
     $   Call ErrQuit('CCSD: Failed to destroy basis', Basis)
c
c
c
       ccsd = .TRUE.
       call ecce_print_module_exit('ccsd', 'ok')
       Return
       end


