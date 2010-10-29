      Logical Function ccsd(rtdb)
       implicit none
#include "errquit.fh"
C     $Id$
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "sym.fh"
c
c
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       logical status, oconverged
       character*255 title
       character*255 movecs_in
       double precision threshscf
       logical moint_driver, scf
       external moint_driver, scf
c
       Call Util_Print_Push
       Call Util_Print_RTDB_Load(RTDB, 'ccsd')
       call ecce_print_module_entry('ccsd')
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c
c     Make sure the SCF is being converged accurately enough
c
      if (.not. rtdb_get(rtdb,'scf:thresh',mt_dbl,1,threshscf))
     $     threshscf = 1.0d-4
      if (threshscf .gt. 1.0d-6) then
         threshscf = 1.0d-6
         if (.not. rtdb_put(rtdb,'scf:thresh',mt_dbl,1,threshscf))
     $        call errquit('mp2_grad: failed setting scf:thresh',0,
     &       RTDB_ERR)
         status = rtdb_delete(rtdb,'scf:converged')
      endif
c
      if (.not. scf(rtdb)) then
         ccsd = .false.
         call ecce_print_module_exit('ccsd', 'failed')
         Call Util_Print_pop
         return
      endif
c
c     Extract high level info from the data-base setting defaults
c
c ------NEW AO CCSD CODE
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('ccsd: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('ccsd: no geometry ', 0, RTDB_ERR)
      if (.not.sym_abelian_group(geom)) then
         if (.not.geom_strip_sym(geom)) call errquit
     &          ('ccsd: geom_strip_sym failed',911,
     &       GEOM_ERR)
         if (ga_nodeid().eq.0) then
            write(Luout,*) 'CCSD: Non-Abelian Symmetry Not Permitted'
            write(Luout,*) '      Symmetry has been lowered to C1'
            write(Luout,*) ' '
         endif
      endif

      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('ccsd: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('ccsd: no ao basis ', 0, RTDB_ERR)
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in)) 
     $     call util_file_name('movecs', .false.,.false.,movecs_in)
c
      call qstat_init()
c
c     init pstat
c	
      call ccsd_pstat_init(rtdb)
c
      Call ccsd_driver( rtdb, geom, basis, movecs_in, oconverged)
C
      call ccsd_pstat_print
c
      If ( .NOT. Geom_Destroy(Geom) )
     $   Call ErrQuit('CCSD: Failed to destroy geometry', Geom,
     &       GEOM_ERR)
      If ( .NOT. Bas_Destroy(Basis) )
     $   Call ErrQuit('CCSD: Failed to destroy basis', Basis,
     &       BASIS_ERR)
c
c
c
       ccsd = oconverged
       if (oconverged) then
          call ecce_print_module_exit('ccsd', 'ok')
          status = rtdb_delete(rtdb,'ccsd:converged')
       else
          call ecce_print_module_exit('ccsd', 'failed')
       endif
c
       call util_print_pop
c
       end


