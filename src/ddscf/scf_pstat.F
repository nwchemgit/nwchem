      block data dscfps
#include "cscfps.fh"
      data oscfps /.false./     ! So that other codes don't break
      end
      subroutine scf_pstat_init()
      implicit none
#include "pstat.fh"
#include "cscfps.fh"
      double precision junk
c
      call pstat_init(20, 1, ' ')
c
      if (.not. pstat_allocate('scf: solve', pstat_qstat, 0,junk,
     $     ps_solve)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: search', pstat_qstat, 0, junk,
     $     ps_search)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: fock-2e', pstat_qstat, 0, junk,
     $     ps_fock_2e)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: int-1e', pstat_qstat, 0, junk,
     $     ps_int_1e)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: int-2e', pstat_qstat, 0, junk,
     $     ps_int_2e)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: fock-add', pstat_qstat, 0, junk,
     $     ps_fock_add)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: matexp', pstat_qstat, 0, junk,
     $     ps_matexp)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: diag', pstat_qstat, 0, junk,
     $     ps_diag)) call errquit('scf: scf_pstat_init', 0)
c
      if (.not. pstat_allocate('scf: total', pstat_qstat, 0, junk,
     $     ps_total)) call errquit('scf: scf_pstat_init', 0)
c
      call pstat_on(ps_total)
      oscfps = .true.
c
      end
      subroutine scf_pstat_print()
      implicit none
#include "pstat.fh"
#include "cscfps.fh"
c
      if (oscfps) then
         call pstat_off(ps_total)
         call pstat_print_all
         call pstat_terminate
         oscfps = .false.
      endif
c
      end
