      subroutine ccsd_pstat_init(rtdb)
c  $Id: ccsd_pstat.F,v 2.1 2004-12-16 02:40:47 edo Exp $
      implicit none
#include "pstat.fh"
#include "ccsdps.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer junk
c
      if (.not. rtdb_get(rtdb, 'ccsd:pstat',
     $     MT_LOG, 1, occsdps)) occsdps = .false.
c
      if (.not. occsdps) return
c
      if (.not. pstat_allocate('ccsd:aomo', pstat_qstat, 0, junk,
     $     ps_aomo)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:aoccsd', pstat_qstat, 0, junk,
     $     ps_aoccsd)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:tripls', pstat_qstat, 0, junk,
     $     ps_tripls)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:graddrv', pstat_qstat, 0, junk,
     $     ps_graddrv)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:itm', pstat_qstat, 0, junk,
     $     ps_itm)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:pampt', pstat_qstat, 0, junk,
     $     ps_pampt)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:idx1', pstat_qstat, 0, junk,
     $     ps_idx1)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:idx2', pstat_qstat, 0, junk,
     $     ps_idx2)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:idx34', pstat_qstat, 0, junk,
     $     ps_idx34)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:t2eri', pstat_qstat, 0, junk,
     $     ps_t2eri)) call errquit('ccsd: ccsd_pstat_init', 0,0)
c       
      call pstat_on(ps_ccsd)
      occsdps = .true.
c
      end
      subroutine ccsd_pstat_print()
      implicit none
#include "pstat.fh"
#include "ccsdps.fh"
c     
      if (occsdps) then
         call pstat_off(ps_ccsd)
         call pstat_print_all
         occsdps = .false.
c
c     delete the SCF timers
c
         if(.not.pstat_free(ps_aomo))call errquit('ccsd_pstat?', 0,0)
         if(.not.pstat_free(ps_aoccsd))call errquit('ccsd_pstat?', 0,0)
         if(.not.pstat_free(ps_tripls))call errquit('ccsd_pstat?',0,0)
         if(.not.pstat_free(ps_graddrv))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_itm))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_pampt))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_idx1))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_idx2))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_idx34))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_t2eri))call errquit('ccsd_pstat',0,0)
c
      endif
c
      end
