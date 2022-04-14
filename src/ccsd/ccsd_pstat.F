      subroutine ccsd_pstat_init(rtdb)
c  $Id$
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
      if (.not. pstat_allocate('ccsd:trpdrv', pstat_qstat, 0, junk,
     $     ps_trpdrv)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:doooo', pstat_qstat, 0, junk,
     $     ps_doooo)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:dovvv', pstat_qstat, 0, junk,
     $     ps_dovvv)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:doxxx', pstat_qstat, 0, junk,
     $     ps_doxxx)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:gpumove', pstat_qstat, 0, junk,
     $     ps_gpumove)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:accwait', pstat_qstat, 0, junk,
     $     ps_accwait)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:z2pm', pstat_qstat, 0, junk,
     $     ps_z2pm)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:hz2pm', pstat_qstat, 0, junk,
     $     ps_hz2pm)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:zvecs', pstat_qstat, 0, junk,
     $     ps_zvecs)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:pdiis', pstat_qstat, 0, junk,
     $     ps_pdiis)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:ht2pm', pstat_qstat, 0, junk,
     $     ps_ht2pm)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:t2pm', pstat_qstat, 0, junk,
     $     ps_t2pm)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:sxy', pstat_qstat, 0, junk,
     $     ps_sxy)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:itrdrv', pstat_qstat, 0, junk,
     $     ps_itrdrv)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('ccsd:tengy', pstat_qstat, 0, junk,
     $     ps_tengy)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('t2eriw', pstat_qstat, 0, junk,
     $     ps_t2eriw)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('t2erin', pstat_qstat, 0, junk,
     $     ps_t2erin)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('t2erih', pstat_qstat, 0, junk,
     $     ps_t2erih)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('rdtrpo', pstat_qstat, 0, junk,
     $     ps_rdtrpo)) call errquit('ccsd: ccsd_pstat_init', 0,0)
      if (.not. pstat_allocate('trpmos', pstat_qstat, 0, junk,
     $     ps_trpmos)) call errquit('ccsd: ccsd_pstat_init', 0,0)
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
         if(.not.pstat_free(ps_trpdrv))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_accwait))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_gpumove))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_doxxx))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_dovvv))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_doooo))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_z2pm))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_hz2pm))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_zvecs))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_pdiis))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_ht2pm))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_t2pm))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_sxy))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_itrdrv))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_tengy))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_t2eriw))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_t2erin))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_t2erih))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_rdtrpo))call errquit('ccsd_pstat',0,0)
         if(.not.pstat_free(ps_trpmos))call errquit('ccsd_pstat',0,0)
c
      endif
c
      end
