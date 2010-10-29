      block data Duccsdtps
#include "cuccsdtps.fh"
C$Id$
      data ouccsdtps /.false./     ! So that other codes do not break
c
      end
      subroutine uccsdt_pstat_init(rtdb)
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "cuccsdtps.fh"
#include "cuccsdtP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer junk
      external duccsdtps ! For T3D linker
c
      if (.not. rtdb_get(rtdb, 'uccsdt:pstat',
     $     MT_LOG, 1, ouccsdtps)) ouccsdtps = .false.
c
      if (.not. ouccsdtps) return
c
      if (.not. pstat_allocate('uccsdt: triples', 
     $     pstat_qstat, 0, junk, ps_triples)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: triples_super_blocking', 
     $     pstat_qstat, 0, junk, ps_triples_super_blocking)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: triples_sub_blocking', 
     $     pstat_qstat, 0, junk, ps_triples_sub_blocking)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: triples_moints', 
     $     pstat_qstat, 0, junk, ps_triples_moints)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: moints_shpair_one_index', 
     $     pstat_qstat, 0, junk, ps_moints_shpair_one_index)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: moints_second_index', 
     $     pstat_qstat, 0, junk, ps_moints_second_index)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_1a', 
     $     pstat_qstat, 0, junk, ps_mo_fin_1a)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_1b', 
     $     pstat_qstat, 0, junk, ps_mo_fin_1b)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_2a', 
     $     pstat_qstat, 0, junk, ps_mo_fin_2a)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_2b', 
     $     pstat_qstat, 0, junk, ps_mo_fin_2b)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_3a', 
     $     pstat_qstat, 0, junk, ps_mo_fin_3a)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_3b', 
     $     pstat_qstat, 0, junk, ps_mo_fin_3b)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_4a', 
     $     pstat_qstat, 0, junk, ps_mo_fin_4a)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_4b', 
     $     pstat_qstat, 0, junk, ps_mo_fin_4b)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_5a', 
     $     pstat_qstat, 0, junk, ps_mo_fin_5a)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_5b', 
     $     pstat_qstat, 0, junk, ps_mo_fin_5b)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_6', 
     $     pstat_qstat, 0, junk, ps_mo_fin_6)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_7', 
     $     pstat_qstat, 0, junk, ps_mo_fin_7)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_8', 
     $     pstat_qstat, 0, junk, ps_mo_fin_8)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_9', 
     $     pstat_qstat, 0, junk, ps_mo_fin_9)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_10', 
     $     pstat_qstat, 0, junk, ps_mo_fin_10)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_11', 
     $     pstat_qstat, 0, junk, ps_mo_fin_11)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_12', 
     $     pstat_qstat, 0, junk, ps_mo_fin_12)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: mo_fin_13', 
     $     pstat_qstat, 0, junk, ps_mo_fin_13)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: moints_triples_final', 
     $     pstat_qstat, 0, junk, ps_moints_triples_final)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_ga_dist', 
     $     pstat_qstat, 0, junk, ps_final_ga_dist)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_ga_get', 
     $     pstat_qstat, 0, junk, ps_final_ga_get)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_dgemm1', 
     $     pstat_qstat, 0, junk, ps_final_dgemm1)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_dgemm2', 
     $     pstat_qstat, 0, junk, ps_final_dgemm2)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_u_trans', 
     $     pstat_qstat, 0, junk, ps_final_u_trans)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_dgemm3', 
     $     pstat_qstat, 0, junk, ps_final_dgemm3)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_dgemm4', 
     $     pstat_qstat, 0, junk, ps_final_dgemm4)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: final_ga_acc', 
     $     pstat_qstat, 0, junk, ps_final_ga_acc)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: triples_amplitudes', 
     $     pstat_qstat, 0, junk, ps_triples_amplitudes)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: triples_sub_blocked', 
     $     pstat_qstat, 0, junk, ps_triples_sub_blocked)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: triples_sub_blocking_end', 
     $     pstat_qstat, 0, junk, ps_triples_sub_blocking_end)) 
     $     call errquit('uccsdt: uccsdt_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('uccsdt: total', pstat_qstat, 0, junk,
     $     ps_total)) call errquit('uccsdt: uccsdt_pstat_init', 0,
     &       MEM_ERR)
c
      call pstat_on(ps_total)
      ouccsdtps = .true.
c
      end
      subroutine uccsdt_pstat_print()
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "cuccsdtps.fh"
#include "cuccsdtP.fh"
c     
      if (ouccsdtps) then
         call pstat_off(ps_total)
         call pstat_print_all
         ouccsdtps = .false.
c
c     delete the UCCSDT timers
c
         if (.not. pstat_free(ps_triples))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_triples_super_blocking))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_triples_sub_blocking))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_triples_moints))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_moints_shpair_one_index))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_moints_second_index))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_1a))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_1b))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_2a))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_2b))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_3a))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_3b))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_4a))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_4b))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_5a))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_5b))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_6))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_7))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_8))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_9))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_10))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_11))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_12))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_mo_fin_13))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_moints_triples_final))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_ga_dist))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_ga_get))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_dgemm1))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_dgemm2))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_u_trans))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_dgemm3))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_dgemm4))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_final_ga_acc))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_triples_amplitudes))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_triples_sub_blocked))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_triples_sub_blocking_end))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
         if (.not. pstat_free(ps_total))
     &        call errquit('uccsdt_pstat?', 0, MEM_ERR)
c
      endif
c
      end
