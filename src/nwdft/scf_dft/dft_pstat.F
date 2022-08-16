      subroutine dft_pstat_init(rtdb)
      implicit none
#include "pstat.fh"
#include "dftps.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer junk
c
      if (.not. rtdb_get(rtdb, 'dft:pstat',
     $     MT_LOG, 1, odftps)) odftps = .false.
c
      if (.not. odftps) return
c
      if (.not. pstat_allocate('dft: 1-e', pstat_qstat, 0, junk,
     $     ps_onel)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft: guess', pstat_qstat, 0, junk,
     $     ps_guess)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:incore', pstat_qstat, 0, junk,
     $     ps_incore)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:chol  ', pstat_qstat, 0, junk,
     $     ps_chol)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:inv   ', pstat_qstat, 0, junk,
     $     ps_inv)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:2ceri', pstat_qstat, 0, junk,
     $     ps_2cer)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:  xc', pstat_qstat, 0, junk,
     $     ps_xc)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:fck2e', pstat_qstat, 0, junk,
     $     ps_f2e)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:gvxc', pstat_qstat, 0, junk,
     $     ps_getvxc)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:qv0b', pstat_qstat, 0, junk,
     $     ps_qv0b)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:xmrge', pstat_qstat, 0, junk,
     $     ps_xmrge)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:xcopy', pstat_qstat, 0, junk,
     $     ps_xcopy)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:ssw', pstat_qstat, 0, junk,
     $     ps_ssw)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:xcrho', pstat_qstat, 0, junk,
     $     ps_xcrho)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:tabcd', pstat_qstat, 0, junk,
     $     ps_tabcd)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:ebf', pstat_qstat, 0, junk,
     $     ps_ebf)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:excf', pstat_qstat, 0, junk,
     $     ps_excf)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:xcfit', pstat_qstat, 0, junk,
     $     ps_fitv)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:diag', pstat_qstat, 0, junk,
     $     ps_diag)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:fitcd', pstat_qstat, 0, junk,
     $     ps_fitcd)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:vcoul', pstat_qstat, 0, junk,
     $     ps_vcoul)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:bld12', pstat_qstat, 0, junk,
     $     ps_bld12)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:diis', pstat_qstat, 0, junk,
     $     ps_diis)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:fockb', pstat_qstat, 0, junk,
     $     ps_fockb)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:dgemm', pstat_qstat, 0, junk,
     $     ps_dgemm)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:bdgem', pstat_qstat, 0, junk,
     $     ps_bdgemm)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:scfend', pstat_qstat, 0, junk,
     $     ps_scfend)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:scf', pstat_qstat, 0, junk,
     $     ps_scfdft)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:cpks_coef', pstat_qstat, 0, junk,
     $     ps_cpkscoef)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('dft:total', pstat_qstat, 0, junk,
     $     ps_dft)) call errquit('dft: dft_pstat_init', 0,0)
c       
      call pstat_on(ps_dft)
      odftps = .true.
c
      end
      subroutine dft_pstat_print()
      implicit none
#include "pstat.fh"
#include "dftps.fh"
c     
      if (odftps) then
         call pstat_off(ps_dft)
         call pstat_print_all
         odftps = .false.
c
c     delete the SCF timers
c
         if (.not. pstat_free(ps_onel))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_guess))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_xc))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_f2e))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_getvxc))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_qv0b))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_diis))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_fockb))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_diag))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_dgemm))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_bdgemm))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_fitcd))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_incore))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_chol))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_inv))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_2cer))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_vcoul))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_scfdft))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_bld12))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_ssw))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_xcrho))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_tabcd))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_ebf))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_excf))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_dft))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_xmrge))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_xcopy))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_scfend))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_fitv))call errquit('dft_pstat?', 0,0)
c
      endif
c
      end
c $Id$
