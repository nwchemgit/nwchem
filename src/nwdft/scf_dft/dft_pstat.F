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
      if (.not. pstat_allocate('dft:  xc', pstat_qstat, 0, junk,
     $     ps_xc)) call errquit('dft: dft_pstat_init', 0)
      if (.not. pstat_allocate('dft:diag', pstat_qstat, 0, junk,
     $     ps_diag)) call errquit('dft: dft_pstat_init', 0)
      if (.not. pstat_allocate('dft:dgemm', pstat_qstat, 0, junk,
     $     ps_dgemm)) call errquit('dft: dft_pstat_init', 0)
      if (.not. pstat_allocate('dft:fitcd', pstat_qstat, 0, junk,
     $     ps_fitcd)) call errquit('dft: dft_pstat_init', 0)
      if (.not. pstat_allocate('dft:vcoul', pstat_qstat, 0, junk,
     $     ps_vcoul)) call errquit('dft: dft_pstat_init', 0)
      if (.not. pstat_allocate('dft:total', pstat_qstat, 0, junk,
     $     ps_dft)) call errquit('dft: dft_pstat_init', 0)
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
         if (.not. pstat_free(ps_xc))call errquit('dft_pstat?', 0)
         if (.not. pstat_free(ps_diag))call errquit('dft_pstat?', 0)
         if (.not. pstat_free(ps_dgemm))call errquit('dft_pstat?', 0)
         if (.not. pstat_free(ps_fitcd))call errquit('dft_pstat?', 0)
         if (.not. pstat_free(ps_vcoul))call errquit('dft_pstat?', 0)
         if (.not. pstat_free(ps_dft))call errquit('dft_pstat?', 0)
c
      endif
c
      end
