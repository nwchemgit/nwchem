      subroutine ccsd_t_pstat_init(rtdb)
      implicit none
#include "pstat.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "ccsd_t_ps.fh"
      integer rtdb
      integer junk
c
      if (.not. rtdb_get(rtdb, 'tce:pstat',
     $     MT_LOG, 1, otceps)) otceps = .false.
c
      call util_flush(6)
      if (.not. otceps) return
c
      if (.not. pstat_allocate('tce:cctdbl', pstat_qstat, 0, junk,
     $     ps_cctdbl)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('tce:sdkrnl', pstat_qstat, 0, junk,
     $     ps_sdkrnl)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('tce:sxoffl', pstat_qstat, 0, junk,
     $     ps_sxoff)) call errquit('dft: dft_pstat_init', 0,0)
      if (.not. pstat_allocate('tce:vtoffl', pstat_qstat, 0, junk,
     $     ps_vtoff)) call errquit('dft: dft_pstat_init', 0,0)
c       
c
      end
      subroutine ccsd_t_pstat_print()
      implicit none
#include "pstat.fh"
#include "ccsd_t_ps.fh"
c     
      if (otceps) then
         call pstat_off(ps_tce)
         call pstat_print_all
         otceps = .false.
c
c     delete the TCE timers
c
         if (.not. pstat_free(ps_cctdbl))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_sdkrnl))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_sxoff))call errquit('dft_pstat?', 0,0)
         if (.not. pstat_free(ps_vtoff))call errquit('dft_pstat?', 0,0)
c
      endif
c
      end
c $Id$
      subroutine ccsd_t_fdist_init()
      implicit none
#include "dist.fh"
#include "ccsd_t_ps.fh"
#include "rtdb.fh"
      if (.not. otceps) return
c
      fdist_sdkrnl11 = getdist( 'sd_1_1_totl')
      fdist_sdcomp11 = getdist( 'sd_1_1_comp')
      fdist_sdkrnl12 = getdist( 'sd_1_2_totl')
      fdist_sdcomp12 = getdist( 'sd_1_2_comp')
      fdist_sdkrnl13 = getdist( 'sd_1_3_totl')
      fdist_sdcomp13 = getdist( 'sd_1_3_comp')
      fdist_sdkrnl14 = getdist( 'sd_1_4_totl')
      fdist_sdcomp14 = getdist( 'sd_1_4_comp')
      fdist_sdkrnl15 = getdist( 'sd_1_5_totl')
      fdist_sdcomp15 = getdist( 'sd_1_5_comp')
      fdist_sdkrnl16 = getdist( 'sd_1_6_totl')
      fdist_sdcomp16 = getdist( 'sd_1_6_comp')
      fdist_sdkrnl17 = getdist( 'sd_1_7_totl')
      fdist_sdcomp17 = getdist( 'sd_1_7_comp')
      fdist_sdkrnl18 = getdist( 'sd_1_8_totl')
      fdist_sdcomp18 = getdist( 'sd_1_8_comp')
      fdist_sdkrnl19 = getdist( 'sd_1_9_totl')
      fdist_sdcomp19 = getdist( 'sd_1_9_comp')
      fdist_sdkrnl21 = getdist( 'sd_2_1_totl')
      fdist_sdcomp21 = getdist( 'sd_2_1_comp')
      fdist_sdkrnl22 = getdist( 'sd_2_2_totl')
      fdist_sdcomp22 = getdist( 'sd_2_2_comp')
      fdist_sdkrnl23 = getdist( 'sd_2_3_totl')
      fdist_sdcomp23 = getdist( 'sd_2_3_comp')
      fdist_sdkrnl24 = getdist( 'sd_2_4_totl')
      fdist_sdcomp24 = getdist( 'sd_2_4_comp')
      fdist_sdkrnl25 = getdist( 'sd_2_5_totl')
      fdist_sdcomp25 = getdist( 'sd_2_5_comp')
      fdist_sdkrnl26 = getdist( 'sd_2_6_totl')
      fdist_sdcomp26 = getdist( 'sd_2_6_comp')
      fdist_sdkrnl27 = getdist( 'sd_2_7_totl')
      fdist_sdcomp27 = getdist( 'sd_2_7_comp')
      fdist_sdkrnl28 = getdist( 'sd_2_8_totl')
      fdist_sdcomp28 = getdist( 'sd_2_8_comp')
      fdist_sdkrnl29 = getdist( 'sd_2_9_totl')
      fdist_sdcomp29 = getdist( 'sd_2_9_comp')
      fdist_vtoff = getdist( 'v_and_t_offload')
      fdist_sxoff = getdist( 'triplesx_offload')
c
      end
c $Id$
      subroutine tce_sdtime0(t0)
      implicit none
#include "ccsd_t_ps.fh"
#include "util.fh"
      double precision t0
#ifdef USE_OPENMP
      double precision omp_get_wtime
      external omp_get_wtime
#endif
      if (otceps) then
#ifdef USE_OPENMP
         t0=omp_get_wtime()
#else
         t0=util_wallsec()
#endif
         call pstat_on(ps_sdkrnl)
      endif
      return
      end
c
      subroutine tce_sdtime(t0,deltat_mic,
     ,     fdist_sdkrnl_in,fdist_sdcomp_in)
      implicit none
#include "ccsd_t_ps.fh"
#include "util.fh"
      double precision t0
      integer deltat_mic
      integer fdist_sdkrnl_in,fdist_sdcomp_in
c
      integer deltat_tot
#ifdef USE_OPENMP
      double precision omp_get_wtime
      external omp_get_wtime
#endif
      if(.not.otceps) return
      call pstat_off(ps_sdkrnl)
#ifdef USE_OPENMP
      deltat_tot=(omp_get_wtime() - t0)*1.D6
#else
      deltat_tot=(util_wallsec() - t0)*1.D6
#endif
         call util_flush(6)
      if(deltat_tot.lt.deltat_mic) 
     D    deltat_tot=deltat_mic
      call updist(fdist_sdcomp_in, deltat_mic)
      call updist(fdist_sdkrnl_in, deltat_tot)
      call updist(fdist_vtoff, (deltat_tot-deltat_mic))
      return
      end
