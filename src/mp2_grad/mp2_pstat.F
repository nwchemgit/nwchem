      subroutine mp2_pstat_init()
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "cmp2ps.fh"
c
      integer junk
c
      if (.not. pstat_allocate('mp2: moints', pstat_qstat, 0,junk,
     $     ps_moints)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: maket', pstat_qstat, 0,junk,
     $     ps_maket)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: pijab', pstat_qstat, 0,junk,
     $     ps_pijab)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: wijab', pstat_qstat, 0,junk,
     $     ps_wijab)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: lai', pstat_qstat, 0,junk,
     $     ps_lai)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: backt', pstat_qstat, 0,junk,
     $     ps_backt)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: nonsep', pstat_qstat, 0,junk,
     $     ps_nonsep)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: laifock', pstat_qstat, 0,junk,
     $     ps_laifock)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: cphf', pstat_qstat, 0,junk,
     $     ps_cphf)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: wijfock', pstat_qstat, 0,junk,
     $     ps_wijfock)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: sep', pstat_qstat, 0,junk,
     $     ps_sep)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('mp2: total', pstat_qstat, 0,junk,
     $     ps_total)) call errquit('mp2_pstat_init', 0, MEM_ERR)
c
      call pstat_on(ps_total)
c
      end
      subroutine mp2_pstat_tidy()
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "cmp2ps.fh"
#include "util.fh"
c
      call pstat_off(ps_total)
c
      if (util_print('profile', print_default)) 
     $     call pstat_print_all()
c
      if (.not. pstat_free(ps_moints)) call errquit('mp2_pstat_tidy',0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_maket)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_pijab)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_wijab)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_lai)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_backt)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_nonsep)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_laifock)) call errquit('mp2_pstat_tidy',0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_cphf)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_wijfock)) call errquit('mp2_pstat_tidy',0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_sep)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
      if (.not. pstat_free(ps_total)) call errquit('mp2_pstat_tidy', 0,
     &       MEM_ERR)
c
      end

      

