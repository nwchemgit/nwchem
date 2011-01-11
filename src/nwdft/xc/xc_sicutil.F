      subroutine xc_sicinit(rtdb, test_sic, condfukui, exact_pot,
     &                      l_degen, i_degen, noc, act_levels)
      implicit none
#include "errquit.fh"
c     
c     $Id$
c     
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
      integer l_degen, i_degen(2), noc(2), tot_oc, act_levels
      integer rtdb ! [in]
      integer test_sic,condfukui,exact_pot ! [out]
c
c
      if (.not. rtdb_get(rtdb, 'dft:test_sic', mt_int, 1,
     &     test_sic))test_sic = 0
      if (.not. rtdb_get(rtdb, 'dft:condfukui', mt_int, 1,
     &   condfukui)) condfukui = 0
      if (.not. rtdb_get(rtdb, 'dft:exact_pot', mt_int, 1,
     &   exact_pot)) exact_pot = 0
      if (test_sic.ne.0) then
        if (.not. rtdb_get(rtdb, 'dft:act_levels', mt_int, 1,
     &      act_levels)) act_levels = 10000
        tot_oc = noc(1) + noc(2)
        if (.not.MA_Push_Get(MT_Int, tot_oc, 'degen_homo',
     &     l_degen, i_degen(1)))
     &     call errquit('dft_scf: cannot allocate degen homo',0, MA_ERR)
      end if
      return
      end
      subroutine xc_sicdeg(i_degen, evals, n_levels, act_levels)
      implicit none
#include "mafdecls.fh"
#include "cdft.fh"
c
      integer i_degen(2), evals(2), n_levels(2)
c
      integer isp, j, aux_levels, n_orbitals, act_levels
      double precision start
      do isp=1,ipol
         if (isp.eq.2) i_degen(2) = i_degen(1) + noc(1)
         if (noc(isp).ne.0) then
           aux_levels = 1
           n_orbitals = 1
           Int_MB(i_degen(isp) + aux_levels - 1) = n_orbitals
           do j = noc(isp),2,-1
             start = (dbl_mb(evals(isp) + j - 1) -
     &                dbl_mb(evals(isp) + j - 2))
             if (start.le.1.0d-04) then
               n_orbitals = n_orbitals + 1
             else
               n_orbitals = 1
               aux_levels = aux_levels + 1
             end if
             Int_MB(i_degen(isp) + aux_levels - 1) = n_orbitals
           end do
           n_levels(isp) = aux_levels
         else
           n_levels(isp) = 0
         end if
      end do
      do isp = 1, ipol
        if (n_levels(isp).ne.0) then
          if(n_levels(isp).gt.act_levels) n_levels(isp) = act_levels
        end if
      end do
      return
      end
