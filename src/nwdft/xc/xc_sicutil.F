      subroutine xc_sicinit(rtdb, test_sic, condfukui, exact_pot,
     &                      l_degen, i_degen, noc)
      implicit none
c     
c     $Id: xc_sicutil.F,v 1.3 2000-04-25 09:25:47 jolguin Exp $
c     
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
      integer l_degen, i_degen(2), noc(2), tot_oc
      integer rtdb ! [in]
      logical oprint ! [in]
      integer test_sic,condfukui,exact_pot ! [out]
c
      integer me
c
      me=ga_nodeid()
      oprint =util_print('SIC information', print_high)
      if (.not. rtdb_get(rtdb, 'dft:test_sic', mt_int, 1,
     &     test_sic))test_sic = 0
      if (.not. rtdb_get(rtdb, 'dft:condfukui', mt_int, 1,
     &   condfukui)) condfukui = 0
      if (.not. rtdb_get(rtdb, 'dft:exact_pot', mt_int, 1,
     &   exact_pot)) exact_pot = 0
      if (me.eq.0.and.oprint)then
         if (test_sic.eq.1) then
            write(LuOut,*) ' SIC perturbative approximation is used'
         elseif (test_sic.eq.2) then
            write(LuOut,*) ' SIC+OEP and NO localization'
         elseif (test_sic.eq.4) then
            write(LuOut,*) ' SIC+OEP and localization in the SCF'
         else
            write(LuOut,*) ' unidentified SIC/OEP option '
         end if
      end if
      if (test_sic.ne.0) then
        tot_oc = noc(1) + noc(2)
        if (.not.MA_Push_Get(MT_Int, tot_oc, 'degen_homo',
     &     l_degen, i_degen(1)))
     &     call errquit('dft_scf: cannot allocate degen homo',0)
      end if
      return
      end
      subroutine xc_sicdeg(i_degen, evals, noc)
      implicit none
#include "mafdecls.fh"
#include "cdft.fh"
c
      integer i_degen(2), evals(2)
      integer noc(*) ! [in]
c
      integer isp, j
      double precision homo_test,start
      do isp=1,ipol
         if (isp.eq.2) i_degen(2) = i_degen(1) + noc(1)
         if (noc(isp).ne.0) then
            homo_test = dbl_mb(evals(isp) + noc(isp) - 1)
            do j = 1,noc(isp)
               start = dbl_mb(evals(isp) + j - 1)
               if (dabs(start-homo_test).le.1.0d-04) then
                  Int_MB(i_degen(isp) + j - 1) = 1
               else
                  Int_MB(i_degen(isp) + j -1) = 0
               end if
            end do
         end if
      end do
      return
      end
