      subroutine xc_sicinit(rtdb,test_sic,condfukui,exact_pot)
      implicit none
c     
c     $Id: xc_sicutil.F,v 1.1 2000-04-11 21:34:45 edo Exp $
c     
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
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
         elseif (test_sic.eq.3) then
            write(LuOut,*) ' SIC+OEP and localization after SCF'
         elseif (test_sic.eq.4) then
            write(LuOut,*) ' SIC+OEP and localization in the SCF'
         else
            write(LuOut,*) ' unidentified SIC/OEP option '
         end if
      end if
      return
      end
      subroutine xc_sicdeg(l_degen,i_degen,evals,noc)
      implicit none
#include "mafdecls.fh"
#include "cdft.fh"
c
      integer l_degen,i_degen(2)
      integer noc(*) ! [in]
      double precision evals(*) ! [in]
c
      integer isp,j,tot_oc
      double precision homo_test,start
c
      tot_oc = noc(1) + noc(2)
      if (.not.MA_Push_Get(MT_Int, tot_oc, 'degen_homo',
     &     l_degen, i_degen(1)))
     &     call errquit('dft_scf: cannot allocate degen homo',0)
      do isp=1,ipol
         if (isp.eq.2) i_degen(2) = i_degen(1) + noc(1)
         if (noc(isp).ne.0) then
            homo_test = evals(isp+noc(isp))
            do j = 1,noc(isp)
               start = evals(isp+j)
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
