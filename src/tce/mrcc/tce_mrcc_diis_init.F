      subroutine tce_mrcc_diis_init
      implicit none
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
      integer i,iref

      do iref=1,maxref

      t1_pointerm(iref) = 0
      t2_pointerm(iref) = 0
      t3_pointerm(iref) = 0
      t4_pointerm(iref) = 0
      r1_pointerm(iref) = 0
      r2_pointerm(iref) = 0
      r3_pointerm(iref) = 0
      r4_pointerm(iref) = 0
      do i = 1,maxdiism
        t1_diis_existm(i,iref) = .false.
        r1_diis_existm(i,iref) = .false.
        t2_diis_existm(i,iref) = .false.
        r2_diis_existm(i,iref) = .false.
        t3_diis_existm(i,iref) = .false.
        r3_diis_existm(i,iref) = .false.
        t4_diis_existm(i,iref) = .false.
        r4_diis_existm(i,iref) = .false.
      enddo
      enddo

      return
      end
c $Id$
