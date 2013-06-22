      subroutine tce_mrcc_diis_tidy
c
c     Terminates the DIIS
c
      implicit none
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
      integer i,iref
c
      do iref=1,maxref
      do i = 1,maxdiism
        if (t1_diis_existm(i,iref))
     1    call deletefile(t1_diism(i,iref)) 
        if (r1_diis_existm(i,iref))
     1    call deletefile(r1_diism(i,iref)) 
        if (t2_diis_existm(i,iref))
     1    call deletefile(t2_diism(i,iref)) 
        if (r2_diis_existm(i,iref))
     1    call deletefile(r2_diism(i,iref)) 
        if (t3_diis_existm(i,iref))
     1    call deletefile(t3_diism(i,iref)) 
        if (r3_diis_existm(i,iref))
     1    call deletefile(r3_diism(i,iref)) 
        if (t4_diis_existm(i,iref))
     1    call deletefile(t4_diism(i,iref)) 
        if (r4_diis_existm(i,iref))
     1    call deletefile(r4_diism(i,iref)) 
      enddo
      enddo

      return
      end
c $Id$
