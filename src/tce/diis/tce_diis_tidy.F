      subroutine tce_diis_tidy
c
c $Id$
c
c     Terminates the DIIS
c
      implicit none
#include "tce_diis.fh"
      integer i
c
      do i = 1,maxdiis
        if (t1_diis_exist(i))
     1    call deletefile(t1_diis(i)) 
        if (r1_diis_exist(i))
     1    call deletefile(r1_diis(i)) 
        if (t2_diis_exist(i))
     1    call deletefile(t2_diis(i)) 
        if (r2_diis_exist(i))
     1    call deletefile(r2_diis(i)) 
        if (t3_diis_exist(i))
     1    call deletefile(t3_diis(i)) 
        if (r3_diis_exist(i))
     1    call deletefile(r3_diis(i)) 
        if (t4_diis_exist(i))
     1    call deletefile(t4_diis(i)) 
        if (r4_diis_exist(i))
     1    call deletefile(r4_diis(i)) 
      enddo
      return
      end
