      subroutine cis_xtidy
c
c $Id$
c
c     Delete files for trial and product vectors
c
      implicit none
#include "tce.fh"
#include "tce_diis.fh"
#include "tce_main.fh"
      integer i
c
      do i = 1,maxtrials
c        if(x1_ini_exist(i)) call deletefile(x1_ini(i))
c        if(x2_ini_exist(i)) call deletefile(x2_ini(i))
        if (x1_exist(i))  call deletefile(x1(i))
        if (x2_exist(i))  call deletefile(x2(i))
        if (x3_exist(i))  call deletefile(x3(i))
        if (x4_exist(i))  call deletefile(x4(i))
        if (xp1_exist(i)) call deletefile(xp1(i))
        if (xp2_exist(i)) call deletefile(xp2(i))
        if (xp3_exist(i)) call deletefile(xp3(i))
        if (xp4_exist(i)) call deletefile(xp4(i))
      enddo
      return
      end
