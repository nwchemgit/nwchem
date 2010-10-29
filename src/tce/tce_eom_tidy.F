      subroutine tce_eom_tidy
c
c $Id$
c
c     Delete files for trial and product vectors
c
      implicit none
#include "tce_diis.fh"
      integer i
c
      do i = 1,maxtrials
        if (xc1_exist(i))  call deletefile(xc1(i))
        if (xc2_exist(i))  call deletefile(xc2(i))
        if (xc3_exist(i))  call deletefile(xc3(i))
        if (xc4_exist(i))  call deletefile(xc4(i))
        if (yc1_exist(i))  call deletefile(yc1(i))
        if (yc2_exist(i))  call deletefile(yc2(i))
        if (yc3_exist(i))  call deletefile(yc3(i))
        if (yc4_exist(i))  call deletefile(yc4(i))
      enddo
      return
      end



      subroutine tce_eom_xtidy
c
c $Id$
c
c     Delete files for trial and product vectors
c
      implicit none
#include "tce_diis.fh"
      integer i
c
      do i = 1,maxtrials
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



      subroutine tce_eom_ytidy
c
c $Id$
c
c     Delete files for trial and product vectors
c
      implicit none
#include "tce_diis.fh"
      integer i
c
      do i = 1,maxtrials
        if (y1_exist(i))  call deletefile(y1(i))
        if (y2_exist(i))  call deletefile(y2(i))
        if (y3_exist(i))  call deletefile(y3(i))
        if (y4_exist(i))  call deletefile(y4(i))
        if (yp1_exist(i)) call deletefile(yp1(i))
        if (yp2_exist(i)) call deletefile(yp2(i))
        if (yp3_exist(i)) call deletefile(yp3(i))
        if (yp4_exist(i)) call deletefile(yp4(i))
      enddo
      return
      end
