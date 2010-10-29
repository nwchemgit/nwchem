      subroutine tce_eom_init
c
c $Id$
c
c     Initializes file handle holders
c
      implicit none
#include "tce_diis.fh"
      integer i
c
      do i = 1,maxtrials
        xc1_exist(i)  = .false.
        yc1_exist(i)  = .false.
        xc2_exist(i)  = .false.
        yc2_exist(i)  = .false.
        xc3_exist(i)  = .false.
        yc3_exist(i)  = .false.
        xc4_exist(i)  = .false.
        yc4_exist(i)  = .false.
        x1_exist(i)  = .false.
        y1_exist(i)  = .false.
        x2_exist(i)  = .false.
        y2_exist(i)  = .false.
        x3_exist(i)  = .false.
        y3_exist(i)  = .false.
        x4_exist(i)  = .false.
        y4_exist(i)  = .false.
        xp1_exist(i) = .false.
        yp1_exist(i) = .false.
        xp2_exist(i) = .false.
        yp2_exist(i) = .false.
        xp3_exist(i) = .false.
        yp3_exist(i) = .false.
        xp4_exist(i) = .false.
        yp4_exist(i) = .false.
      enddo
      return
      end
