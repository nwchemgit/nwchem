      subroutine tce_eom_init
c
c $Id: tce_eom_init.F,v 1.1 2003-07-25 17:56:56 sohirata Exp $
c
c     Initializes file handle holders
c
      implicit none
#include "tce.fh"
      integer i
c
      do i = 1,maxtrials
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
