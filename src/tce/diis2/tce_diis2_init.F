      subroutine tce_diis2_init
c
c $Id$
c
c     Initializes file handle holders
c
      implicit none
#include "tce_diis.fh"
      integer i
c
      t1_pointer = 0
      t2_pointer = 0
      t3_pointer = 0
      t4_pointer = 0
      r1_pointer = 0
      r2_pointer = 0
      r3_pointer = 0
      r4_pointer = 0
      do i = 1,maxdiis
        t1_diis_exist(i) = .false.
        r1_diis_exist(i) = .false.
        t2_diis_exist(i) = .false.
        r2_diis_exist(i) = .false.
        t3_diis_exist(i) = .false.
        r3_diis_exist(i) = .false.
        t4_diis_exist(i) = .false.
        r4_diis_exist(i) = .false.
      enddo
      return
      end
