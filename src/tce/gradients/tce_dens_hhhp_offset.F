      SUBROUTINE tce_dens_hhhp_offset(l_a_offset,k_a_offset,size)
C     dens_hhhp ( h1 h2 h3 p4 )
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER l_a_offset
      INTEGER k_a_offset
      INTEGER size
      INTEGER length
      INTEGER addr
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER p4b
      length = 0
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = 1,noab
      DO p4b = noab+1,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1) .eq. int_mb(k_spin
     &+h3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),ieor(int_m
     &b(k_sym+h3b-1),int_mb(k_sym+p4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b
     &-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+p4b-1).ne.8)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('ccsd_t2_2_2_2_1',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = 1,noab
      DO p4b = noab+1,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1) .eq. int_mb(k_spin
     &+h3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),ieor(int_m
     &b(k_sym+h3b-1),int_mb(k_sym+p4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b
     &-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+p4b-1).ne.8)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = p4b - noab - 1 + nvab * (h3b - 1 + noab
     &* (h2b - 1 + noab * (h1b - 1)))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * in
     &t_mb(k_range+h3b-1) * int_mb(k_range+p4b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
c $Id$
