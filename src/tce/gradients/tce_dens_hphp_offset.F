      SUBROUTINE tce_dens_hphp_offset(l_a_offset,k_a_offset,size)
C     dens_hphp ( h10 p3 h1 p5 )
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
      INTEGER p2b
      INTEGER h3b
      INTEGER p4b
      length = 0
      DO h1b = 1,noab
      DO p2b = noab+1,noab+nvab
      DO h3b = 1,noab
      DO p4b = noab+1,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+
     &h3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb
     &(k_sym+h3b-1),int_mb(k_sym+p4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+p2b-
     &1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+p4b-1).ne.8)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('tce_dens_hphp_offset',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO h1b = 1,noab
      DO p2b = noab+1,noab+nvab
      DO h3b = 1,noab
      DO p4b = noab+1,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+
     &h3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb
     &(k_sym+h3b-1),int_mb(k_sym+p4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+p2b-
     &1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+p4b-1).ne.8)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = p4b - noab - 1 + nvab * (h3b - 1 + noab
     &* (p2b - noab - 1 + nvab * (h1b - 1)))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+p2b-1) * int_mb(k_range+h1b-1) * int
     &_mb(k_range+h3b-1) * int_mb(k_range+p4b-1)
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
