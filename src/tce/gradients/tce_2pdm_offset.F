      SUBROUTINE tce_2pdm_offset(l_a_offset,k_a_offset,size)
C     dens_hphp ( g1 g2 g3 g4 )
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
      INTEGER g1b
      INTEGER g2b
      INTEGER g3b
      INTEGER g4b
      length = 0
      DO g1b = 1,noab+nvab
      DO g2b = 1,noab+nvab
      DO g3b = 1,noab+nvab
      DO g4b = 1,noab+nvab
      IF (int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1) .eq. int_mb(k_spin+
     &g3b-1)+int_mb(k_spin+g4b-1)) THEN
      IF (ieor(int_mb(k_sym+g1b-1),ieor(int_mb(k_sym+g2b-1),ieor(int_mb
     &(k_sym+g3b-1),int_mb(k_sym+g4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-
     &1)+int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1).ne.8)) THEN
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
      DO g1b = 1,noab+nvab
      DO g2b = 1,noab+nvab
      DO g3b = 1,noab+nvab
      DO g4b = 1,noab+nvab
      IF (int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1) .eq. int_mb(k_spin+
     &g3b-1)+int_mb(k_spin+g4b-1)) THEN
      IF (ieor(int_mb(k_sym+g1b-1),ieor(int_mb(k_sym+g2b-1),ieor(int_mb
     &(k_sym+g3b-1),int_mb(k_sym+g4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-
     &1)+int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1).ne.8)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = g4b - 1 + (noab+nvab) * (g3b - 1 + 
     1 (noab+nvab) * (g2b - 1 + (noab+nvab) * (g1b - 1)))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+g1b-1) * int_mb(k_range+g2b-1) * int
     &_mb(k_range+g3b-1) * int_mb(k_range+g4b-1)
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
