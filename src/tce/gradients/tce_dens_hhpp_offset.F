      SUBROUTINE tce_dens_hhpp_offset(l_a_offset,k_a_offset,size)
c     written by Peng-Dong Fan
C     dens_hhpp (h1 h2 p3 p4 )
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
      INTEGER p3b
      INTEGER p4b
      length = 0
      do h1b = 1,noab
      do h2b = h1b,noab
      DO p3b = 1+noab,noab+nvab
      DO p4b = p3b,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1) .eq. int_mb(k_spin+p
     &3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),ieor(int_mb(
     &k_sym+p3b-1),int_mb(k_sym+p4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1
     &)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).ne.8)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('tce_dens_hhpp_offset',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO p3b = 1+noab,noab+nvab
      DO p4b = p3b,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1) .eq. int_mb(k_spin+p
     &3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),ieor(int_mb(
     &k_sym+p3b-1),int_mb(k_sym+p4b-1)))) .eq. irrep_e) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1
     &)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).ne.8)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = p4b - noab - 1 + nvab * (p3b - noab - 1 
     1  + nvab * (h2b- 1 + noab * (h1b - 1)))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * int_
     &mb(k_range+p3b-1) * int_mb(k_range+p4b-1)
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
