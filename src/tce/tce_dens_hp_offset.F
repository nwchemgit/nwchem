      SUBROUTINE tce_dens_hp_offset(l_a_offset,k_a_offset,size)
C     written by Peng-Dong Fan, 2007/07/25
C     dens_ph ( h5 p6 )
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
      INTEGER h5b
      INTEGER p6b
      length = 0
      do h5b = 1,noab
      do p6b = noab+1, noab+nvab
      IF (int_mb(k_spin+h5b-1) .eq. int_mb(k_spin+p6b-1)) THEN
      IF (ieor(int_mb(k_sym+h5b-1),int_mb(k_sym+p6b-1)) .eq. irrep_t) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+h5b-1)+int_mb(k_spin+p6b-1
     &).ne.4)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'denshp',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('tce_dens_ph_offset',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO h5b = 1,noab
      DO p6b = noab+1,noab+nvab
      IF (int_mb(k_spin+h5b-1) .eq. int_mb(k_spin+p6b-1)) THEN
      IF (ieor(int_mb(k_sym+h5b-1),int_mb(k_sym+p6b-1)) .eq. irrep_t) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+h5b-1)+int_mb(k_spin+p6b-1
     &).ne.4)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = (p6b-noab-1) + nvab * (h5b - 1)
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+h5b-1) * int_mb(k_range+p6b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      RETURN
      END
c $Id$
