      SUBROUTINE tce_dens_pp_offset(l_a_offset,k_a_offset,size)
c     written by Peng-Dong Fan on 07/25/2007
c     dens_pp ( p5 p6 )
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
      INTEGER p5b
      INTEGER p6b
      length = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+p6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+p6b-1)) .eq. irrep_t) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &).ne.4)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'denspp',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('tce_dens_pp_offset',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+p6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+p6b-1)) .eq. irrep_t) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &).ne.4)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = (p6b-noab-1) + nvab * (p5b - noab - 1)
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      RETURN
      END
c $Id$
