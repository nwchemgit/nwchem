      subroutine x3_restart_save(d_x3,k_x3_offset,size_x3,state,
     1                           hand_x3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_x3,k_x3_offset,size_x3
      integer state,irrep
      INTEGER hand_x3
#if !defined(CCSDT)
      CALL ERRQUIT('CCSDT not compiled (tce_energy.F)',0,CAPMIS_ERR)
#else
      INTEGER chunk_x3
      INTEGER l_aux,k_aux
      INTEGER p2b,p3b,p4b,h1b,h5b,h6b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_x3=((tile_dim)**6)
c
      if(state.eq.0) then
        call util_file_name('t3_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x3_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x3_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x3_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x3_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x3_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x3_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x3_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x3_8_copy',.false.,.false.,filename)
      else
        call errquit('x3_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x3_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_x3),
     1   dble(bytes)*dble(size_x3),dble(chunk_x3),hand_x3,
     2   X_SUFFIX)
     1   .ne.0)
     2   call errquit('x3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+p4b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h5b-1),int
     &_mb(k_sym+h6b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1
     &)+int_mb(k_spin+p4b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+i
     &nt_mb(k_spin+h6b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) * 
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * 
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3_restart_save: MA l_aux',0,MA_ERR)
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab * 
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab * 
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_x3_offset),key,offset_aux)
      call ga_get(d_x3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_x3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x3_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x3_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x3_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'x3_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine y3_restart_save(d_y3,k_y3_offset,size_y3,state,
     1                           hand_y3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_y3,k_y3_offset,size_y3
      INTEGER chunk_y3,hand_y3
      INTEGER l_aux,k_aux
      INTEGER h6b,h7b,h8b,p1b,p2b,p3b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer state,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_y3=((tile_dim)**6)
      if(state.eq.0) then
        call util_file_name('lambda3_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y3_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y3_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y3_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y3_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y3_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y3_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y3_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y3_8_copy',.false.,.false.,filename)
      else
        call errquit('y3_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y3_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_y3),
     1   dble(bytes)*dble(size_y3),dble(chunk_y3),hand_y3,
     2   Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     & .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),int
     &_mb(k_sym+p3b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+i
     &nt_mb(k_spin+p3b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('y3_restart_save: MA l_aux',0,MA_ERR)
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_y3_offset),key,offset_aux)
      call ga_get(d_y3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_y3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y3_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('y3_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('y3_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'y3_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
      subroutine tr3_restart_save(d_tr3,k_tr3_offset,size_tr3,axis,
     1                           hand_tr3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_tr3,k_tr3_offset,size_tr3
      INTEGER chunk_tr3,hand_tr3
      INTEGER l_aux,k_aux
      INTEGER p2b,p3b,p4b,h1b,h5b,h6b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer axis,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_tr3=((tile_dim)**6)
c
      if(axis.eq.1) then
        call util_file_name('tr3_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr3_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr3_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr3_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr3_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr3_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr3_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr3_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr3_9_copy',.false.,.false.,filename)
      else
        call errquit('tr3_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr3_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_tr3),
     1   dble(bytes)*dble(size_tr3),dble(chunk_tr3),hand_tr3,
     2   TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+p4b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h5b-1),int
     &_mb(k_sym+h6b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1
     &)+int_mb(k_spin+p4b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+i
     &nt_mb(k_spin+h6b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) * 
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * 
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr3_restart_save: MA l_aux',0,MA_ERR)
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab * 
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab * 
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_tr3_offset),key,offset_aux)
      call ga_get(d_tr3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_tr3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr3_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('tr3_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tr3_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'tr3_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine yr3_restart_save(d_yr3,k_yr3_offset,size_yr3,axis,
     1                           hand_yr3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_yr3,k_yr3_offset,size_yr3
      INTEGER chunk_yr3,hand_yr3
      INTEGER l_aux,k_aux
      INTEGER h6b,h7b,h8b,p1b,p2b,p3b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer axis,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_yr3=((tile_dim)**6)
      if(axis.eq.1) then
        call util_file_name('yr3_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr3_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr3_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr3_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr3_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr3_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr3_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr3_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr3_9_copy',.false.,.false.,filename)
      else
        call errquit('yr3_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr3_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_yr3),
     1   dble(bytes)*dble(size_yr3),dble(chunk_yr3),hand_yr3,
     2   YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     & .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),int
     &_mb(k_sym+p3b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+i
     &nt_mb(k_spin+p3b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('yr3_restart_save: MA l_aux',0,MA_ERR)
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_yr3_offset),key,offset_aux)
      call ga_get(d_yr3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_yr3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr3_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('yr3_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('yr3_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'yr3_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
      subroutine x3_restart_read(d_x3,k_x3_offset,size_x3,state,
     1                           hand_x3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_x3,k_x3_offset,size_x3
      INTEGER chunk_x3,hand_x3
      INTEGER l_aux,k_aux
      INTEGER p2b,p3b,p4b,h1b,h5b,h6b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer state,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_x3=((tile_dim)**6)
      if(state.eq.0) then
        call util_file_name('t3_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x3_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x3_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x3_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x3_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x3_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x3_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x3_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x3_8_copy',.false.,.false.,filename)
      else
        call errquit('x3_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x3_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_x3),
     1   dble(bytes)*dble(size_x3),dble(chunk_x3),hand_x3,
     2   X_SUFFIX)
     1   .ne.0)
     2   call errquit('x3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+p4b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h5b-1),int
     &_mb(k_sym+h6b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1
     &)+int_mb(k_spin+p4b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+i
     &nt_mb(k_spin+h6b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) *
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) *
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3_restart_read: MA l_aux',0,MA_ERR)
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab *
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_x3_offset),key,offset_aux)
        ierrcode1=sf_read(hand_x3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x3_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x3_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_x3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x3_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine y3_restart_read(d_y3,k_y3_offset,size_y3,state,
     1                           hand_y3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_y3,k_y3_offset,size_y3
      INTEGER chunk_y3,hand_y3
      INTEGER l_aux,k_aux
      INTEGER h6b,h7b,h8b,p1b,p2b,p3b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer state,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_y3=((tile_dim)**6)
      if(state.eq.0) then
        call util_file_name('lambda3_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y3_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y3_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y3_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y3_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y3_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y3_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y3_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y3_8_copy',.false.,.false.,filename)
      else
        call errquit('y3_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y3_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_y3),
     1   dble(bytes)*dble(size_y3),dble(chunk_y3),hand_y3,
     2   Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     & .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),int
     &_mb(k_sym+p3b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+i
     &nt_mb(k_spin+p3b-1).ne.12)) THEN
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('y3_restart_read: MA l_aux',0,MA_ERR)
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_y3_offset),key,offset_aux)
        ierrcode1=sf_read(hand_y3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y3_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('y3_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_y3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('y3_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine tr3_restart_read(d_tr3,k_tr3_offset,size_tr3,axis,
     1                           hand_tr3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_tr3,k_tr3_offset,size_tr3
      INTEGER chunk_tr3,hand_tr3
      INTEGER l_aux,k_aux
      INTEGER p2b,p3b,p4b,h1b,h5b,h6b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer axis,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_tr3=((tile_dim)**6)
      if(axis.eq.1) then
        call util_file_name('tr3_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr3_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr3_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr3_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr3_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr3_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr3_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr3_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr3_9_copy',.false.,.false.,filename)
      else
        call errquit('tr3_restart_read: too many read',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr3_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_tr3),
     1   dble(bytes)*dble(size_tr3),dble(chunk_tr3),hand_tr3,
     2   TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+p4b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h5b-1),int
     &_mb(k_sym+h6b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1
     &)+int_mb(k_spin+p4b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+i
     &nt_mb(k_spin+h6b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) *
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) *
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr3_restart_read: MA l_aux',0,MA_ERR)
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab *
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_tr3_offset),key,offset_aux)
        ierrcode1=sf_read(hand_tr3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr3_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('tr3_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_tr3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tr3_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine yr3_restart_read(d_yr3,k_yr3_offset,size_yr3,axis,
     1                           hand_yr3,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_yr3,k_yr3_offset,size_yr3
      INTEGER chunk_yr3,hand_yr3
      INTEGER l_aux,k_aux
      INTEGER h6b,h7b,h8b,p1b,p2b,p3b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer axis,irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_yr3=((tile_dim)**6)
      if(axis.eq.1) then
        call util_file_name('yr3_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr3_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr3_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr3_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr3_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr3_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr3_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr3_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr3_9_copy',.false.,.false.,filename)
      else
        call errquit('yr3_restart_read: too many read',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr3_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_yr3),
     1   dble(bytes)*dble(size_yr3),dble(chunk_yr3),hand_yr3,
     2   YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     & .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),int
     &_mb(k_sym+p3b-1)))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+i
     &nt_mb(k_spin+p3b-1).ne.12)) THEN
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('yr3_restart_read: MA l_aux',0,MA_ERR)
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_yr3_offset),key,offset_aux)
        ierrcode1=sf_read(hand_yr3,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr3_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('yr3_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_yr3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('yr3_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine x3a_restart_save(d_x3a,k_x3a_offset,size_x3a,state,
     1                           hand_x3a,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_x3a,k_x3a_offset,size_x3a
      INTEGER chunk_x3a,hand_x3a
      INTEGER l_aux,k_aux
      INTEGER p5b,p6b,p7b,h3b,h4b,h9b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer state,irrep
      logical acolo
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_x3a=((tile_dim)**6)
c
      if(state.eq.0) then
        call util_file_name('t3a_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x3a_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x3a_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x3a_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x3a_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x3a_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x3a_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x3a_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x3a_8_copy',.false.,.false.,filename)
      else
        call errquit('x3a_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x3a_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_x3a),
     1   dble(bytes)*dble(size_x3a),dble(chunk_x3a),hand_x3a,
     2   X_SUFFIX)
     1   .ne.0)
     2   call errquit('x3a_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      DO h9b = h4b,noab
!
      IF (next.eq.count) THEN
!
      IF(acolo(p5b,p6b,p7b,h3b,h4b,h9b)) THEN 
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     & .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),int
     &_mb(k_sym+h9b-1)))))) .eq. irrep_x) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+i
     &nt_mb(k_spin+h9b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1)
     1         * int_mb(k_range+p7b-1) * int_mb(k_range+h3b-1)
     2         * int_mb(k_range+h4b-1) * int_mb(k_range+h9b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3a_restart_save: MA l_aux',0,MA_ERR)
!
      key = h9b - 1 + noab * (h4b - 1 + noab * (h3b - 1 + noab
     1    * (p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab
     2    * (p5b - noab - 1)))))
!
      call tce_hash(int_mb(k_x3a_offset),key,offset_aux)
      call ga_get(d_x3a,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_x3a,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x3a_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x3a_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x3a_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine x3a_restart_read(d_x3a,k_x3a_offset,size_x3a,state,
     1                           hand_x3a,irrep)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_restart.fh"
      INTEGER d_x3a,k_x3a_offset,size_x3a
      INTEGER chunk_x3a,hand_x3a
      INTEGER l_aux,k_aux
      INTEGER p5b,p6b,p7b,h3b,h4b,h9b
      INTEGER key
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer state,irrep
      logical acolo
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_x3a=((tile_dim)**6)
      if(state.eq.0) then
        call util_file_name('t3a_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x3a_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x3a_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x3a_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x3a_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x3a_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x3a_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x3a_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x3a_8_copy',.false.,.false.,filename)
      else
        call errquit('x3a_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x3a_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_x3a),
     1   dble(bytes)*dble(size_x3a),dble(chunk_x3a),hand_x3a,
     2   X_SUFFIX)
     1   .ne.0)
     2   call errquit('x3a_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      DO h9b = h4b,noab
!
      IF (next.eq.count) THEN
!
      IF(acolo(p5b,p6b,p7b,h3b,h4b,h9b)) THEN 
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     & .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),int
     &_mb(k_sym+h9b-1)))))) .eq. irrep_x) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+i
     &nt_mb(k_spin+h9b-1).ne.12)) THEN
!
      size_aux = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1)
     1         * int_mb(k_range+p7b-1) * int_mb(k_range+h3b-1)
     2         * int_mb(k_range+h4b-1) * int_mb(k_range+h9b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3a_restart_save: MA l_aux',0,MA_ERR)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3a_restart_save: MA l_aux',0,MA_ERR)
!
      key = h9b - 1 + noab * (h4b - 1 + noab * (h3b - 1 + noab
     1    * (p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab
     2    * (p5b - noab - 1)))))
!
      call tce_hash(int_mb(k_x3a_offset),key,offset_aux)
        ierrcode1=sf_read(hand_x3a,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x3a_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x3a_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_x3a,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x3a_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
#endif
      return 
      end 
c $Id$
