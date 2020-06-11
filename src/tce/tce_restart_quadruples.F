        subroutine x4_restart_save(d_x4,k_x4_offset,size_x4,state,
     1                           hand_x4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_x4,k_x4_offset,size_x4
      INTEGER chunk_x4,hand_x4
      INTEGER l_aux,k_aux
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER h9b
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
        chunk_x4=((tile_dim)**8)
c
      if(state.eq.0) then
        call util_file_name('t4_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x4_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x4_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x4_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x4_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x4_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x4_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x4_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x4_8_copy',.false.,.false.,filename)
      else
        call errquit('x4_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x4_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_x4),
     1   dble(bytes)*dble(size_x4),dble(chunk_x4),hand_x4,
     2   X_SUFFIX)
     1   .ne.0)
     2   call errquit('x4_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      DO h4b = h3b,noab
      DO h9b = h4b,noab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h2b-1),ieo
     &r(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),int_mb(k_sym+h9b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+p5b-1) * int_mb(k_range+h2b-1) *
     1           int_mb(k_range+p6b-1) * int_mb(k_range+h3b-1) *
     1           int_mb(k_range+p7b-1) * int_mb(k_range+h4b-1) *
     2           int_mb(k_range+p8b-1) * int_mb(k_range+h9b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x4_restart_save: MA l_aux',0,MA_ERR)
!
!      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
!      &noab - 1 + nvab * (p1b - noab - 1)))
!
      key = h9b - 1 + noab * (h4b - 1 + noab * (h3b - 1 + noab *
     &    ( h2b - 1 + noab * (p8b - noab - 1 + nvab *
     &    ( p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab *
     &    ( p5b - noab - 1)))))))
!
      call tce_hash(int_mb(k_x4_offset),key,offset_aux)
      call ga_get(d_x4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_x4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x4_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('x4_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x4_restart_save: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end
c
c
      subroutine y4_restart_save(d_y4,k_y4_offset,size_y4,state,
     1                           hand_y4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_y4,k_y4_offset,size_y4
      INTEGER chunk_y4,hand_y4
      INTEGER l_aux,k_aux
      INTEGER h6b
      INTEGER h7b
      INTEGER h8b
      INTEGER h9b
      INTEGER p1b
      INTEGER p2b
      INTEGER p3b
      INTEGER p4b
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
        chunk_y4=((tile_dim)**8)
c
      if(state.eq.0) then
        call util_file_name('lambda4_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y4_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y4_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y4_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y4_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y4_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y4_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y4_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y4_8_copy',.false.,.false.,filename)
      else
        call errquit('y4_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y4_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_y4),
     1   dble(bytes)*dble(size_y4),dble(chunk_y4),hand_y4,
     2   Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y4_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO h9b = h8b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     &+int_mb(k_spin+h9b-1) .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-
     &1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+h9b-1),ieor(int_mb(k_sym+p1b-1),ieo
     &r(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),int_mb(k_sym+p4b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+h9b-1)+int_mb(k_spin+p1b-1)+i
     &nt_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+h7b-1) * int_mb(k_range+p2b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p3b-1)
     &         * int_mb(k_range+h9b-1) * int_mb(k_range+p4b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('y4_restart_save: MA l_aux',0,MA_ERR)
!
      key = p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     &    ( p2b - noab - 1 + nvab * (p1b - noab - 1 + nvab *
     &    ( h9b - 1 + noab * (h8b - 1 + noab * (h7b - 1 + noab *
     &    ( h6b - 1)))))))
!
      call tce_hash(int_mb(k_y4_offset),key,offset_aux)
      call ga_get(d_y4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_y4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y4_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('y4_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('y4_restart_save: MA l_aux',0,MA_ERR)
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
        subroutine tr4_restart_save(d_tr4,k_tr4_offset,size_tr4,axis,
     1                           hand_tr4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_tr4,k_tr4_offset,size_tr4
      INTEGER chunk_tr4,hand_tr4
      INTEGER l_aux,k_aux
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER h9b
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
        chunk_tr4=((tile_dim)**8)
c
      if(axis.eq.1) then
        call util_file_name('tr4_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr4_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr4_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr4_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr4_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr4_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr4_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr4_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr4_9_copy',.false.,.false.,filename)
      else
        call errquit('tr4_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr4_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_tr4),
     1   dble(bytes)*dble(size_tr4),dble(chunk_tr4),hand_tr4,
     2   TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr4_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      DO h4b = h3b,noab
      DO h9b = h4b,noab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h2b-1),ieo
     &r(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),int_mb(k_sym+h9b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+p5b-1) * int_mb(k_range+h2b-1) *
     1           int_mb(k_range+p6b-1) * int_mb(k_range+h3b-1) *
     1           int_mb(k_range+p7b-1) * int_mb(k_range+h4b-1) *
     2           int_mb(k_range+p8b-1) * int_mb(k_range+h9b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr4_restart_save: MA l_aux',0,MA_ERR)
!
!      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
!      &noab - 1 + nvab * (p1b - noab - 1)))
!
      key = h9b - 1 + noab * (h4b - 1 + noab * (h3b - 1 + noab *
     &    ( h2b - 1 + noab * (p8b - noab - 1 + nvab *
     &    ( p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab *
     &    ( p5b - noab - 1)))))))
!
      call tce_hash(int_mb(k_tr4_offset),key,offset_aux)
      call ga_get(d_tr4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_tr4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr4_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('tr4_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tr4_restart_save: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end
c
c
      subroutine yr4_restart_save(d_yr4,k_yr4_offset,size_yr4,axis,
     1                           hand_yr4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_yr4,k_yr4_offset,size_yr4
      INTEGER chunk_yr4,hand_yr4
      INTEGER l_aux,k_aux
      INTEGER h6b
      INTEGER h7b
      INTEGER h8b
      INTEGER h9b
      INTEGER p1b
      INTEGER p2b
      INTEGER p3b
      INTEGER p4b
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
        chunk_yr4=((tile_dim)**8)
c
      if(axis.eq.1) then
        call util_file_name('yr4_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr4_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr4_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr4_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr4_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr4_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr4_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr4_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr4_9_copy',.false.,.false.,filename)
      else
        call errquit('yr4_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr4_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_yr4),
     1   dble(bytes)*dble(size_yr4),dble(chunk_yr4),hand_yr4,
     2   YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr4_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO h9b = h8b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     &+int_mb(k_spin+h9b-1) .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-
     &1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+h9b-1),ieor(int_mb(k_sym+p1b-1),ieo
     &r(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),int_mb(k_sym+p4b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+h9b-1)+int_mb(k_spin+p1b-1)+i
     &nt_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+h7b-1) * int_mb(k_range+p2b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p3b-1)
     &         * int_mb(k_range+h9b-1) * int_mb(k_range+p4b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('yr4_restart_save: MA l_aux',0,MA_ERR)
!
      key = p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     &    ( p2b - noab - 1 + nvab * (p1b - noab - 1 + nvab *
     &    ( h9b - 1 + noab * (h8b - 1 + noab * (h7b - 1 + noab *
     &    ( h6b - 1)))))))
!
      call tce_hash(int_mb(k_yr4_offset),key,offset_aux)
      call ga_get(d_yr4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(hand_yr4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr4_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('yr4_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('yr4_restart_save: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end













      subroutine x4_restart_read(d_x4,k_x4_offset,size_x4,state,
     1                           hand_x4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_x4,k_x4_offset,size_x4
      INTEGER chunk_x4,hand_x4
      INTEGER l_aux,k_aux
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER h9b
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
        chunk_x4=((tile_dim)**8)
      if(state.eq.0) then
        call util_file_name('t4_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x4_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x4_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x4_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x4_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x4_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x4_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x4_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x4_8_copy',.false.,.false.,filename)
      else
        call errquit('x4_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x4_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_x4),
     1   dble(bytes)*dble(size_x4),dble(chunk_x4),hand_x4,
     2   X_SUFFIX)
     1   .ne.0)
     2   call errquit('x4_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      DO h4b = h3b,noab
      DO h9b = h4b,noab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h2b-1),ieo
     &r(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),int_mb(k_sym+h9b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+p5b-1) * int_mb(k_range+h2b-1) *
     1           int_mb(k_range+p6b-1) * int_mb(k_range+h3b-1) *
     1           int_mb(k_range+p7b-1) * int_mb(k_range+h4b-1) *
     2           int_mb(k_range+p8b-1) * int_mb(k_range+h9b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x4_restart_read: MA l_aux',0,MA_ERR)
!
      key = h9b - 1 + noab * (h4b - 1 + noab * (h3b - 1 + noab *
     &    ( h2b - 1 + noab * (p8b - noab - 1 + nvab *
     &    ( p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab *
     &    ( p5b - noab - 1)))))))
!
      call tce_hash(int_mb(k_x4_offset),key,offset_aux)
        ierrcode1=sf_read(hand_x4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x4_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('x4_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_x4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x4_restart_read: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end



      subroutine y4_restart_read(d_y4,k_y4_offset,size_y4,state,
     1                           hand_y4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_y4,k_y4_offset,size_y4
      INTEGER chunk_y4,hand_y4
      INTEGER l_aux,k_aux
      INTEGER h6b
      INTEGER h7b
      INTEGER h8b
      INTEGER h9b
      INTEGER p1b
      INTEGER p2b
      INTEGER p3b
      INTEGER p4b
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
        chunk_y4=((tile_dim)**8)
      if(state.eq.0) then
        call util_file_name('lambda4_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y4_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y4_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y4_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y4_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y4_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y4_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y4_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y4_8_copy',.false.,.false.,filename)
      else
        call errquit('y4_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y4_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_y4),
     1   dble(bytes)*dble(size_y4),dble(chunk_y4),hand_y4,
     2   Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y4_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO h9b = h8b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     &+int_mb(k_spin+h9b-1) .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-
     &1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+h9b-1),ieor(int_mb(k_sym+p1b-1),ieo
     &r(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),int_mb(k_sym+p4b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+h9b-1)+int_mb(k_spin+p1b-1)+i
     &nt_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+h7b-1) * int_mb(k_range+p2b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p3b-1)
     &         * int_mb(k_range+h9b-1) * int_mb(k_range+p4b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('y4_restart_save: MA l_aux',0,MA_ERR)
!
      key = p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     &    ( p2b - noab - 1 + nvab * (p1b - noab - 1 + nvab *
     &    ( h9b - 1 + noab * (h8b - 1 + noab * (h7b - 1 + noab *
     &    ( h6b - 1)))))))
!
      call tce_hash(int_mb(k_y4_offset),key,offset_aux)
        ierrcode1=sf_read(hand_y4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y4_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('y4_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_y4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('y4_restart_read: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end

      subroutine tr4_restart_read(d_tr4,k_tr4_offset,size_tr4,axis,
     1                           hand_tr4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_tr4,k_tr4_offset,size_tr4
      INTEGER chunk_tr4,hand_tr4
      INTEGER l_aux,k_aux
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER h9b
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
        chunk_tr4=((tile_dim)**8)
      if(axis.eq.1) then
        call util_file_name('tr4_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr4_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr4_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr4_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr4_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr4_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr4_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr4_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr4_9_copy',.false.,.false.,filename)
      else
        call errquit('tr4_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr4_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_tr4),
     1   dble(bytes)*dble(size_tr4),dble(chunk_tr4),hand_tr4,
     2   TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr4_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      DO h4b = h3b,noab
      DO h9b = h4b,noab
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h2b-1),ieo
     &r(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),int_mb(k_sym+h9b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)+int_mb(k_spin+h9b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+p5b-1) * int_mb(k_range+h2b-1) *
     1           int_mb(k_range+p6b-1) * int_mb(k_range+h3b-1) *
     1           int_mb(k_range+p7b-1) * int_mb(k_range+h4b-1) *
     2           int_mb(k_range+p8b-1) * int_mb(k_range+h9b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr4_restart_read: MA l_aux',0,MA_ERR)
!
      key = h9b - 1 + noab * (h4b - 1 + noab * (h3b - 1 + noab *
     &    ( h2b - 1 + noab * (p8b - noab - 1 + nvab *
     &    ( p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab *
     &    ( p5b - noab - 1)))))))
!
      call tce_hash(int_mb(k_tr4_offset),key,offset_aux)
        ierrcode1=sf_read(hand_tr4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr4_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('tr4_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_tr4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tr4_restart_read: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end



      subroutine yr4_restart_read(d_yr4,k_yr4_offset,size_yr4,axis,
     1                           hand_yr4,irrep)
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
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      INTEGER d_yr4,k_yr4_offset,size_yr4
      INTEGER chunk_yr4,hand_yr4
      INTEGER l_aux,k_aux
      INTEGER h6b
      INTEGER h7b
      INTEGER h8b
      INTEGER h9b
      INTEGER p1b
      INTEGER p2b
      INTEGER p3b
      INTEGER p4b
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
        chunk_yr4=((tile_dim)**8)
      if(axis.eq.1) then
        call util_file_name('yr4_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr4_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr4_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr4_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr4_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr4_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr4_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr4_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr4_9_copy',.false.,.false.,filename)
      else
        call errquit('yr4_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr4_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dble(bytes)*dble(size_yr4),
     1   dble(bytes)*dble(size_yr4),dble(chunk_yr4),hand_yr4,
     2   YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr4_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO h9b = h8b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
!
      IF (next.eq.count) THEN
!
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)
     &+int_mb(k_spin+h9b-1) .eq. int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-
     &1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h8b-1),ieor(int_mb(k_sym+h9b-1),ieor(int_mb(k_sym+p1b-1),ieo
     &r(int_mb(k_sym+p2b-1),ieor(int_mb(k_sym+p3b-1),int_mb(k_sym+p4b-1)
     &))))))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h8b-1)+int_mb(k_spin+h9b-1)+int_mb(k_spin+p1b-1)+i
     &nt_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).ne.1
     &6)) THEN
!
      size_aux = int_mb(k_range+h6b-1) * int_mb(k_range+p1b-1)
     &         * int_mb(k_range+h7b-1) * int_mb(k_range+p2b-1)
     &         * int_mb(k_range+h8b-1) * int_mb(k_range+p3b-1)
     &         * int_mb(k_range+h9b-1) * int_mb(k_range+p4b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('yr4_restart_save: MA l_aux',0,MA_ERR)
!
      key = p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     &    ( p2b - noab - 1 + nvab * (p1b - noab - 1 + nvab *
     &    ( h9b - 1 + noab * (h8b - 1 + noab * (h7b - 1 + noab *
     &    ( h6b - 1)))))))
!
      call tce_hash(int_mb(k_yr4_offset),key,offset_aux)
        ierrcode1=sf_read(hand_yr4,
     2                   dble(bytes)*dble(offset_aux),
     2    dble(bytes)*dble(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr4_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('yr4_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_yr4,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('yr4_restart_read: MA l_aux',0,MA_ERR)
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
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return
      end
c $Id$
