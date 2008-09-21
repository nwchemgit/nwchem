c
c TO DO:
c
c 1. initialize all handles to 0 or -1 in tce_energy.F
c 2. in every restart_save call, add logic so that if a handle is 0 or -1 
c (whatever is chosen so that this value is always an invalid SF handle) then
c call sf_create, otherwise just do status check (what is this call in SF?) and write.
c
c motivation: handles are repeatedly created despite just writing over existing
c file and in response calculations this adds up to more than the number of
c valid handles (100) SF supports by default.
c
      subroutine v2_restart_save(d_v2,k_v2_alpha_offset_s,size_2e,
     1                           handle_v2)
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
      INTEGER d_v2,k_v2_alpha_offset_s,size_2e
      INTEGER chunk_v2,handle_v2
      INTEGER l_aux,k_aux
      INTEGER g3b,g4b,g1b,g2b
      INTEGER IROW,ICOL,IRES
      INTEGER INDEX_PAIR
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
        chunk_v2=((tile_dim)**4)
        call util_file_name('v2_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,1)'v2_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_2e),
     1   dfloat(bytes)*dfloat(size_2e),chunk_v2,handle_v2,V2_SUFFIX)
     1   .ne.0)
     2   call errquit('v2_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO g3b = 1,noa+nva   !j
      DO g4b = g3b,noa+nva !i
      DO g1b = 1,noa+nva   !l
      DO g2b = g1b,noa+nva !k
      IF (next.eq.count) THEN
      IF (int_mb(k_spin_alpha+g3b-1)+int_mb(k_spin_alpha+g4b-1).eq.
     &int_mb(k_spin_alpha+g1b-1)+int_mb(k_spin_alpha+g2b-1)) THEN
      IF (ieor(int_mb(k_sym_alpha+g3b-1),ieor(int_mb(k_sym_alpha+g4b-1),
     &ieor(int_mb(k_sym_alpha+g1b-1),int_mb(k_sym_alpha+g2b-1)))) .eq.
     &irrep_v) THEN
      IROW=INDEX_PAIR(g4b,g3b)
      ICOL=INDEX_PAIR(g2b,g1b)
      IF(IROW.GE.ICOL) THEN
      IRES=INDEX_PAIR(IROW,ICOL)
      size_aux = int_mb(k_range_alpha+g3b-1) *
     &int_mb(k_range_alpha+g4b-1) * int_mb(k_range_alpha+g1b-1) *
     &int_mb(k_range_alpha+g2b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('v2_restart_save: MA l_aux',0,MA_ERR)
      call tce_hash_n(int_mb(k_v2_alpha_offset_s),IRES,offset_aux)
      call ga_get(d_v2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_v2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(1))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('v2_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(1))) then
         call errquit('v2_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('v2_restart_save: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'v2_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
      subroutine f1_restart_save(d_f1,k_f1_offset_s,size_1e,
     1                           handle_f1)
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
      INTEGER d_f1,k_f1_offset_s,size_1e
      INTEGER chunk_f1,handle_f1
      INTEGER l_aux,k_aux
      INTEGER g5b,g6b
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
        chunk_f1=((tile_dim)**2)
        call util_file_name('f1_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,1)'f1_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_1e),
     1   dfloat(bytes)*dfloat(size_1e),chunk_f1,handle_f1,F1_SUFFIX)
     1   .ne.0)
     2   call errquit('f1_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO g6b = 1,noab+nvab
      DO g5b = 1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+g6b-1) .eq. int_mb(k_spin+g5b-1)) THEN
      IF (ieor(int_mb(k_sym+g6b-1),int_mb(k_sym+g5b-1)) .eq. irrep_f) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+g6b-1)+int_mb(k_spin+g5b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+g5b-1)*int_mb(k_range+g6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('f1_restart_save: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_f1_offset_s),key,offset_aux)
      call ga_get(d_f1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_f1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(1))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('f1_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(1))) then
         call errquit('f1_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('f1_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'f1_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
      subroutine d1_restart_save(d_d1,k_d1_offset_s,size_d1,axis,
     &                          handle_d1,irrep)
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
      INTEGER d_d1,k_d1_offset_s,size_d1
      INTEGER chunk_d1,handle_d1
      INTEGER l_aux,k_aux
      INTEGER g5b,g6b
      INTEGER key,axis
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_d1=((tile_dim)**2)
      if(axis.eq.1) then
        call util_file_name('d1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('d1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('d1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('q1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('q1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('q1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('q1_4_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('q1_5_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('q1_6_copy',.false.,.false.,filename)
      endif
c
        if(nodezero) write(6,1)'d1_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_d1),
     1   dfloat(bytes)*dfloat(size_d1),chunk_d1,handle_d1,D1_SUFFIX)
     1   .ne.0)
     2   call errquit('d1_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO g6b = 1,noab+nvab
      DO g5b = 1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+g6b-1) .eq. int_mb(k_spin+g5b-1)) THEN
      IF (ieor(int_mb(k_sym+g6b-1),int_mb(k_sym+g5b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+g6b-1)+int_mb(k_spin+g5b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+g5b-1)*int_mb(k_range+g6b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('d1_restart_save: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_d1_offset_s),key,offset_aux)
      call ga_get(d_d1,offset_aux+1,offset_aux+size_aux,1,1,
     &     dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_d1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(1))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('d1_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(1))) then
         call errquit('d1_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('d1_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'d1_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
      subroutine x1_restart_save(d_x1,k_x1_offset,size_x1,state,
     1                           handle_x1,irrep)
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
      INTEGER d_x1,k_x1_offset,size_x1
      INTEGER chunk_x1,handle_x1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
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
        chunk_x1=((tile_dim)**2)
      if(state.eq.0) then
        call util_file_name('t1_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x1_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x1_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x1_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x1_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x1_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x1_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x1_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x1_8_copy',.false.,.false.,filename)
      else
        call errquit('x1_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x1_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x1),
     1   dfloat(bytes)*dfloat(size_x1),chunk_x1,handle_x1,X_SUFFIX)
     1   .ne.0)
     2   call errquit('x1_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('x1_restart_save: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_x1_offset),key,offset_aux)
      call ga_get(d_x1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_x1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x1_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x1_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('x1_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'x1_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine tr1_restart_save(d_tr1,k_tr1_offset,size_tr1,axis,
     1                            handle_tr1,irrep)
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
      INTEGER d_tr1,k_tr1_offset,size_tr1
      INTEGER chunk_tr1,handle_tr1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
      INTEGER key
      INTEGER axis
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_tr1=((tile_dim)**2)
c
      if(axis.eq.1) then
        call util_file_name('tr1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr1_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr1_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr1_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr1_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr1_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr1_9_copy',.false.,.false.,filename)
      else
        call errquit('tr1_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr1_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_tr1),
     1   dfloat(bytes)*dfloat(size_tr1),chunk_tr1,handle_tr1,TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr1_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tr1_restart_save: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_tr1_offset),key,offset_aux)
      call ga_get(d_tr1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_tr1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr1_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('tr1_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tr1_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'tr1_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine y1_restart_save(d_y1,k_y1_offset,size_y1,state,
     1                           handle_y1,irrep)
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
      INTEGER d_y1,k_y1_offset,size_y1
      INTEGER chunk_y1,handle_y1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
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
        chunk_y1=((tile_dim)**2)
      if(state.eq.0) then
        call util_file_name('lambda1_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y1_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y1_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y1_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y1_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y1_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y1_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y1_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y1_8_copy',.false.,.false.,filename)
      else
        call errquit('y1_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y1_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_y1),
     1   dfloat(bytes)*dfloat(size_y1),chunk_y1,handle_y1,Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y1_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO p5b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),int_mb(k_sym+p5b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('y1_restart_save: MA l_aux',0,MA_ERR)
      key=p5b - noab - 1 + nvab * (h6b - 1)
      call tce_hash(int_mb(k_y1_offset),key,offset_aux)
      call ga_get(d_y1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_y1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y1_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('y1_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('y1_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'y1_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine yr1_restart_save(d_yr1,k_yr1_offset,size_yr1,axis,
     1                           handle_yr1,irrep)
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
      INTEGER d_yr1,k_yr1_offset,size_yr1
      INTEGER chunk_yr1,handle_yr1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
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
        chunk_yr1=((tile_dim)**2)
      if(axis.eq.1) then
        call util_file_name('yr1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr1_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr1_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr1_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr1_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr1_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr1_9_copy',.false.,.false.,filename)
      else
        call errquit('yr1_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr1_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_yr1),
     1   dfloat(bytes)*dfloat(size_yr1),chunk_yr1,handle_yr1,YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr1_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO p5b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),int_mb(k_sym+p5b-1)) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('yr1_restart_save: MA l_aux',0,MA_ERR)
      key=p5b - noab - 1 + nvab * (h6b - 1)
      call tce_hash(int_mb(k_yr1_offset),key,offset_aux)
      call ga_get(d_yr1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_yr1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr1_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('yr1_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('yr1_restart_save: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'yr1_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
      subroutine x2_restart_save(d_x2,k_x2_offset,size_x2,state,
     1                           handle_x2,irrep)
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
      INTEGER d_x2,k_x2_offset,size_x2
      INTEGER chunk_x2,handle_x2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
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
        chunk_x2=((tile_dim)**4)
      if(state.eq.0) then
        call util_file_name('t2_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x2_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x2_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x2_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x2_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x2_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x2_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x2_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x2_8_copy',.false.,.false.,filename)
      else
        call errquit('x2_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x2_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x2),
     1   dfloat(bytes)*dfloat(size_x2),chunk_x2,handle_x2,X_SUFFIX)
     1   .ne.0)
     2   call errquit('x2_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x2_restart_save: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_x2_offset),key,offset_aux)
      call ga_get(d_x2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_x2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x2_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x2_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x2_restart_save: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'x2_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine tr2_restart_save(d_tr2,k_tr2_offset,size_tr2,axis,
     1                            handle_tr2,irrep)
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
      INTEGER d_tr2,k_tr2_offset,size_tr2
      INTEGER chunk_tr2,handle_tr2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
      INTEGER key,axis
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_tr2=((tile_dim)**4)
      if(axis.eq.1) then
        call util_file_name('tr2_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr2_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr2_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr2_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr2_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr2_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr2_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr2_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr2_9_copy',.false.,.false.,filename)
      else
        call errquit('tr2_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr2_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_tr2),
     1   dfloat(bytes)*dfloat(size_tr2),chunk_tr2,handle_tr2,TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr2_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr2_restart_save: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_tr2_offset),key,offset_aux)
      call ga_get(d_tr2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_tr2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr2_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('tr2_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tr2_restart_save: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'tr2_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine y2_restart_save(d_y2,k_y2_offset,size_y2,state,
     1                           handle_y2,irrep)
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
      INTEGER d_y2,k_y2_offset,size_y2
      INTEGER chunk_y2,handle_y2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
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
        chunk_y2=((tile_dim)**4)
      if(state.eq.0) then
        call util_file_name('lambda2_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y2_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y2_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y2_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y2_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y2_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y2_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y2_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y2_8_copy',.false.,.false.,filename)
      else
        call errquit('y2_restart_save: too many saved',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y2_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_y2),
     1   dfloat(bytes)*dfloat(size_y2),chunk_y2,handle_y2,Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y2_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h3b = 1,noab
      DO h4b = h3b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
     &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('y2_restart_save: MA l_aux',0,MA_ERR)
      key=p2b - noab - 1 + nvab * (p1b - noab - 1
     &+ nvab * (h4b - 1 + noab * (h3b - 1)))
      call tce_hash(int_mb(k_y2_offset),key,offset_aux)
      call ga_get(d_y2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_y2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y2_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('y2_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('y2_restart_save: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'y2_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
      subroutine yr2_restart_save(d_yr2,k_yr2_offset,size_yr2,axis,
     1                           handle_yr2,irrep)
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
      INTEGER d_yr2,k_yr2_offset,size_yr2
      INTEGER chunk_yr2,handle_yr2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
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
        chunk_yr2=((tile_dim)**4)
      if(axis.eq.1) then
        call util_file_name('yr2_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr2_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr2_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr2_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr2_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr2_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr2_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr2_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr2_9_copy',.false.,.false.,filename)
      else
        call errquit('yr2_restart_save: too many saved',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr2_restart_save',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_yr2),
     1   dfloat(bytes)*dfloat(size_yr2),chunk_yr2,handle_yr2,YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr2_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h3b = 1,noab
      DO h4b = h3b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
     &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('yr2_restart_save: MA l_aux',0,MA_ERR)
      key=p2b - noab - 1 + nvab * (p1b - noab - 1
     &+ nvab * (h4b - 1 + noab * (h3b - 1)))
      call tce_hash(int_mb(k_yr2_offset),key,offset_aux)
      call ga_get(d_yr2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_yr2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr2_restart_save: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('yr2_restart_save: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('yr2_restart_save: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
        if(nodezero) write(6,2)'yr2_restart_save finished'
    2 format(1x,A80)
        call util_flush(6)
c
      return 
      end 
c
c
c
      subroutine x3_restart_save(d_x3,k_x3_offset,size_x3,state,
     1                           handle_x3,irrep)
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
      INTEGER d_x3,k_x3_offset,size_x3
      INTEGER chunk_x3,handle_x3
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x3),
     1   dfloat(bytes)*dfloat(size_x3),chunk_x3,handle_x3,X_SUFFIX)
     1   .ne.0)
     2   call errquit('x3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
!       DO p1b = noab+1,noab+nvab
!       DO p2b = p1b,noab+nvab
!       DO h3b = 1,noab
!       DO h4b = h3b,noab
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
!
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
!      &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
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
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) * 
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * 
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3_restart_save: MA l_aux',0,MA_ERR)
!
!      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
!      &noab - 1 + nvab * (p1b - noab - 1)))
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab * 
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab * 
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_x3_offset),key,offset_aux)
      call ga_get(d_x3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_x3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_y3,irrep)
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
      INTEGER d_y3,k_y3_offset,size_y3
      INTEGER chunk_y3,handle_y3
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_y3),
     1   dfloat(bytes)*dfloat(size_y3),chunk_y3,handle_y3,Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!       DO h3b = 1,noab
!       DO h4b = h3b,noab
!       DO p1b = noab+1,noab+nvab
!       DO p2b = p1b,noab+nvab
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
!
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
!      &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
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
     1   call errquit('y3_restart_save: MA l_aux',0,MA_ERR)
!
!       key=p2b - noab - 1 + nvab * (p1b - noab - 1
!      &+ nvab * (h4b - 1 + noab * (h3b - 1)))
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_y3_offset),key,offset_aux)
      call ga_get(d_y3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_y3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_tr3,irrep)
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
      INTEGER d_tr3,k_tr3_offset,size_tr3
      INTEGER chunk_tr3,handle_tr3
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_tr3),
     1   dfloat(bytes)*dfloat(size_tr3),chunk_tr3,handle_tr3,X_SUFFIX)
     1   .ne.0)
     2   call errquit('tr3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!
!       DO p1b = noab+1,noab+nvab
!       DO p2b = p1b,noab+nvab
!       DO h3b = 1,noab
!       DO h4b = h3b,noab
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
!
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
!      &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
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
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) * 
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * 
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr3_restart_save: MA l_aux',0,MA_ERR)
!
!      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
!      &noab - 1 + nvab * (p1b - noab - 1)))
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab * 
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab * 
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_tr3_offset),key,offset_aux)
      call ga_get(d_tr3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_tr3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_yr3,irrep)
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
      INTEGER d_yr3,k_yr3_offset,size_yr3
      INTEGER chunk_yr3,handle_yr3
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_yr3),
     1   dfloat(bytes)*dfloat(size_yr3),chunk_yr3,handle_yr3,Y_SUFFIX)
     1   .ne.0)
     2   call errquit('yr3_restart_save: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!       DO h3b = 1,noab
!       DO h4b = h3b,noab
!       DO p1b = noab+1,noab+nvab
!       DO p2b = p1b,noab+nvab
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
!
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
!      &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
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
     1   call errquit('yr3_restart_save: MA l_aux',0,MA_ERR)
!
!       key=p2b - noab - 1 + nvab * (p1b - noab - 1
!      &+ nvab * (h4b - 1 + noab * (h3b - 1)))
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_yr3_offset),key,offset_aux)
      call ga_get(d_yr3,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_yr3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
c
      subroutine v2_restart_read(d_v2,k_v2_alpha_offset_s,size_2e,
     1                           handle_v2)
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
      INTEGER d_v2,k_v2_alpha_offset_s,size_2e
      INTEGER chunk_v2,handle_v2
      INTEGER l_aux,k_aux
      INTEGER g3b,g4b,g1b,g2b
      INTEGER IROW,ICOL,IRES
      INTEGER INDEX_PAIR
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
        chunk_v2=((tile_dim)**4)
        call util_file_name('v2_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,1)'v2_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_2e),
     1   dfloat(bytes)*dfloat(size_2e),chunk_v2,handle_v2,V2_SUFFIX)
     1   .ne.0)
     2   call errquit('v2 restart read: sf_problem',0,DISK_ERR)
c
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO g3b = 1,noa+nva   !j
      DO g4b = g3b,noa+nva !i
      DO g1b = 1,noa+nva   !l
      DO g2b = g1b,noa+nva !k
      IF (next.eq.count) THEN
      IF (int_mb(k_spin_alpha+g3b-1)+int_mb(k_spin_alpha+g4b-1).eq.
     &int_mb(k_spin_alpha+g1b-1)+int_mb(k_spin_alpha+g2b-1)) THEN
      IF (ieor(int_mb(k_sym_alpha+g3b-1),ieor(int_mb(k_sym_alpha+g4b-1),
     &ieor(int_mb(k_sym_alpha+g1b-1),int_mb(k_sym_alpha+g2b-1)))) .eq.
     &irrep_v) THEN
      IROW=INDEX_PAIR(g4b,g3b)
      ICOL=INDEX_PAIR(g2b,g1b)
      IF(IROW.GE.ICOL) THEN
      IRES=INDEX_PAIR(IROW,ICOL)
      size_aux = int_mb(k_range_alpha+g3b-1) *
     &int_mb(k_range_alpha+g4b-1) * int_mb(k_range_alpha+g1b-1) *
     &int_mb(k_range_alpha+g2b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_v2: MA l_aux',0,MA_ERR)
      call tce_hash_n(int_mb(k_v2_alpha_offset_s),IRES,offset_aux)
        ierrcode1=sf_read(handle_v2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(1))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('zones put: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(1))) then
         call errquit('zones put: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_v2,offset_aux+1,offset_aux+size_aux,1,1,
     &       dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_v2: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine f1_restart_read(d_f1,k_f1_offset_s,size_1e,handle_f1)
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
      INTEGER d_f1,k_f1_offset_s,size_1e
      INTEGER chunk_f1,handle_f1
      INTEGER l_aux,k_aux
      INTEGER g5b,g6b
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
        chunk_f1=((tile_dim)**2)
        call util_file_name('f1_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,1)'f1_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_1e),
     1   dfloat(bytes)*dfloat(size_1e),chunk_f1,handle_f1,F1_SUFFIX)
     1   .ne.0)
     2   call errquit('f1 restart read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO g6b = 1,noab+nvab
      DO g5b = 1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+g6b-1) .eq. int_mb(k_spin+g5b-1)) THEN
      IF (ieor(int_mb(k_sym+g6b-1),int_mb(k_sym+g5b-1)) .eq. irrep_f) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+g6b-1)+int_mb(k_spin+g5b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+g5b-1)*int_mb(k_range+g6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_f1: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_f1_offset_s),key,offset_aux)
        ierrcode1=sf_read(handle_f1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(1))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('f1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(1))) then
         call errquit('f1: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_f1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_f1: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
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
      subroutine d1_restart_read(d_d1,k_d1_offset_s,size_d1,axis,
     1                           handle_d1,irrep)
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
      INTEGER d_d1,k_d1_offset_s,size_d1
      INTEGER chunk_d1,handle_d1
      INTEGER l_aux,k_aux
      INTEGER g5b,g6b
      INTEGER key,axis
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_d1=((tile_dim)**2)
      if(axis.eq.1) then
        call util_file_name('d1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('d1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('d1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('q1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('q1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('q1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('q1_4_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('q1_5_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('q1_6_copy',.false.,.false.,filename)
      endif
c
        if(nodezero) write(6,1)'d1_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_d1),
     1   dfloat(bytes)*dfloat(size_d1),chunk_d1,handle_d1,D1_SUFFIX)
     1   .ne.0)
     2   call errquit('d1_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO g6b = 1,noab+nvab
      DO g5b = 1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+g6b-1) .eq. int_mb(k_spin+g5b-1)) THEN
      IF (ieor(int_mb(k_sym+g6b-1),int_mb(k_sym+g5b-1)) .eq. irrep_d) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+g6b-1)+int_mb(k_spin+g5b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+g5b-1)*int_mb(k_range+g6b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('d1_restart_read: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_d1_offset_s),key,offset_aux)
        ierrcode1=sf_read(handle_d1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(1))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('d1_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(1))) then
         call errquit('d1_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_d1,offset_aux+1,offset_aux+size_aux,1,1,
     &     dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('d1_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
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
c
c
c
c
c
c T,Lambda,T_R - read
c
c
c
c
c
c
      subroutine x1_restart_read(d_x1,k_x1_offset,size_x1,state,
     1                           handle_x1,irrep)
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
      INTEGER d_x1,k_x1_offset,size_x1
      INTEGER chunk_x1,handle_x1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
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
        chunk_x1=((tile_dim)**2)
      if(state.eq.0) then
        call util_file_name('t1_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x1_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x1_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x1_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x1_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x1_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x1_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x1_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x1_8_copy',.false.,.false.,filename)
      else
        call errquit('x1_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x1_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x1),
     1   dfloat(bytes)*dfloat(size_x1),chunk_x1,handle_x1,X_SUFFIX)
     1   .ne.0)
     2   call errquit('x1_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('x1_restart_read: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_x1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_x1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x1_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x1_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_x1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('x1_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
      subroutine tr1_restart_read(d_tr1,k_tr1_offset,size_tr1,axis,
     1                           handle_tr1,irrep)
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
      INTEGER d_tr1,k_tr1_offset,size_tr1
      INTEGER chunk_tr1,handle_tr1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
      INTEGER key
      INTEGER axis
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_tr1=((tile_dim)**2)
c
      if(axis.eq.1) then
        call util_file_name('tr1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr1_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr1_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr1_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr1_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr1_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr1_9_copy',.false.,.false.,filename)
      else
        call errquit('tr1_restart_read: cannot exist',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr1_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_tr1),
     1   dfloat(bytes)*dfloat(size_tr1),chunk_tr1,handle_tr1,TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr1_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tr1_restart_read: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_tr1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_tr1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr1_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('tr1_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_tr1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tr1_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
      subroutine y1_restart_read(d_y1,k_y1_offset,size_y1,state,
     1                           handle_y1,irrep)
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
      INTEGER d_y1,k_y1_offset,size_y1
      INTEGER chunk_y1,handle_y1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
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
        chunk_y1=((tile_dim)**2)
      if(state.eq.0) then
        call util_file_name('lambda1_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y1_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y1_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y1_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y1_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y1_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y1_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y1_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y1_8_copy',.false.,.false.,filename)
      else
        call errquit('y1_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y1_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_y1),
     1   dfloat(bytes)*dfloat(size_y1),chunk_y1,handle_y1,Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y1_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO p5b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),int_mb(k_sym+p5b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('y1_restart_read: MA l_aux',0,MA_ERR)
      key=p5b - noab - 1 + nvab * (h6b - 1)
      call tce_hash(int_mb(k_y1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_y1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y1_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('y1_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_y1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('y1_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
      subroutine yr1_restart_read(d_yr1,k_yr1_offset,size_yr1,axis,
     1                           handle_yr1,irrep)
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
      INTEGER d_yr1,k_yr1_offset,size_yr1
      INTEGER chunk_yr1,handle_yr1
      INTEGER l_aux,k_aux
      INTEGER p5b,h6b
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
        chunk_yr1=((tile_dim)**2)
      if(axis.eq.1) then
        call util_file_name('yr1_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr1_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr1_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr1_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr1_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr1_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr1_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr1_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr1_9_copy',.false.,.false.,filename)
      else
        call errquit('yr1_restart_read: cannot exist',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr1_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_yr1),
     1   dfloat(bytes)*dfloat(size_yr1),chunk_yr1,handle_yr1,YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr1_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO p5b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),int_mb(k_sym+p5b-1)) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('yr1_restart_read: MA l_aux',0,MA_ERR)
      key=p5b - noab - 1 + nvab * (h6b - 1)
      call tce_hash(int_mb(k_yr1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_yr1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr1_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('yr1_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_yr1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('yr1_restart_read: MA l_aux',0,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
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
      subroutine x2_restart_read(d_x2,k_x2_offset,size_x2,state,
     1                           handle_x2,irrep)
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
      INTEGER d_x2,k_x2_offset,size_x2
      INTEGER chunk_x2,handle_x2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
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
        chunk_x2=((tile_dim)**4)
      if(state.eq.0) then
        call util_file_name('t2_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_x2_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_x2_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_x2_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_x2_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_x2_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_x2_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_x2_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_x2_8_copy',.false.,.false.,filename)
      else
        call errquit('x2_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'x2_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x2),
     1   dfloat(bytes)*dfloat(size_x2),chunk_x2,handle_x2,X_SUFFIX)
     1   .ne.0)
     2   call errquit('x2_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x2_restart_read: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_x2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_x2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('x2_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('x2_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_x2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('x2_restart_read: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
      subroutine tr2_restart_read(d_tr2,k_tr2_offset,size_tr2,axis,
     1                           handle_tr2,irrep)
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
      INTEGER d_tr2,k_tr2_offset,size_tr2
      INTEGER chunk_tr2,handle_tr2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
      INTEGER key,axis
      INTEGER size_aux,offset_aux
      INTEGER ierrcode1,ierrcode2
      INTEGER request
      logical nodezero
      character*255 filename
      INTEGER NXTASK
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer irrep
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_tr2=((tile_dim)**4)
      if(axis.eq.1) then
        call util_file_name('tr2_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('tr2_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('tr2_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('tr2_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('tr2_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('tr2_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('tr2_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('tr2_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('tr2_9_copy',.false.,.false.,filename)
      else
        call errquit('tr2_restart_read: cannot exist',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'tr2_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_tr2),
     1   dfloat(bytes)*dfloat(size_tr2),chunk_tr2,handle_tr2,TR_SUFFIX)
     1   .ne.0)
     2   call errquit('tr2_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr2_restart_read: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_tr2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_tr2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('tr2_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('tr2_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_tr2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tr2_restart_read: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
      subroutine y2_restart_read(d_y2,k_y2_offset,size_y2,state,
     1                           handle_y2,irrep)
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
      INTEGER d_y2,k_y2_offset,size_y2
      INTEGER chunk_y2,handle_y2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
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
        chunk_y2=((tile_dim)**4)
      if(state.eq.0) then
        call util_file_name('lambda2_copy',.false.,.false.,filename)
      elseif(state.eq.1) then
        call util_file_name('eom_y2_1_copy',.false.,.false.,filename)
      elseif(state.eq.2) then
        call util_file_name('eom_y2_2_copy',.false.,.false.,filename)
      elseif(state.eq.3) then
        call util_file_name('eom_y2_3_copy',.false.,.false.,filename)
      elseif(state.eq.4) then
        call util_file_name('eom_y2_4_copy',.false.,.false.,filename)
      elseif(state.eq.5) then
        call util_file_name('eom_y2_5_copy',.false.,.false.,filename)
      elseif(state.eq.6) then
        call util_file_name('eom_y2_6_copy',.false.,.false.,filename)
      elseif(state.eq.7) then
        call util_file_name('eom_y2_7_copy',.false.,.false.,filename)
      elseif(state.eq.8) then
        call util_file_name('eom_y2_8_copy',.false.,.false.,filename)
      else
        call errquit('y2_restart_read: cannot exist',state,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'y2_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_y2),
     1   dfloat(bytes)*dfloat(size_y2),chunk_y2,handle_y2,Y_SUFFIX)
     1   .ne.0)
     2   call errquit('y2_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h3b = 1,noab
      DO h4b = h3b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
     &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('y2_restart_read: MA l_aux',0,MA_ERR)
      key=p2b - noab - 1 + nvab * (p1b - noab - 1
     &+ nvab * (h4b - 1 + noab * (h3b - 1)))
      call tce_hash(int_mb(k_y2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_y2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(2))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y2_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(2))) then
         call errquit('y2_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_y2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('y2_restart_read: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
      subroutine yr2_restart_read(d_yr2,k_yr2_offset,size_yr2,axis,
     1                           handle_yr2,irrep)
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
      INTEGER d_yr2,k_yr2_offset,size_yr2
      INTEGER chunk_yr2,handle_yr2
      INTEGER l_aux,k_aux
      INTEGER p1b,p2b,h3b,h4b
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
        chunk_yr2=((tile_dim)**4)
      if(axis.eq.1) then
        call util_file_name('yr2_1_copy',.false.,.false.,filename)
      elseif(axis.eq.2) then
        call util_file_name('yr2_2_copy',.false.,.false.,filename)
      elseif(axis.eq.3) then
        call util_file_name('yr2_3_copy',.false.,.false.,filename)
      elseif(axis.eq.4) then
        call util_file_name('yr2_4_copy',.false.,.false.,filename)
      elseif(axis.eq.5) then
        call util_file_name('yr2_5_copy',.false.,.false.,filename)
      elseif(axis.eq.6) then
        call util_file_name('yr2_6_copy',.false.,.false.,filename)
      elseif(axis.eq.7) then
        call util_file_name('yr2_7_copy',.false.,.false.,filename)
      elseif(axis.eq.8) then
        call util_file_name('yr2_8_copy',.false.,.false.,filename)
      elseif(axis.eq.9) then
        call util_file_name('yr2_9_copy',.false.,.false.,filename)
      else
        call errquit('yr2_restart_read: cannot exist',axis,DISK_ERR)
      endif
c
        if(nodezero) write(6,1)'yr2_restart_read',filename
    1 format(1x,A20,' filename: ',A80)
        call util_flush(6)
c
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_yr2),
     1   dfloat(bytes)*dfloat(size_yr2),chunk_yr2,handle_yr2,YR_SUFFIX)
     1   .ne.0)
     2   call errquit('yr2_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h3b = 1,noab
      DO h4b = h3b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
     &3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
     &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('yr2_restart_read: MA l_aux',0,MA_ERR)
      key=p2b - noab - 1 + nvab * (p1b - noab - 1
     &+ nvab * (h4b - 1 + noab * (h3b - 1)))
      call tce_hash(int_mb(k_yr2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_yr2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if ((ierrcode1.ne.0).and.(.not.restart_ignore_failure(3))) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('yr2_restart_read: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if ((ierrcode2.ne.0).and.(.not.restart_ignore_failure(3))) then
         call errquit('yr2_restart_read: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_yr2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('yr2_restart_read: MA l_aux',0,MA_ERR)
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
      next = NXTASK(-nprocs, 1)
      call GA_SYNC()
c
      return 
      end 
c
c
c
      subroutine x3_restart_read(d_x3,k_x3_offset,size_x3,state,
     1                           handle_x3,irrep)
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
      INTEGER d_x3,k_x3_offset,size_x3
      INTEGER chunk_x3,handle_x3
      INTEGER l_aux,k_aux
c      INTEGER p1b,p2b,h3b,h4b
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x3),
     1   dfloat(bytes)*dfloat(size_x3),chunk_x3,handle_x3,X_SUFFIX)
     1   .ne.0)
     2   call errquit('x3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!      DO p1b = noab+1,noab+nvab
!      DO p2b = p1b,noab+nvab
!      DO h3b = 1,noab
!      DO h4b = h3b,noab
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
!      &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
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
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) *
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) *
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('x3_restart_read: MA l_aux',0,MA_ERR)
!
!      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
!      &noab - 1 + nvab * (p1b - noab - 1)))
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab *
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_x3_offset),key,offset_aux)
        ierrcode1=sf_read(handle_x3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_y3,irrep)
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
      INTEGER d_y3,k_y3_offset,size_y3
      INTEGER chunk_y3,handle_y3
      INTEGER l_aux,k_aux
c      INTEGER p1b,p2b,h3b,h4b
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
        call util_file_name('t3_copy',.false.,.false.,filename)
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_y3),
     1   dfloat(bytes)*dfloat(size_y3),chunk_y3,handle_y3,X_SUFFIX)
     1   .ne.0)
     2   call errquit('y3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!      DO p1b = noab+1,noab+nvab
!      DO p2b = p1b,noab+nvab
!      DO h3b = 1,noab
!      DO h4b = h3b,noab
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
!      &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
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
!       key=p2b - noab - 1 + nvab * (p1b - noab - 1
!      &+ nvab * (h4b - 1 + noab * (h3b - 1)))
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_y3_offset),key,offset_aux)
        ierrcode1=sf_read(handle_y3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_tr3,irrep)
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
      INTEGER d_tr3,k_tr3_offset,size_tr3
      INTEGER chunk_tr3,handle_tr3
      INTEGER l_aux,k_aux
c      INTEGER p1b,p2b,h3b,h4b
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_tr3),
     1   dfloat(bytes)*dfloat(size_tr3),chunk_tr3,handle_tr3,X_SUFFIX)
     1   .ne.0)
     2   call errquit('tr3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!      DO p1b = noab+1,noab+nvab
!      DO p2b = p1b,noab+nvab
!      DO h3b = 1,noab
!      DO h4b = h3b,noab
      DO p2b = noab+1,noab+nvab
      DO p3b = p2b,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h1b = 1,noab
      DO h5b = h1b,noab
      DO h6b = h5b,noab
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
!      &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
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
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
!
      size_aux = int_mb(k_range+p2b-1) * int_mb(k_range+p3b-1) *
     1           int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) *
     2           int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
!
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tr3_restart_read: MA l_aux',0,MA_ERR)
!
!      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
!      &noab - 1 + nvab * (p1b - noab - 1)))
!
      key = h6b - 1 + noab * (h5b - 1 + noab * (h1b - 1 + noab *
     1     (p4b - noab - 1 + nvab * (p3b - noab - 1 + nvab *
     2     (p2b - noab - 1)))))
!
      call tce_hash(int_mb(k_tr3_offset),key,offset_aux)
        ierrcode1=sf_read(handle_tr3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_yr3,irrep)
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
      INTEGER d_yr3,k_yr3_offset,size_yr3
      INTEGER chunk_yr3,handle_yr3
      INTEGER l_aux,k_aux
c      INTEGER p1b,p2b,h3b,h4b
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_yr3),
     1   dfloat(bytes)*dfloat(size_yr3),chunk_yr3,handle_yr3,X_SUFFIX)
     1   .ne.0)
     2   call errquit('yr3_restart_read: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
!      DO p1b = noab+1,noab+nvab
!      DO p2b = p1b,noab+nvab
!      DO h3b = 1,noab
!      DO h4b = h3b,noab
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h8b = h7b,noab
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      IF (next.eq.count) THEN
!
!       IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. int_mb(k_spin+h
!      &3b-1)+int_mb(k_spin+h4b-1)) THEN
!       IF (ieor(int_mb(k_sym+h3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(
!      &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep) THEN
!       IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
!      &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
!       size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
!      &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
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
!       key=p2b - noab - 1 + nvab * (p1b - noab - 1
!      &+ nvab * (h4b - 1 + noab * (h3b - 1)))
!
      key = p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * 
     1     (p1b - noab - 1 + nvab * (h8b - 1 + noab *
     2     (h7b - 1 + noab * (h6b - 1)))))
!
      call tce_hash(int_mb(k_yr3_offset),key,offset_aux)
        ierrcode1=sf_read(handle_yr3,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_x3a,irrep)
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
      INTEGER d_x3a,k_x3a_offset,size_x3a
      INTEGER chunk_x3a,handle_x3a
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x3a),
     1   dfloat(bytes)*dfloat(size_x3a),chunk_x3a,handle_x3a,X_SUFFIX)
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
        ierrcode1=sf_write(handle_x3a,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
     1                           handle_x3a,irrep)
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
      INTEGER d_x3a,k_x3a_offset,size_x3a
      INTEGER chunk_x3a,handle_x3a
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
        if (sf_create_suffix(filename,dfloat(bytes)*dfloat(size_x3a),
     1   dfloat(bytes)*dfloat(size_x3a),chunk_x3a,handle_x3a,X_SUFFIX)
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
        ierrcode1=sf_read(handle_x3a,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
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
      return 
      end 
