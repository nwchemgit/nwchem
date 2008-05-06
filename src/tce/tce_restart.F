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
        if(nodezero) write(6,*)'v2 restart write file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_2e),
     1   dfloat(bytes)*dfloat(size_2e),chunk_v2,handle_v2)
     1   .ne.0)
     2   call errquit('v2 restart: sf_problem',0,DISK_ERR)
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
      call ga_get(d_v2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_v2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('zones put: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('zones put: sf problem31-b',2,DISK_ERR)
        end if
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
        if(nodezero) write(6,*)'f1 restart write file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_1e),
     1   dfloat(bytes)*dfloat(size_1e),chunk_f1,handle_f1)
     1   .ne.0)
     2   call errquit('f1 restart: sf_problem',0,DISK_ERR)
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
     1  call errquit('tce_rest_f1_o: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_f1_offset_s),key,offset_aux)
      call ga_get(d_f1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_f1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('f1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('f1: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_f1_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
c
      subroutine d1_restart_save(d_f1,k_f1_offset_s,size_1e,axis,
     &                          handle_f1)
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
      INTEGER d_f1,k_f1_offset_s,size_1e
      INTEGER chunk_f1,handle_f1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_f1=((tile_dim)**2)
      if(axis.eq.1) 
     &   call util_file_name('d1_1_copy',.false.,.false.,filename)
      if(axis.eq.2) 
     &   call util_file_name('d1_2_copy',.false.,.false.,filename)
      if(axis.eq.3) 
     &   call util_file_name('d1_3_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'d1 restart write file:',filename,axis
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_1e),
     1   dfloat(bytes)*dfloat(size_1e),chunk_f1,handle_f1)
     1   .ne.0)
     2   call errquit('f1 restart: sf_problem',0,DISK_ERR)
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
     1   call errquit('tce_rest_f1_o: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_f1_offset_s),key,offset_aux)
      call ga_get(d_f1,offset_aux+1,offset_aux+size_aux,1,1,
     &     dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_f1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('f1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('f1: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_f1_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
c
      subroutine t1_restart_save(d_t1,k_t1_offset,size_t1,
     1                           handle_t1)
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
      INTEGER d_t1,k_t1_offset,size_t1
      INTEGER chunk_t1,handle_t1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t1=((tile_dim)**2)
        call util_file_name('t1_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t1 restart write file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t1),
     1   dfloat(bytes)*dfloat(size_t1),chunk_t1,handle_t1)
     1   .ne.0)
     2   call errquit('t1 restart: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep_t) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_t1_o: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_t1_offset),key,offset_aux)
      call ga_get(d_t1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_t1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('t1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('t1: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_t1_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine lrt1_restart_save(d_t1,k_t1_offset,size_t1,
     1                           handle_t1,axis)
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
      INTEGER d_t1,k_t1_offset,size_t1
      INTEGER chunk_t1,handle_t1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t1=((tile_dim)**2)
c
      if(axis.eq.1) 
     &   call util_file_name('t1r_1_copy',.false.,.false.,filename)
      if(axis.eq.2) 
     &   call util_file_name('t1r_2_copy',.false.,.false.,filename)
      if(axis.eq.3) 
     &   call util_file_name('t1r_3_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t1r restart write file:',filename,axis
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t1),
     1   dfloat(bytes)*dfloat(size_t1),chunk_t1,handle_t1)
     1   .ne.0)
     2   call errquit('t1r restart: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep_x) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_t1r_o: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_t1_offset),key,offset_aux)
      call ga_get(d_t1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_t1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('t1r: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('t1r: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_t1r_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine y1_restart_save(d_t1,k_t1_offset,size_t1,
     1                           handle_t1)
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
      INTEGER d_t1,k_t1_offset,size_t1
      INTEGER chunk_t1,handle_t1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t1=((tile_dim)**2)
        call util_file_name('y1_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'y1 restart write file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t1),
     1   dfloat(bytes)*dfloat(size_t1),chunk_t1,handle_t1)
     1   .ne.0)
     2   call errquit('y1 restart: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO p5b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),int_mb(k_sym+p5b-1)) .eq. irrep_y) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_y1_o: MA l_aux',0,MA_ERR)
      key=p5b - noab - 1 + nvab * (h6b - 1)
      call tce_hash(int_mb(k_t1_offset),key,offset_aux)
      call ga_get(d_t1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_t1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('f1: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_y1_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
c
      subroutine t2_restart_save(d_t2,k_t2_offset,size_t2,
     1                           handle_t2)
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
      INTEGER d_t2,k_t2_offset,size_t2
      INTEGER chunk_t2,handle_t2
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t2=((tile_dim)**4)
        call util_file_name('t2_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t2 restart write file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t2),
     1   dfloat(bytes)*dfloat(size_t2),chunk_t2,handle_t2)
     1   .ne.0)
     2   call errquit('t2 restart: sf_problem',0,DISK_ERR)
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
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_t2: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_t2_offset),key,offset_aux)
      call ga_get(d_t2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_t2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('copy t2: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('copy t2: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_t2: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine lrt2_restart_save(d_t2,k_t2_offset,size_t2,
     1                           handle_t2,axis)
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
      INTEGER d_t2,k_t2_offset,size_t2
      INTEGER chunk_t2,handle_t2
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t2=((tile_dim)**4)
        call util_file_name('t2r_copy',.false.,.false.,filename)
c
      if(axis.eq.1) 
     &   call util_file_name('t2r_1_copy',.false.,.false.,filename)
      if(axis.eq.2) 
     &   call util_file_name('t2r_2_copy',.false.,.false.,filename)
      if(axis.eq.3) 
     &   call util_file_name('t2r_3_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t2r restart write file:',filename,axis
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t2),
     1   dfloat(bytes)*dfloat(size_t2),chunk_t2,handle_t2)
     1   .ne.0)
     2   call errquit('t2r restart: sf_problem',0,DISK_ERR)
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
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep_x) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_t2r: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_t2_offset),key,offset_aux)
      call ga_get(d_t2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_t2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('copy t2r: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('copy t2r: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_t2r: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine y2_restart_save(d_t2,k_t2_offset,size_t2,
     1                           handle_t2)
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
      INTEGER d_t2,k_t2_offset,size_t2
      INTEGER chunk_t2,handle_t2
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t2=((tile_dim)**4)
        call util_file_name('y2_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'y2 restart write file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t2),
     1   dfloat(bytes)*dfloat(size_t2),chunk_t2,handle_t2)
     1   .ne.0)
     2   call errquit('y2 restart: sf_problem',0,DISK_ERR)
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
     &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep_y) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_y2: MA l_aux',0,MA_ERR)
      key=p2b - noab - 1 + nvab * (p1b - noab - 1
     &+ nvab * (h4b - 1 + noab * (h3b - 1)))
      call tce_hash(int_mb(k_t2_offset),key,offset_aux)
      call ga_get(d_t2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        ierrcode1=sf_write(handle_t2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('copy y2: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('copy y2: sf problem31-b',2,DISK_ERR)
        end if
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_y2: MA l_aux',0,MA_ERR)
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
        if(nodezero) write(6,*)'v2 restart read file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_2e),
     1   dfloat(bytes)*dfloat(size_2e),chunk_v2,handle_v2)
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
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('zones put: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
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
        if(nodezero) write(6,*)'f1 restart read file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_1e),
     1   dfloat(bytes)*dfloat(size_1e),chunk_f1,handle_f1)
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
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('f1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
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
      return 
      end 
c
c
c
      subroutine d1_restart_read(d_f1,k_f1_offset_s,size_1e,axis,
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
      INTEGER d_f1,k_f1_offset_s,size_1e
      INTEGER chunk_f1,handle_f1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_f1=((tile_dim)**2)
      if(axis.eq.1) 
     &   call util_file_name('d1_1_copy',.false.,.false.,filename)
      if(axis.eq.2) 
     &   call util_file_name('d1_2_copy',.false.,.false.,filename)
      if(axis.eq.3) 
     &   call util_file_name('d1_3_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'d1 restart read file:',filename,axis
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_1e),
     1   dfloat(bytes)*dfloat(size_1e),chunk_f1,handle_f1)
     1   .ne.0)
     2   call errquit('f1 restart: sf_problem',0,DISK_ERR)
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
     1   call errquit('tce_rest_f1: MA l_aux',0,MA_ERR)
      key=g5b - 1 + (noab+nvab) * (g6b - 1)
      call tce_hash_n(int_mb(k_f1_offset_s),key,offset_aux)
        ierrcode1=sf_read(handle_f1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('f1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('f1: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_f1,offset_aux+1,offset_aux+size_aux,1,1,
     &     dbl_mb(k_aux),1)
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
      subroutine t1_restart_read(d_t1,k_t1_offset,size_t1,
     1                           handle_t1)
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
      INTEGER d_t1,k_t1_offset,size_t1
      INTEGER chunk_t1,handle_t1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t1=((tile_dim)**2)
        call util_file_name('t1_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t1 restart read file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t1),
     1   dfloat(bytes)*dfloat(size_t1),chunk_t1,handle_t1)
     1   .ne.0)
     2   call errquit('t1 restart: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep_t) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_t1_o: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_t1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_t1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('t1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('t1: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_t1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_t1_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine lrt1_restart_read(d_t1,k_t1_offset,size_t1,
     1                           handle_t1,axis)
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
      INTEGER d_t1,k_t1_offset,size_t1
      INTEGER chunk_t1,handle_t1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t1=((tile_dim)**2)
c
      if(axis.eq.1) 
     &   call util_file_name('t1r_1_copy',.false.,.false.,filename)
      if(axis.eq.2) 
     &   call util_file_name('t1r_2_copy',.false.,.false.,filename)
      if(axis.eq.3) 
     &   call util_file_name('t1r_3_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t1r restart read file:',filename,axis
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t1),
     1   dfloat(bytes)*dfloat(size_t1),chunk_t1,handle_t1)
     1   .ne.0)
     2   call errquit('t1r restart: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep_x) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_t1r_o: MA l_aux',0,MA_ERR)
      key=h6b - 1 + noab * (p5b - noab - 1)
      call tce_hash(int_mb(k_t1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_t1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('t1r: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('t1r: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_t1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_t1r_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine y1_restart_read(d_t1,k_t1_offset,size_t1,
     1                           handle_t1)
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
      INTEGER d_t1,k_t1_offset,size_t1
      INTEGER chunk_t1,handle_t1
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t1=((tile_dim)**2)
        call util_file_name('y1_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'y1 restart read file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t1),
     1   dfloat(bytes)*dfloat(size_t1),chunk_t1,handle_t1)
     1   .ne.0)
     2   call errquit('y1 restart: sf_problem',0,DISK_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
      DO h6b = 1,noab
      DO p5b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),int_mb(k_sym+p5b-1)) .eq. irrep_y) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
      size_aux = int_mb(k_range+p5b-1)*int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1  l_aux,k_aux))
     1  call errquit('tce_rest_y1_o: MA l_aux',0,MA_ERR)
      key=p5b - noab - 1 + nvab * (h6b - 1)
      call tce_hash(int_mb(k_t1_offset),key,offset_aux)
        ierrcode1=sf_read(handle_t1,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('y1: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('f1: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_t1,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1   call errquit('tce_rest_y1_c: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
c
      subroutine t2_restart_read(d_t2,k_t2_offset,size_t2,
     1                           handle_t2)
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
      INTEGER d_t2,k_t2_offset,size_t2
      INTEGER chunk_t2,handle_t2
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t2=((tile_dim)**4)
        call util_file_name('t2_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t2 restart read file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t2),
     1   dfloat(bytes)*dfloat(size_t2),chunk_t2,handle_t2)
     1   .ne.0)
     2   call errquit('t2 restart: sf_problem',0,DISK_ERR)
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
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_t2: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_t2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_t2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('copy t2: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('copy t2: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_t2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_t2: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine lrt2_restart_read(d_t2,k_t2_offset,size_t2,
     1                           handle_t2,axis)
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
      INTEGER d_t2,k_t2_offset,size_t2
      INTEGER chunk_t2,handle_t2
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t2=((tile_dim)**4)
        call util_file_name('t2r_copy',.false.,.false.,filename)
c
      if(axis.eq.1) 
     &   call util_file_name('t2r_1_copy',.false.,.false.,filename)
      if(axis.eq.2) 
     &   call util_file_name('t2r_2_copy',.false.,.false.,filename)
      if(axis.eq.3) 
     &   call util_file_name('t2r_3_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'t2r restart read file:',filename,axis
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t2),
     1   dfloat(bytes)*dfloat(size_t2),chunk_t2,handle_t2)
     1   .ne.0)
     2   call errquit('t2r restart: sf_problem',0,DISK_ERR)
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
     &k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. irrep_x) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_t2r: MA l_aux',0,MA_ERR)
      key=h4b - 1 + noab * (h3b - 1 + noab * (p2b-
     &noab - 1 + nvab * (p1b - noab - 1)))
      call tce_hash(int_mb(k_t2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_t2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('copy t2r: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('copy t2r: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_t2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_t2r: MA l_aux',0,MA_ERR)
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
      return 
      end 
c
c
      subroutine y2_restart_read(d_t2,k_t2_offset,size_t2,
     1                           handle_t2)
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
      INTEGER d_t2,k_t2_offset,size_t2
      INTEGER chunk_t2,handle_t2
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
      EXTERNAL NXTASK
c
      nodezero=(ga_nodeid().eq.0)
c
        chunk_t2=((tile_dim)**4)
        call util_file_name('y2_copy',.false.,.false.,filename)
c
        if(nodezero) write(6,*)'y2 restart read file:',filename
        call util_flush(6)
c
        if (sf_create(filename,dfloat(bytes)*dfloat(size_t2),
     1   dfloat(bytes)*dfloat(size_t2),chunk_t2,handle_t2)
     1   .ne.0)
     2   call errquit('y2 restart: sf_problem',0,DISK_ERR)
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
     &k_sym+p1b-1),int_mb(k_sym+p2b-1)))) .eq. irrep_y) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) THEN
      size_aux = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_
     &mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
        if (.not.ma_push_get(mt_dbl,size_aux,'l_aux',
     1   l_aux,k_aux))
     1   call errquit('tce_rest_y2: MA l_aux',0,MA_ERR)
      key=p2b - noab - 1 + nvab * (p1b - noab - 1
     &+ nvab * (h4b - 1 + noab * (h3b - 1)))
      call tce_hash(int_mb(k_t2_offset),key,offset_aux)
        ierrcode1=sf_read(handle_t2,
     2                   dfloat(bytes)*dfloat(offset_aux),
     2    dfloat(bytes)*dfloat(size_aux),dbl_mb(k_aux),request)
        if(ierrcode1.ne.0) then
         write(6,*)'ierrcode1=',ierrcode1
         call util_flush(6)
         call errquit('copy y2: sf problem21-b',1,DISK_ERR)
        end if
        ierrcode2=sf_wait(request)
        if(ierrcode2.ne.0) then
         call errquit('copy y2: sf problem31-b',2,DISK_ERR)
        end if
      call ga_put(d_t2,offset_aux+1,offset_aux+size_aux,1,1,
     &            dbl_mb(k_aux),1)
        if (.not.ma_pop_stack(l_aux))
     1    call errquit('tce_rest_y2: MA l_aux',0,MA_ERR)
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
      return 
      end 
