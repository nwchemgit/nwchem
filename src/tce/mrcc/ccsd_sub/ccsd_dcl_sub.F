      SUBROUTINE xcl_t1t1_1(d_a,k_a_offset,d_b,k_b_offset,d_c,k_c_offse
     &t)
C     $Id$
C     i0 ( p3 p4 h1 h2 )_tt + = 1/2 * P( 4 ) * t ( p3 h1 )_t * t ( p4 h2 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER NXTASKsub
      INTEGER nxt
      INTEGER nprocs
      INTEGER count
      INTEGER p3b
      INTEGER p4b
      INTEGER h1b
      INTEGER h2b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p3b_1
      INTEGER h1b_1
      INTEGER p4b_2
      INTEGER h2b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      INTEGER l_c
      INTEGER k_c
      EXTERNAL NXTASKsub
      nprocs=GA_pgroup_NNODES(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      count = 0
      nxt=NXTASKsub(nprocs,1,int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      DO p3b = noab+1,noab+nvab
      DO p4b = noab+1,noab+nvab
      DO h1b = 1,noab
      DO h2b = 1,noab
      IF (nxt.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      IF (int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p3b-1),ieor(int_mb(k_sym+p4b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. ieor(irrep_t,irrep_t)) TH
     &EN
      dimc = int_mb(k_range+p3b-1) * int_mb(k_range+p4b-1) * int_mb(k_ra
     &nge+h1b-1) * int_mb(k_range+h2b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('dcl_t1t1_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p3b-1) .eq. int_mb(k_spin+h1b-1)) THEN
      IF (ieor(int_mb(k_sym+p3b-1),int_mb(k_sym+h1b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p3b,h1b,p3b_1,h1b_1)
      CALL TCE_RESTRICTED_2(p4b,h2b,p4b_2,h2b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p3b-1) * int_mb(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p4b-1) * int_mb(k_range+h2b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('dcl_t1t1_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &dcl_t1t1_1',2,MA_ERR)
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
      CALL GET_HASH_BLOCK_MA(dbl_mb(d_a),dbl_mb(k_a),dima,
     + int_mb(k_a_offset),(h1b_1 - 1 + noab * (p3b_1 - noab - 1)))
#else
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p3b_1 - noab - 1)))
#endif
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p3b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('dcl_t1t1_1',3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('dcl_t1t1_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &dcl_t1t1_1',5,MA_ERR)
#ifdef MRCC_LOCAL_T1
ckbn @did : localize t1amp
      CALL GET_HASH_BLOCK_MA(dbl_mb(d_b),dbl_mb(k_b),dimb,
     + int_mb(k_b_offset),(h2b_2 - 1 + noab * (p4b_2 - noab - 1)))
#else
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h2b_2
     & - 1 + noab * (p4b_2 - noab - 1)))
#endif
      CALL TCE_SORT_2(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p4b-1)
     &,int_mb(k_range+h2b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('dcl_t1t1_1',6,MA_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('dcl_t1t1_1',7,MA_E
     &RR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('dcl_t1t1_1',8,MA_E
     &RR)
      END IF
      END IF
      END IF
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &dcl_t1t1_1',9,MA_ERR)
      IF ((p3b .le. p4b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p4b-1),int_mb(k_range+h1b-1),int_mb(k_range+p3b-1)
     &,4,2,3,1,-mrccshift)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
     & 1 + noab * (h1b - 1 + noab * (p4b - noab - 1 + nvab * (p3b - noab
     & - 1)))))
      END IF
      IF ((p3b .le. p4b) .and. (h2b .le. h1b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p4b-1),int_mb(k_range+h1b-1),int_mb(k_range+p3b-1)
     &,4,2,1,3,mrccshift)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h2b - 1 + noab * (p4b - noab - 1 + nvab * (p3b - noab
     & - 1)))))
      END IF
c      IF ((p4b .le. p3b) .and. (h1b .le. h2b)) THEN
c      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
c     &,int_mb(k_range+p4b-1),int_mb(k_range+h1b-1),int_mb(k_range+p3b-1)
c     &,2,4,3,1,-1.0d0/1.0d0)
c      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
c     & 1 + noab * (h1b - 1 + noab * (p3b - noab - 1 + nvab * (p4b - noab
c     & - 1)))))
c      END IF
c      IF ((p4b .le. p3b) .and. (h2b .le. h1b)) THEN
c      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
c     &,int_mb(k_range+p4b-1),int_mb(k_range+h1b-1),int_mb(k_range+p3b-1)
c     &,2,4,1,3,1.0d0/1.0d0)
c      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
c     & 1 + noab * (h2b - 1 + noab * (p3b - noab - 1 + nvab * (p4b - noab
c     & - 1)))))
c      END IF
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('half_t1t1_1',10,MA_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('half_t1t1_1',11,MA_
     &ERR)
      END IF
      END IF
      END IF
      nxt=NXTASKsub(nprocs,1,int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      nxt=NXTASKsub(-nprocs,1,int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      call GA_Pgroup_SYNC(int_mb(k_innodes+ga_nnodes()+ga_nodeid()))
      RETURN
      END


