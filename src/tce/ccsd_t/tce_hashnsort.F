      logical function tce_hashnsort(dryrun,
     d     d_a, d_b,
     I     k_a_offset,k_b_offset,
     I     p4b,p5b,h1b,h7b,
     I     p6b,h2b,h3b,
     O     t2sub,v2sub
     P     )
      implicit none
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "offl.fh"
      integer k_a_offset ! [in]
      integer k_b_offset ! [in]
      logical dryrun ! [in]
      integer     p4b,p5b,h1b,h7b,
     I     p6b,h2b,h3b ! [in]
      double precision d_a(*),d_b(*) ! [in]
      double precision t2sub(*),v2sub(*) ! [out]
c
      integer p4b_1,p5b_1,h1b_1,h7b_1
      integer p6b_2,h7b_2,h2b_2,h3b_2
      integer l_a,k_a
      integer dima,dimb,dima_sort,dimb_sort,dim_common
c

      IF ((int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1).eq.
     Q     int_mb(k_spin+h1b-1)+int_mb(k_spin+h7b-1)).and.
     A     (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p5b-1),
     I     ieor(int_mb(k_sym+h1b-1),int_mb(k_sym+h7b-1)))) .eq. 
     I     irrep_t)) THEN
         dim_common = int_mb(k_range+h7b-1)
         dima_sort = int_mb(k_range+p4b-1) * int_mb(k_range+p5b-1) * 
     *        int_mb(k_range+h1b-1)
         dima = dim_common * dima_sort
         dimb_sort = int_mb(k_range+p6b-1) * int_mb(k_range+h2b-1) * 
     *        int_mb(k_range+h3b-1)
         dimb = dim_common * dimb_sort
         IF ((.not.dryrun).and.(dima.gt.0) .and. (dimb.gt.0)) THEN
            CALL TCE_RESTRICTED_4(p4b,p5b,h1b,h7b,
     P           p4b_1,p5b_1,h1b_1,h7b_1)
            CALL TCE_RESTRICTED_4(p6b,h7b,h2b,h3b,
     P           p6b_2,h7b_2,h2b_2,h3b_2)
            IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL 
     E           ERRQUIT('ccsd_t_doubles_1',2,MA_ERR)
            IF ((h7b .lt. h1b)) THEN
               CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,
     I              int_mb(k_a_offset),(h1b_1 - 1 + noab * (h7b_1 - 1 + 
     N              noab * (p5b_1 - noab - 1 + nvab * 
     *              (p4b_1 - noab - 1)))))
               CALL TCE_SORT_4(dbl_mb(k_a),t2sub,int_mb(k_range+p4b-1),
     &              int_mb(k_range+p5b-1),int_mb(k_range+h7b-1),
     I              int_mb(k_range+h1b-1),4,2,1,3,-1.0d0)
            END IF
            IF ((h1b .le. h7b)) THEN
               CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,
     I              int_mb(k_a_offset),(h7b_1 - 1 + noab * (h1b_1 - 1 +
     N              noab * (p5b_1 - noab - 1 + nvab * (p4b_1 - 
     N              noab - 1)))))
               CALL TCE_SORT_4(dbl_mb(k_a),t2sub,int_mb(k_range+p4b-1),
     &              int_mb(k_range+p5b-1),int_mb(k_range+h1b-1),
     I              int_mb(k_range+h7b-1),3,2,1,4,1.0d0)
            END IF
            IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsd_t_doubles_1',
     3           3,MA_ERR)
            IF ((h7b .le. p6b)) THEN
               if(.not.intorb) then
                  CALL GET_HASH_BLOCK(d_b,v2sub,dimb,int_mb(k_b_offset),
     H                 (h3b_2 - 1 + (noab+nvab) * (h2b_2 - 1 + 
     N                 (noab+nvab) * (p6b_2 - 1 + (noab
     &                 +nvab) * (h7b_2 - 1)))))
               else
                  CALL GET_HASH_BLOCK_I(d_b,v2sub,
     &                 dimb,int_mb(k_b_offset),
     &                 (h3b_2 - 1 + (noab+nvab) * (h2b_2 - 1 + 
     N                 (noab+nvab) * (p6b_2 - 1 + (noab
     &                 +nvab) * (h7b_2 - 1)))),h3b_2,h2b_2,p6b_2,h7b_2)
               end if
            endif
         endif
         tce_hashnsort=.true.
      else      
         tce_hashnsort=.false.
      END IF
      return
      end
      logical function tce_hashnsort_2(dryrun,
     d     d_a, d_b,
     I     k_a_offset,k_b_offset,
     I     p4b,p7b,h1b,h2b,
     I     p5b,p6b,h3b,
     O     t2sub,v2sub
     P     )
      implicit none
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      integer k_a_offset ! [in]
      integer k_b_offset ! [in]
      logical dryrun ! [in]
      integer     p4b,p7b,h1b,h2b,
     I     p5b,p6b,h3b ! [in]
      double precision d_a(*),d_b(*) ! [in]
      double precision t2sub(*),v2sub(*) ! [out]
c
      integer p4b_1,p7b_1,h1b_1,h2b_1
      integer p5b_2,p6b_2,h3b_2,p7b_2
      integer l_a,k_a
      integer dima,dimb,dima_sort,dimb_sort,dim_common
c
      IF ((int_mb(k_spin+p4b-1)+int_mb(k_spin+p7b-1) .eq. 
     I     int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)).and.
     A (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p7b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t)) THEN
      dim_common = int_mb(k_range+p7b-1)
      dima_sort = int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * int_mb
     &(k_range+h2b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb
     &(k_range+h3b-1)
      dimb = dim_common * dimb_sort
      IF ((.not.dryrun).and.(dima .gt. 0) .and. (dimb .gt. 0)) THEN
      CALL TCE_RESTRICTED_4(p4b,p7b,h1b,h2b,p4b_1,p7b_1,h1b_1,h2b_1)
      CALL TCE_RESTRICTED_4(p5b,p6b,h3b,p7b,p5b_2,p6b_2,h3b_2,p7b_2)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsd_t_doubles_2',2,MA_ERR)
      IF ((p7b .lt. p4b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p4b_1 - noab - 1 + nvab * (p7b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),t2sub,int_mb(k_range+p7b-1)
     &,int_mb(k_range+p4b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,2,1,-1.0d0)
      END IF
      IF ((p4b .le. p7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p7b_1 - noab - 1 + nvab * (p4b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),t2sub,int_mb(k_range+p4b-1)
     &,int_mb(k_range+p7b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,1,2,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsd_t_doubles_2',3,MA_E
     &RR)
ccx      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
ccx     &ccsd_t_doubles_2',5,MA_ERR)
      IF ((h3b .le. p7b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,v2sub,dimb,
     &int_mb(k_b_offset),(p7b_2
     & - 1 + (noab+nvab) * (h3b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (p5b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,v2sub,
     &dimb,int_mb(k_b_offset),
     &(p7b_2
     & - 1 + (noab+nvab) * (h3b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (p5b_2 - 1)))),p7b_2,h3b_2,p6b_2,p5b_2)
      end if
ccx      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p5b-1)
ccx     &,int_mb(k_range+p6b-1),int_mb(k_range+h3b-1),int_mb(k_range+p7b-1)
ccx     &,3,2,1,4,1.0d0)
      END IF
      END IF

         tce_hashnsort_2=.true.
      else      
         tce_hashnsort_2=.false.
      END IF
      return
      end
c $Id$
