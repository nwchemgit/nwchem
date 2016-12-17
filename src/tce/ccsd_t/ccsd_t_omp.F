      SUBROUTINE ccsd_t_omp(d_t1,k_t1_offset,
     &                      d_t2,k_t2_offset,
     &                      d_v2,k_v2_offset,
     &                      energy1,energy2,size_t1)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "hybrid.fh"
#include "stdio.fh"
      integer, intent(in) :: d_t1, d_t2, d_v2
      integer, intent(in) :: k_t1_offset, k_t2_offset, k_v2_offset
      integer h1b
      integer h2b
      integer h3b
      integer p4b
      integer p5b
      integer p6b
      double precision, allocatable :: f_singles(:),f_doubles(:),f_t1(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_singles, f_doubles, f_t1
#endif
      integer alloc_error
      integer size_t3,i
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer size_t1
      double precision energy1,energy2
      double precision factor
      external nxtask
      logical nodezero
      integer off4,off5,off6
      integer off1,off2,off3
      !
      integer dim1,dim2,dim3,dim4,dim5,dim6
      integer spn1,spn2,spn3,spn4,spn5,spn6
      integer sym1,sym2,sym3,sym4,sym5,sym6
      integer dim12,dim34,dim56
      integer spn12,spn34,spn56
      integer sym12,sym34,sym56
      !
      nodezero=(ga_nodeid().eq.0)
      allocate( f_t1(1:size_t1), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: alloc t1',0,MA_ERR)
      endif
      call ma_zero(f_t1,size_t1)
      call ga_get(d_t1,1,size_t1,1,1,f_t1,size_t1)
      nprocs = GA_NNODES()
      !
      dim4=0
      do p4b = noab+1,noab+nvab
         dim4 = max(dim4,int_mb(k_range+p4b-1))
      enddo
      dim1=0
      do h1b = 1,noab
         dim1 = max(dim1,int_mb(k_range+h1b-1))
      enddo
      size_t3=(dim4**3)*(dim1**3)
      !
      allocate( f_singles(1:size_t3), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: MA error sgl',0,MA_ERR)
      endif
      allocate( f_doubles(1:size_t3), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: MA error dbl',0,MA_ERR)
      endif
      !
      count = 0
      next = nxtask(nprocs,1)
      energy1=0.0d0
      energy2=0.0d0
      do p4b = noab+1,noab+nvab
       dim4=int_mb(k_range+p4b-1)
       spn4=int_mb(k_spin +p4b-1)
       sym4=int_mb(k_sym  +p4b-1)
       off4 = k_evl_sorted+int_mb(k_offset+p4b-1)
       do p5b = p4b,noab+nvab
        dim5=int_mb(k_range+p5b-1)
        spn5=int_mb(k_spin +p5b-1)
        sym5=int_mb(k_sym  +p5b-1)
        off5 = k_evl_sorted+int_mb(k_offset+p5b-1)
        do p6b = p5b,noab+nvab
         dim6=int_mb(k_range+p6b-1)
         spn6=int_mb(k_spin +p6b-1)
         sym6=int_mb(k_sym  +p6b-1)
         off6 = k_evl_sorted+int_mb(k_offset+p6b-1)
         do h1b = 1,noab
          dim1=int_mb(k_range+h1b-1)
          spn1=int_mb(k_spin +h1b-1)
          sym1=int_mb(k_sym  +h1b-1)
          off1 = k_evl_sorted+int_mb(k_offset+h1b-1)
          do h2b = h1b,noab
           dim2=int_mb(k_range+h2b-1)
           spn2=int_mb(k_spin +h2b-1)
           sym2=int_mb(k_sym  +h2b-1)
           off2 = k_evl_sorted+int_mb(k_offset+h2b-1)
           do h3b = h2b,noab
            dim3=int_mb(k_range+h3b-1)
            spn3=int_mb(k_spin +h3b-1)
            sym3=int_mb(k_sym  +h3b-1)
            off3 = k_evl_sorted+int_mb(k_offset+h3b-1)
            if ((spn4+spn5+spn6).eq.(spn1+spn2+spn3)) then
            if ((.not.restricted).or.
     &          (spn4+spn5+spn6+spn1+spn2+spn3.le.8)) then
            if (ieor(sym4,ieor(sym5,ieor(sym6,ieor(sym1,
     &               ieor(sym2,sym3))))).eq.0) then
            if (next.eq.count) then
            size_t3 = dim4 * dim5 * dim6 * dim1 * dim2 * dim3
            call dfill(size_t3, 0.0d0, f_singles, 1)
            call dfill(size_t3, 0.0d0, f_doubles, 1)
            call ccsd_t_singles_omp(f_singles,f_t1,d_v2,
     1                              k_t1_offset,k_v2_offset,
     2                              h1b,h2b,h3b,p4b,p5b,p6b,2)
            call ccsd_t_doubles_omp(f_doubles,d_t2,d_v2,
     1                              k_t2_offset,k_v2_offset,
     2                              h1b,h2b,h3b,p4b,p5b,p6b,2)
            call ccsd_t_dot_omp(f_singles,f_doubles,
     &                          restricted,
     &                          h1b,h2b,h3b,p4b,p5b,p6b,
     &                          dbl_mb(off1),dbl_mb(off2),
     &                          dbl_mb(off3),dbl_mb(off4),
     &                          dbl_mb(off5),dbl_mb(off6),
     &                          dim1,dim2,dim3,dim4,dim5,dim6,
     &                          energy1,energy2)
            next = nxtask(nprocs,1)
            endif
            count = count + 1
            endif
            endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_dgop(1975,energy1,1,'+')
      call ga_dgop(1976,energy2,1,'+')
      deallocate( f_t1, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t t1',3,MA_ERR)
      endif
      deallocate( f_doubles, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t doubles',3,MA_ERR)
      endif
      deallocate( f_singles, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t singles',4,MA_ERR)
      endif
      return
      end
      subroutine ccsd_t_dot_omp(a_singles, a_doubles,
     &                          restricted,
     &                          h1b,h2b,h3b,p4b,p5b,p6b,
     &                          o_h1,o_h2,o_h3,
     &                          o_p4,o_p5,o_p6,
     &                          r_h1,r_h2,r_h3,
     &                          r_p4,r_p5,r_p6,
     &                          energy1,energy2)
      implicit none
      integer h1b, h2b, h3b, p4b, p5b, p6b
      double precision o_h1(*),o_h2(*),o_h3(*)
      double precision o_p4(*),o_p5(*),o_p6(*)
      integer r_h1,r_h2,r_h3
      integer r_p4,r_p5,r_p6
      double precision a_singles(r_h3,r_h2,r_h1,r_p6,r_p5,r_p4)
      double precision a_doubles(r_h3,r_h2,r_h1,r_p6,r_p5,r_p4)
      logical restricted
      double precision energy1,energy2
      integer h1, h2, h3, p4, p5, p6
      double precision factor,denom,sing,doub,denom_0
      double precision d_p4,d_p5,d_p6
      double precision d_h1,d_h2,d_h3
!
! factor = [ 1/36, 1/18, 1/12, 1/6, 1/4, 1/3, 1/2, 1, 2]
!
      if (restricted) then
        factor = 2.0d0
      else
        factor = 1.0d0
      endif
      if ((p4b.eq.p5b).and.(p5b.eq.p6b)) then
        factor = factor / 6.0d0
      else if ((p4b.eq.p5b).or.(p5b.eq.p6b)) then
        factor = factor / 2.0d0
      endif
      if ((h1b.eq.h2b).and.(h2b.eq.h3b)) then
        factor = factor / 6.0d0
      else if ((h1b.eq.h2b).or.(h2b.eq.h3b)) then
        factor = factor / 2.0d0
      endif
!$omp  parallel do collapse(3)
!$omp& schedule(static)
!$omp& default(shared)
!$omp& private(p4,p5,p6,h1,h2,h3,sing,doub,denom)
!$omp& private(d_p4,d_p5,d_p6,d_h1,d_h2,d_h3)
!$omp& reduction(+:energy1) reduction(+:energy2)
      do p4 = 1, r_p4
       do p5 = 1, r_p5
        do p6 = 1, r_p6
         d_p4 = o_p4(p4)
         d_p5 = o_p5(p5)
         d_p6 = o_p6(p6)
         denom_0=-(d_p4+d_p5+d_p6)
         do h1 = 1, r_h1
          d_h1 = o_h1(h1)
          do h2 = 1, r_h2
           d_h2 = o_h2(h2)
!dir$ simd
           do h3 = 1, r_h3
            d_h3 = o_h3(h3)
            sing = a_singles(h3,h2,h1,p6,p5,p4)
            doub = a_doubles(h3,h2,h1,p6,p5,p4)
            denom = doub * factor / ( d_h1+d_h2+d_h3+denom_0)
            energy1 = energy1 + denom*doub
            energy2 = energy2 + denom*(doub+sing)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      return
      end
      SUBROUTINE ccsd_t_singles_omp(a_i0,f_t1,d_v2,k_t1_offset,
     + k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,toggle)
!     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
!     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
!     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1 * P( 9 ) * t ( p4 h1 )_t * v ( p5 p6 h2 h3 )_v
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER t_p4b
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER toggle
      double precision f_t1(*)
      INTEGER k_t1_offset
      INTEGER d_v2
      INTEGER k_v2_offset
      DOUBLE PRECISION a_i0(*)
      IF(toggle .eq. 2) CALL ccsd_t_singles_omp_1(f_t1,k_t1_offset,d_v2,
     + k_v2_offset,a_i0,t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b)
      RETURN
      END
      SUBROUTINE ccsd_t_singles_omp_1(f_a,k_a_offset,d_b,k_b_offset,a_c,
     + t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b)
!     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
!     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
!     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1 * P( 9 ) * t ( p4 h1 )_t * v ( p5 p6 h2 h3 )_v
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      double precision f_a(*)
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER t_p4b
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER p4b
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p4b_1
      INTEGER h1b_1
      INTEGER p5b_2
      INTEGER p6b_2
      INTEGER h2b_2
      INTEGER h3b_2
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
      INTEGER a3(9,6)
      INTEGER ia6,ja6
      DOUBLE PRECISION a_c(*)
!
      a3(1,1)=t_p4b
      a3(1,2)=t_p5b
      a3(1,3)=t_p6b
      a3(1,4)=t_h1b
      a3(1,5)=t_h2b
      a3(1,6)=t_h3b
!
      a3(2,1)=t_p4b
      a3(2,2)=t_p5b
      a3(2,3)=t_p6b
      a3(2,4)=t_h2b
      a3(2,5)=t_h1b
      a3(2,6)=t_h3b
!
      a3(3,1)=t_p4b
      a3(3,2)=t_p5b
      a3(3,3)=t_p6b
      a3(3,4)=t_h3b
      a3(3,5)=t_h1b
      a3(3,6)=t_h2b
!
      a3(4,1)=t_p5b
      a3(4,2)=t_p4b
      a3(4,3)=t_p6b
      a3(4,4)=t_h1b
      a3(4,5)=t_h2b
      a3(4,6)=t_h3b
!
      a3(5,1)=t_p5b
      a3(5,2)=t_p4b
      a3(5,3)=t_p6b
      a3(5,4)=t_h2b
      a3(5,5)=t_h1b
      a3(5,6)=t_h3b
!
      a3(6,1)=t_p5b
      a3(6,2)=t_p4b
      a3(6,3)=t_p6b
      a3(6,4)=t_h3b
      a3(6,5)=t_h1b
      a3(6,6)=t_h2b
!
      a3(7,1)=t_p6b
      a3(7,2)=t_p4b
      a3(7,3)=t_p5b
      a3(7,4)=t_h1b
      a3(7,5)=t_h2b
      a3(7,6)=t_h3b
!
      a3(8,1)=t_p6b
      a3(8,2)=t_p4b
      a3(8,3)=t_p5b
      a3(8,4)=t_h2b
      a3(8,5)=t_h1b
      a3(8,6)=t_h3b
!
      a3(9,1)=t_p6b
      a3(9,2)=t_p4b
      a3(9,3)=t_p5b
      a3(9,4)=t_h3b
      a3(9,5)=t_h1b
      a3(9,6)=t_h2b
!
      do ia6=1,8
       if(a3(ia6,1).ne.0) then
       do ja6=ia6+1,9
        if((a3(ia6,1).eq.a3(ja6,1)).and.(a3(ia6,2).eq.a3(ja6,2))
     &   .and.(a3(ia6,3).eq.a3(ja6,3)).and.(a3(ia6,4).eq.a3(ja6,4))
     &   .and.(a3(ia6,5).eq.a3(ja6,5)).and.(a3(ia6,6).eq.a3(ja6,6)))
     &   then
         a3(ja6,1)=0
         a3(ja6,2)=0
         a3(ja6,3)=0
         a3(ja6,4)=0
         a3(ja6,5)=0
         a3(ja6,6)=0
        end if
       enddo
       end if
      enddo
!
      do ia6=1,9
      p4b=a3(ia6,1)
      p5b=a3(ia6,2)
      p6b=a3(ia6,3)
      h1b=a3(ia6,4)
      h2b=a3(ia6,5)
      h3b=a3(ia6,6)

      if((p5b.le.p6b).and.(h2b.le.h3b).and.(p4b.ne.0)) then
      !IF (.true.) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1
     &)+int_mb(k_spin+p6b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      IF (int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb(
     &k_sym+p6b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_v,irrep_t)) THEN
      dimc = int_mb(k_range+p4b-1) * int_mb(k_range+p5b-1) * int_mb(k_ra
     &nge+p6b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * int_m
     &b(k_range+h3b-1)
      IF (int_mb(k_spin+p4b-1) .eq. int_mb(k_spin+h1b-1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),int_mb(k_sym+h1b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p4b,h1b,p4b_1,h1b_1)
      CALL TCE_RESTRICTED_4(p5b,p6b,h2b,h3b,p5b_2,p6b_2,h2b_2,h3b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb
     &(k_range+h2b-1) * int_mb(k_range+h3b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsd_t_singles_omp_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsd_t_singles_omp_1',2,MA_ERR)
      CALL GET_HASH_BLOCK_MA(f_a,dbl_mb(k_a),dima,
     & int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p4b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p4b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsd_t_singles_omp_1',3,
     +MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsd_t_singles_omp_1',4,MA_ERR)
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b_sort),dimb,
     &int_mb(k_b_offset),(h3b_2
     & - 1 + (noab+nvab) * (h2b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (p5b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b_sort),dimb,
     &int_mb(k_b_offset),
     &(h3b_2
     & - 1 + (noab+nvab) * (h2b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (p5b_2 - 1)))),h3b_2,h2b_2,p6b_2,p5b_2)
      end if
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_s1_1(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_s1_2(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN

      call sd_t_s1_3(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_s1_4(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_s1_5(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      call sd_t_s1_6(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_s1_7(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_s1_8(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      call sd_t_s1_9(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF (.not.MA_POP_STACK(l_b_sort))
     + CALL ERRQUIT('ccsd_t_singles_omp_1',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort))
     + CALL ERRQUIT('ccsd_t_singles_omp_1',8,MA_ERR)
!
      END IF
      END IF
      END IF
!
      END IF
      END IF
      END IF
      !END IF
!
      END IF
      END DO
      RETURN
      END
      SUBROUTINE ccsd_t_doubles_omp(a_i0,d_t2,d_v2,k_t2_offset,
     + k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,toggle)
!     $Id: ccsd_t_doubles_omp.F 24345 2013-06-22 06:19:30Z d3y133 $
!     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
!     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
!     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( h7 ) * t ( p4 p5 h1 h7 )_t * v ( h7 p6 h2 h3 )_v
!     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( p7 ) * t ( p4 p7 h1 h2 )_t * v ( p5 p6 h3 p7 )_v
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "sym.fh"
      INTEGER t_p4b
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER toggle
      INTEGER d_t2
      INTEGER k_t2_offset
      INTEGER d_v2
      INTEGER k_v2_offset
      integer all_t1
      integer all_t2
      DOUBLE PRECISION a_i0(*)
      IF(toggle .eq. 2) CALL ccsd_t_doubles_omp_1(d_t2,k_t2_offset,d_v2,
     + k_v2_offset,a_i0,t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b)
      IF(toggle .eq. 2) CALL ccsd_t_doubles_omp_2(d_t2,k_t2_offset,d_v2,
     + k_v2_offset,a_i0,t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b)
      RETURN
      END
      SUBROUTINE ccsd_t_doubles_omp_1(d_a,k_a_offset,d_b,k_b_offset,a_c,
     + t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b)
!     $Id: ccsd_t_doubles_omp.F 24345 2013-06-22 06:19:30Z d3y133 $
!     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
!     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
!     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( h7 ) * t ( p4 p5 h1 h7 )_t * v ( h7 p6 h2 h3 )_v
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "util.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER t_p4b
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER p4b
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER h7b
      INTEGER p4b_1
      INTEGER p5b_1
      INTEGER h1b_1
      INTEGER h7b_1
      INTEGER p6b_2
      INTEGER h7b_2
      INTEGER h2b_2
      INTEGER h3b_2
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
      INTEGER a3(9,6)
      INTEGER ia6,ja6
      DOUBLE PRECISION a_c(*)
!
      a3(1,1)=t_p4b
      a3(1,2)=t_p5b
      a3(1,3)=t_p6b
      a3(1,4)=t_h1b
      a3(1,5)=t_h2b
      a3(1,6)=t_h3b
!
      a3(2,1)=t_p4b
      a3(2,2)=t_p5b
      a3(2,3)=t_p6b
      a3(2,4)=t_h2b
      a3(2,5)=t_h1b
      a3(2,6)=t_h3b
!
      a3(3,1)=t_p4b
      a3(3,2)=t_p5b
      a3(3,3)=t_p6b
      a3(3,4)=t_h3b
      a3(3,5)=t_h1b
      a3(3,6)=t_h2b
!
      a3(4,1)=t_p5b
      a3(4,2)=t_p6b
      a3(4,3)=t_p4b
      a3(4,4)=t_h1b
      a3(4,5)=t_h2b
      a3(4,6)=t_h3b
!
      a3(5,1)=t_p5b
      a3(5,2)=t_p6b
      a3(5,3)=t_p4b
      a3(5,4)=t_h2b
      a3(5,5)=t_h1b
      a3(5,6)=t_h3b
!
      a3(6,1)=t_p5b
      a3(6,2)=t_p6b
      a3(6,3)=t_p4b
      a3(6,4)=t_h3b
      a3(6,5)=t_h1b
      a3(6,6)=t_h2b
!
      a3(7,1)=t_p4b
      a3(7,2)=t_p6b
      a3(7,3)=t_p5b
      a3(7,4)=t_h1b
      a3(7,5)=t_h2b
      a3(7,6)=t_h3b
!
      a3(8,1)=t_p4b
      a3(8,2)=t_p6b
      a3(8,3)=t_p5b
      a3(8,4)=t_h2b
      a3(8,5)=t_h1b
      a3(8,6)=t_h3b
!
      a3(9,1)=t_p4b
      a3(9,2)=t_p6b
      a3(9,3)=t_p5b
      a3(9,4)=t_h3b
      a3(9,5)=t_h1b
      a3(9,6)=t_h2b
!
      do ia6=1,8
       if(a3(ia6,1).ne.0) then
       do ja6=ia6+1,9
        if((a3(ia6,1).eq.a3(ja6,1)).and.(a3(ia6,2).eq.a3(ja6,2))
     &   .and.(a3(ia6,3).eq.a3(ja6,3)).and.(a3(ia6,4).eq.a3(ja6,4))
     &   .and.(a3(ia6,5).eq.a3(ja6,5)).and.(a3(ia6,6).eq.a3(ja6,6)))
     &   then
         a3(ja6,1)=0
         a3(ja6,2)=0
         a3(ja6,3)=0
         a3(ja6,4)=0
         a3(ja6,5)=0
         a3(ja6,6)=0
        end if
       enddo
       end if
      enddo
!
      do ia6=1,9
      p4b=a3(ia6,1)
      p5b=a3(ia6,2)
      p6b=a3(ia6,3)
      h1b=a3(ia6,4)
      h2b=a3(ia6,5)
      h3b=a3(ia6,6)
      if((p4b.le.p5b).and.(h2b.le.h3b).and.(p4b.ne.0)) then
      IF (.true.) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1
     &)+int_mb(k_spin+p6b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      IF (int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb(
     &k_sym+p6b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_v,irrep_t)) THEN
      dimc = int_mb(k_range+p4b-1) * int_mb(k_range+p5b-1) * int_mb(k_ra
     &nge+p6b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * int_m
     &b(k_range+h3b-1)
!
      DO h7b = 1,noab
      IF (int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h7b-1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h7b-1)))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_4(p4b,p5b,h1b,h7b,p4b_1,p5b_1,h1b_1,h7b_1)
      CALL TCE_RESTRICTED_4(p6b,h7b,h2b,h3b,p6b_2,h7b_2,h2b_2,h3b_2)
      dim_common = int_mb(k_range+h7b-1)
      dima_sort = int_mb(k_range+p4b-1) * int_mb(k_range+p5b-1) * int_mb
     &(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p6b-1) * int_mb(k_range+h2b-1) * int_mb
     &(k_range+h3b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsd_t_doubles_omp_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsd_t_doubles_omp_1',2,MA_ERR)
      IF ((h7b .lt. h1b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (h7b_1 - 1 + noab * (p5b_1 - noab - 1 + nvab * (p4b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p4b-1)
     &,int_mb(k_range+p5b-1),int_mb(k_range+h7b-1),int_mb(k_range+h1b-1)
     &,4,2,1,3,-1.0d0)
      END IF
      IF ((h1b .le. h7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h7b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p5b_1 - noab - 1 + nvab * (p4b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p4b-1)
     &,int_mb(k_range+p5b-1),int_mb(k_range+h1b-1),int_mb(k_range+h7b-1)
     &,3,2,1,4,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsd_t_doubles_omp_1',3,
     + MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsd_t_doubles_omp_1',4,MA_ERR)
      IF ((h7b .le. p6b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b_sort),
     &dimb,int_mb(k_b_offset),(h3b_2
     & - 1 + (noab+nvab) * (h2b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (h7b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b_sort),
     &dimb,int_mb(k_b_offset),
     &(h3b_2
     & - 1 + (noab+nvab) * (h2b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (h7b_2 - 1)))),h3b_2,h2b_2,p6b_2,h7b_2)
      end if
      END IF

      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d1_1(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d1_2(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      call sd_t_d1_3(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p6b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p5b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d1_4(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p6b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p5b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d1_5(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p6b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p5b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      call sd_t_d1_6(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p5b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d1_7(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p5b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d1_8(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p5b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      call sd_t_d1_9(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+h7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
!
      IF (.not.MA_POP_STACK(l_b_sort))
     + CALL ERRQUIT('ccsd_t_doubles_omp_1',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort))
     + CALL ERRQUIT('ccsd_t_doubles_omp_1',8,MA_ERR)
!
      END IF
      END IF
      END IF
      END DO !h7b loop
!
      END IF
      END IF
      END IF
      END IF
!
      END IF
      END DO
      RETURN
      END
      SUBROUTINE ccsd_t_doubles_omp_2(d_a,k_a_offset,d_b,k_b_offset,a_c,
     + t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b)
!     $Id: ccsd_t_doubles_omp.F 24345 2013-06-22 06:19:30Z d3y133 $
!     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
!     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
!     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( p7 ) * t ( p4 p7 h1 h2 )_t * v ( p5 p6 h3 p7 )_v
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "util.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER t_p4b
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER p4b
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p7b
      INTEGER p4b_1
      INTEGER p7b_1
      INTEGER h1b_1
      INTEGER h2b_1
      INTEGER p5b_2
      INTEGER p6b_2
      INTEGER h3b_2
      INTEGER p7b_2
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
      INTEGER a3(9,6)
      INTEGER ia6,ja6
      DOUBLE PRECISION a_c(*)
!
        a3(1,1)=t_p4b
        a3(1,2)=t_p5b
        a3(1,3)=t_p6b
        a3(1,4)=t_h1b
        a3(1,5)=t_h2b
        a3(1,6)=t_h3b
!
        a3(2,1)=t_p4b
        a3(2,2)=t_p5b
        a3(2,3)=t_p6b
        a3(2,4)=t_h2b
        a3(2,5)=t_h3b
        a3(2,6)=t_h1b
!
        a3(3,1)=t_p4b
        a3(3,2)=t_p5b
        a3(3,3)=t_p6b
        a3(3,4)=t_h1b
        a3(3,5)=t_h3b
        a3(3,6)=t_h2b
!
        a3(4,1)=t_p5b
        a3(4,2)=t_p4b
        a3(4,3)=t_p6b
        a3(4,4)=t_h1b
        a3(4,5)=t_h2b
        a3(4,6)=t_h3b
!
        a3(5,1)=t_p5b
        a3(5,2)=t_p4b
        a3(5,3)=t_p6b
        a3(5,4)=t_h2b
        a3(5,5)=t_h3b
        a3(5,6)=t_h1b
!
        a3(6,1)=t_p5b
        a3(6,2)=t_p4b
        a3(6,3)=t_p6b
        a3(6,4)=t_h1b
        a3(6,5)=t_h3b
        a3(6,6)=t_h2b
!
        a3(7,1)=t_p6b
        a3(7,2)=t_p4b
        a3(7,3)=t_p5b
        a3(7,4)=t_h1b
        a3(7,5)=t_h2b
        a3(7,6)=t_h3b
!
        a3(8,1)=t_p6b
        a3(8,2)=t_p4b
        a3(8,3)=t_p5b
        a3(8,4)=t_h2b
        a3(8,5)=t_h3b
        a3(8,6)=t_h1b
!
        a3(9,1)=t_p6b
        a3(9,2)=t_p4b
        a3(9,3)=t_p5b
        a3(9,4)=t_h1b
        a3(9,5)=t_h3b
        a3(9,6)=t_h2b
!
        do ia6=1,8
         if(a3(ia6,1).ne.0) then
         do ja6=ia6+1,9
          if((a3(ia6,1).eq.a3(ja6,1)).and.(a3(ia6,2).eq.a3(ja6,2))
     &     .and.(a3(ia6,3).eq.a3(ja6,3)).and.(a3(ia6,4).eq.a3(ja6,4))
     &     .and.(a3(ia6,5).eq.a3(ja6,5)).and.(a3(ia6,6).eq.a3(ja6,6)))
     &     then
           a3(ja6,1)=0
           a3(ja6,2)=0
           a3(ja6,3)=0
           a3(ja6,4)=0
           a3(ja6,5)=0
           a3(ja6,6)=0
          end if
         enddo 
         end if
        enddo 
!
      do ia6=1,9
      p4b=a3(ia6,1)
      p5b=a3(ia6,2)
      p6b=a3(ia6,3)
      h1b=a3(ia6,4)
      h2b=a3(ia6,5)
      h3b=a3(ia6,6)
      if((p5b.le.p6b).and.(h1b.le.h2b).and.(p4b.ne.0)) then
      !IF (.true.) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1
     &)+int_mb(k_spin+p6b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      IF (int_mb(k_spin+p4b-1)+int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb(
     &k_sym+p6b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_v,irrep_t)) THEN
      dimc = int_mb(k_range+p4b-1) * int_mb(k_range+p5b-1) * int_mb(k_ra
     &nge+p6b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * int_m
     &b(k_range+h3b-1)
      DO p7b = noab+1,noab+nvab
      IF (int_mb(k_spin+p4b-1)+int_mb(k_spin+p7b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p4b-1),ieor(int_mb(k_sym+p7b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_4(p4b,p7b,h1b,h2b,p4b_1,p7b_1,h1b_1,h2b_1)
      CALL TCE_RESTRICTED_4(p5b,p6b,h3b,p7b,p5b_2,p6b_2,h3b_2,p7b_2)
      dim_common = int_mb(k_range+p7b-1)
      dima_sort = int_mb(k_range+p4b-1) * int_mb(k_range+h1b-1) * int_mb
     &(k_range+h2b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb
     &(k_range+h3b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsd_t_doubles_omp_2',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsd_t_doubles_omp_2',2,MA_ERR)
      IF ((p7b .lt. p4b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p4b_1 - noab - 1 + nvab * (p7b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p7b-1)
     &,int_mb(k_range+p4b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,2,1,-1.0d0)
      END IF
      IF ((p4b .le. p7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p7b_1 - noab - 1 + nvab * (p4b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p4b-1)
     &,int_mb(k_range+p7b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,1,2,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a))
     + CALL ERRQUIT('ccsd_t_doubles_omp_2',3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsd_t_doubles_omp_2',4,MA_ERR)
      IF ((h3b .le. p7b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b_sort),dimb,
     &int_mb(k_b_offset),(p7b_2
     & - 1 + (noab+nvab) * (h3b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (p5b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b_sort),
     &dimb,int_mb(k_b_offset),
     &(p7b_2
     & - 1 + (noab+nvab) * (h3b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
     &+nvab) * (p5b_2 - 1)))),p7b_2,h3b_2,p6b_2,p5b_2)
      end if
      END IF

      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d2_1(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &2b)) THEN
      call sd_t_d2_2(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &2b)) THEN
      call sd_t_d2_3(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d2_4(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &2b)) THEN
      call sd_t_d2_5(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &2b)) THEN
      call sd_t_d2_6(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      call sd_t_d2_7(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &2b)) THEN
      call sd_t_d2_8(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &2b)) THEN
      call sd_t_d2_9(int_mb(k_range+h3b-1),int_mb(k_range+h2b-1),
     1               int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),
     2               int_mb(k_range+p5b-1),int_mb(k_range+p4b-1),
     3               int_mb(k_range+p7b-1),
     4               a_c,dbl_mb(k_a_sort),dbl_mb(k_b_sort))
      END IF

      IF (.not.MA_POP_STACK(l_b_sort))
     + CALL ERRQUIT('ccsd_t_doubles_omp_2',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort))
     + CALL ERRQUIT('ccsd_t_doubles_omp_2',8,MA_ERR)


      END IF
      END IF
      END IF
      END DO ! p7b loop
!
      END IF
      END IF
      END IF
      !END IF
!
      END IF
      END DO

      RETURN
      END

