      SUBROUTINE cr_qcisd_t(d_t1,k_t1_offset,size_t1,
     1                      d_t2,k_t2_offset,size_t2,
     2                      d_f1,k_f1_offset,d_v2,k_v2_offset,delta)
c     4                      e_singles,   ! < T_1+ | M_1 >
c     5                      e_doubles,   ! < C_2 | M_2 + M_1 * T_1 >
c     6                      e_trip_c,    ! < R_0[ V_N T_1 ] + R_0[ V_N T_2 ] | M_3>
c     7                      e_trip_u)    ! < R_0[ V_N T_1 ] + R_0[ V_N T_2 ] | M_2 * T_1 + M_1 * C_2 >
c                                         ! where C_2 = T_2 + (1/2) * T_1^2
c
c     $Id$
c
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1,k_t1_offset
      integer d_t2,k_t2_offset
      integer d_m1,d_m2,d_c2,d_d2
      integer d_f1,k_f1_offset
      integer d_v2,k_v2_offset
      integer d_e,k_e_offset,l_e_offset,size_e
      integer size_t1,size_t2
      double precision e_singles
      double precision e_doubles,e_doubles_1,e_doubles_2
      double precision e_trip_u,e_trip_c
      double precision e_trip_u_1_1,e_trip_u_1_2
      double precision e_trip_u_2_1,e_trip_u_2_2
      double precision delta
      double precision den_sing,den_doub,den_trip
      double precision ddotfile
      external ddotfile
      character*255 filename
c
c zero out target components of delta
c
      e_singles    = 0.0d0
      e_doubles    = 0.0d0
      e_doubles_1  = 0.0d0
      e_doubles_2  = 0.0d0
      e_trip_c     = 0.0d0
      e_trip_u     = 0.0d0
      e_trip_u_1_1 = 0.0d0
      e_trip_u_1_2 = 0.0d0
      e_trip_u_2_1 = 0.0d0
      e_trip_u_2_2 = 0.0d0
      den_sing     = 0.0d0
      den_doub     = 0.0d0
      den_trip     = 0.0d0
c
c scalar accumulation target
c
      call tce_filename('delta',filename)
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      call createfile(filename,d_e,size_e)
      call tce_zero(d_e,size_e)
c
c form M_1
c
c      print*,'form M_1'
      call tce_filename('m1',filename)
      call createfile(filename,d_m1,size_t1)
      call ccsd_t1(d_f1,d_m1,d_t1,d_t2,d_v2,
     1             k_f1_offset,k_t1_offset,
     2             k_t1_offset,k_t2_offset,k_v2_offset)
c
c form M_2
c
c      print*,'form M_2'
      call tce_filename('m2',filename)
      call createfile(filename,d_m2,size_t2)
      call ccsd_t2(d_f1,d_m2,d_t1,d_t2,d_v2,
     1             k_f1_offset,k_t2_offset,
     2             k_t1_offset,k_t2_offset,k_v2_offset)
c
c C_2 = T_2 + (1/2) * T_1^2
c
c      print*,'form C_2'
      call tce_filename('c2',filename)
      call createfile(filename,d_c2,size_t2)
      call copyfile(d_t2,d_c2,size_t2)
      call c2f_t2_t12(d_t1,k_t1_offset,d_c2,k_t2_offset) ! plus
c
c < T_1 | M_1 >
c
      e_singles = ddotfile(d_t1,d_m1,size_t1)
c
c < T_1 | T_1 > 
c
      den_sing  = ddotfile(d_t1,d_t1,size_t1)
c
c < C_2 | M_1 * T_1 >
c
c      print*,'calling m2_dot_m1t1'
      call m2_dot_m1t1(d_m1,d_m2,d_e,d_t1,
     1                 k_t1_offset,k_t2_offset,
     2                 k_e_offset,k_t1_offset)
      call reconcilefile(d_e,1)
      call get_block(d_e,e_doubles_1,1,0)
      call tce_zero(d_e,size_e)
c
c < C_2 | M_2 >
c
c      print*,'calling c2_d2'
      call c2_d2(d_c2,d_m2,d_e,k_t2_offset,k_t2_offset,k_e_offset)
      call reconcilefile(d_e,1)
      call get_block(d_e,e_doubles_2,1,0)
      call tce_zero(d_e,size_e)
c
c < C_2 | M_1 * T_1 > + < C_2 | M_2 >
c
      e_doubles = e_doubles_1 + e_doubles_2
c
c < C_2 | C_2 >
c
c      print*,'calling c2_d2'
      call c2_d2(d_c2,d_c2,d_e,k_t2_offset,k_t2_offset,k_e_offset)
      call reconcilefile(d_e,1)
      call get_block(d_e,den_doub,1,0)
      call tce_zero(d_e,size_e)
c
c < R_0[ V_N T_1 ] + R_0[ V_N T_2 ] | M_2 * T_1 + M_1 * C_2 >
c  plus CR-CCSD(T) terms
c
c      print*,'calling cr_qcisd_t_U'
      call cr_qcisd_t_U(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                  d_c2,d_m1,d_m2,
     3                  e_trip_u_1_1,e_trip_u_1_2,
     4                  e_trip_u_2_1,e_trip_u_2_2,
     5                  e_trip_c,den_trip,size_t1)
      e_trip_u = e_trip_u_1_1 + e_trip_u_1_2
     1         + e_trip_u_2_1 + e_trip_u_2_2
c
c sum correction
c
#ifdef DEBUG_PRINT
      write(6,10) 'non-CC singles: < T_1 | M_1 >',e_singles
      write(6,10) 'non-CC doubles: < C_2 | M_1 * T_1 >',e_doubles_1
      write(6,10) 'non-CC doubles: < C_2 | M_2 >',e_doubles_2
      write(6,10) 'non-CC doubles: total',e_doubles
      write(6,10) 'non-CC triples: < R_0[ V_N T_1 ] | M_1 * C_2 >',
     1             e_trip_u_1_1
      write(6,10) 'non-CC triples: < R_0[ V_N T_1 ] | M_2 * T_1 >',
     1             e_trip_u_1_2
      write(6,10) 'non-CC triples: < R_0[ V_N T_2 ] | M_1 * C_2 >',
     1             e_trip_u_2_1
      write(6,10) 'non-CC triples: < R_0[ V_N T_2 ] | M_2 * T_1 >',
     1             e_trip_u_2_2
      write(6,10) 'non-CC triples: total',e_trip_u
      write(6,10) ' CR-CC triples: total',e_trip_c
      write(6,*)
      write(6,11) 'singles: < T_1 | T_1 >',den_sing
      write(6,11) 'doubles: < C_2 | C_2 >',den_doub
      write(6,11) 'triples:   CR-CCSD(T) ',den_trip
   10 format(1x,'energy contribution from',1x,a,1x,f25.15)
   11 format(1x,'denominator from',1x,a,1x,f25.15)
#endif
      delta = ( e_singles + e_doubles + e_trip_u + e_trip_c )
     1      / (   1.0d0   + den_sing  + den_doub + den_trip )
c
c clean up temporaries
c
      call deletefile(d_e)
      call deletefile(d_m1)
      call deletefile(d_m2)
      call deletefile(d_c2)
      if (.not.ma_pop_stack(l_e_offset)) 
     1        call errquit('cr_qcisd_t',1,MA_ERR)
c
      return
      end



      SUBROUTINE cr_qcisd_t_U(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                        d_c2,d_m1,d_m2,
     3                        e1_1,e1_2,e2_1,e2_2,e3_1,den2,size_t1)
c
c     ccsd_t_singles:
c     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1 * P( 9 ) * t ( p4 h1 )_t * v ( p5 p6 h2 h3 )_v
c
c k_singles = R_0[ V_N T_1 ]
c k_doubles = R_0[ V_N T_2 ]
c k_right_1 = M_1 * C_2
c k_right_2 = M_2 * T_1
c
c e1_1 = < R_0[ V_N T_1 ] | M_1 * C_2 >
c e1_2 = < R_0[ V_N T_1 ] | M_2 * T_1 >
c e2_1 = < R_0[ V_N T_2 ] | M_1 * C_2 >
c e2_2 = < R_0[ V_N T_2 ] | M_2 * T_1 >
c
c e3_1 = CR-CCSD(T) numerator
c den  = CR-CCSD(T) denominator
c
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1,k_t1_offset
      integer d_t2,k_t2_offset
      integer d_m1,d_m2,d_c2,d_d2
      integer d_f1,k_f1_offset
      integer d_v2,k_v2_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_singles,l_singles
      integer k_doubles,l_doubles
      integer k_right_1,l_right_1
      integer k_right_2,l_right_2
      integer k_right_3,l_right_3
      integer k_denom,l_denom
      integer d_i1_1,d_i1_2,d_i1_3
      integer k_i1_offset_1,l_i1_offset_1
      integer k_i1_offset_2,l_i1_offset_2
      integer k_i1_offset_3,l_i1_offset_3
      integer size,i
      integer g_energy
      integer nxtask
      integer next
      integer nprocs
      integer count
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
c ---------------------------------------
      double precision e1_1,e1_2,e2_1,e2_2,e3_1,den2
      double precision factor
      external nxtask
c
c - T1/X1 LOCALIZATION ----------
c    opening l_t1_local and l_x1_local
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
c    copy d_t1 ==> l_t1_local
cc        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
      call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
c -------------------------------
c
c     This is used to accumulate contribution from all nodes.
c
      if (.not.ga_create(mt_dbl,3,2,'perturbative',1,1,g_energy))
     1  call errquit('cr_qcisd_t_U: GA problem',0,GA_ERR)
c
      call cr_ccsd_t_N(dbl_mb(k_right_3),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      call cr_ccsd_t_E(dbl_mb(k_denom),d_i1_3,
     1  k_t1_local,d_t2,k_i1_offset_3,k_t1_offset,k_t2_offset,
     2  l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      e1_1 = 0.0d0
      e1_2 = 0.0d0
      e2_1 = 0.0d0
      e2_2 = 0.0d0
      e3_1 = 0.0d0
      den2 = 0.0d0
      do t_p4b = noab+1,noab+nvab
       do t_p5b = t_p4b,noab+nvab
        do t_p6b = t_p5b,noab+nvab
         do t_h1b = 1,noab
          do t_h2b = t_h1b,noab
           do t_h3b = t_h2b,noab
            if (next.eq.count) then                        
            if (int_mb(k_spin+t_p4b-1)
     1         +int_mb(k_spin+t_p5b-1)
     2         +int_mb(k_spin+t_p6b-1)
     3      .eq.int_mb(k_spin+t_h1b-1)
     4         +int_mb(k_spin+t_h2b-1)
     5         +int_mb(k_spin+t_h3b-1)) then
            if ((.not.restricted).or.
     1         (int_mb(k_spin+t_p4b-1)
     1         +int_mb(k_spin+t_p5b-1)
     2         +int_mb(k_spin+t_p6b-1)
     3         +int_mb(k_spin+t_h1b-1)
     4         +int_mb(k_spin+t_h2b-1)
     5         +int_mb(k_spin+t_h3b-1).le.8)) then
            if (ieor(int_mb(k_sym+t_p4b-1),
     1          ieor(int_mb(k_sym+t_p5b-1),
     2          ieor(int_mb(k_sym+t_p6b-1),
     3          ieor(int_mb(k_sym+t_h1b-1),
     4          ieor(int_mb(k_sym+t_h2b-1),
     5               int_mb(k_sym+t_h3b-1)))))).eq.0) then
            size = int_mb(k_range+t_p4b-1)
     1           * int_mb(k_range+t_p5b-1)
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1       k_singles)) call errquit('cr_qcisd_t_U: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1       k_doubles)) call errquit('cr_qcisd_t_U: MA error',2,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) right_1',l_right_1,
     1       k_right_1)) call errquit('cr_qcisd_t_U: MA error',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) right_2',l_right_2,
     1       k_right_2)) call errquit('cr_qcisd_t_U: MA error',4,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) right_3',l_right_3,
     1       k_right_3)) call errquit('cr_qcisd_t_U: MA error',5,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'denominator',l_denom,
     1       k_denom)) call errquit('cr_qcisd_t_U: MA error',6,MA_ERR)
            do i = 1, size ! replace with ma_zero(dbl_mb(k_a),size_a)
             dbl_mb(k_singles+i-1) = 0.0d0
             dbl_mb(k_doubles+i-1) = 0.0d0
             dbl_mb(k_right_1+i-1) = 0.0d0
             dbl_mb(k_right_2+i-1) = 0.0d0
             dbl_mb(k_right_3+i-1) = 0.0d0
             dbl_mb(k_denom+i-1)   = 0.0d0
            enddo
c      print*,'calling ccsd_t_singles'
            call ccsd_t_singles(dbl_mb(k_singles),d_t1,d_v2,
     1        k_t1_offset,k_v2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c      print*,'calling ccsd_t_doubles'
            call ccsd_t_doubles(dbl_mb(k_doubles),d_t2,d_v2,
     1        k_t2_offset,k_v2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c      print*,'calling cr_qcisd_t_right'
            call cr_qcisd_t_right(dbl_mb(k_right_1),
     1        d_m1,k_t1_offset,d_c2,k_t2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b)
c      print*,'calling cr_qcisd_t_right'
            call cr_qcisd_t_right(dbl_mb(k_right_2),
     1        d_t1,k_t1_offset,d_m2,k_t2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b)
c      print*,'calling cr_ccsd_t'
            call cr_ccsd_t_N(dbl_mb(k_right_3),d_f1,d_i1_1,d_i1_2,
     1        k_t1_local,d_t2,d_v2,
     1        k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2        k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3        l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c      print*,'calling cr_ccsd_t_E'
            call cr_ccsd_t_E(dbl_mb(k_denom),d_i1_3,
     1        k_t1_local,d_t2,k_i1_offset_3,k_t1_offset,k_t2_offset,
     2        l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            if (restricted) then
              factor = 2.0d0
            else
              factor = 1.0d0
            endif
            if ((t_p4b.eq.t_p5b).and.(t_p5b.eq.t_p6b)) then
              factor = factor / 6.0d0
            else if ((t_p4b.eq.t_p5b).or.(t_p5b.eq.t_p6b)) then
              factor = factor / 2.0d0
            endif
            if ((t_h1b.eq.t_h2b).and.(t_h2b.eq.t_h3b)) then
              factor = factor / 6.0d0
            else if ((t_h1b.eq.t_h2b).or.(t_h2b.eq.t_h3b)) then
              factor = factor / 2.0d0
            endif
            i = 0
            do t_p4 = 1, int_mb(k_range+t_p4b-1)
             do t_p5 = 1, int_mb(k_range+t_p5b-1)
              do t_p6 = 1, int_mb(k_range+t_p6b-1)
               do t_h1 = 1, int_mb(k_range+t_h1b-1)
                do t_h2 = 1, int_mb(k_range+t_h2b-1)
                 do t_h3 = 1, int_mb(k_range+t_h3b-1)
                  i = i + 1
c
c                 e1_1 = < R_0[ V_N T_1 ] | M_1 * C_2 >
c
                  e1_1 = e1_1 + factor * dbl_mb(k_singles+i-1)
     1                                 * dbl_mb(k_right_1+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c
c                 e1_2 = < R_0[ V_N T_1 ] | M_2 * T_1 >
c
                  e1_2 = e1_2 + factor * dbl_mb(k_singles+i-1)
     1                                 * dbl_mb(k_right_2+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c
c                 e2_1 = < R_0[ V_N T_2 ] | M_1 * C_2 >
c
                  e2_1 = e2_1 + factor * dbl_mb(k_doubles+i-1)
     1                                 * dbl_mb(k_right_1+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c
c                 e2_2 = < R_0[ V_N T_2 ] | M_2 * T_1 >
c
                  e2_2 = e2_2 + factor * dbl_mb(k_doubles+i-1)
     1                                 * dbl_mb(k_right_2+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c
c                 e3_1 = CR-CCSD(T) correction
c
                  e3_1 = e3_1 + factor * dbl_mb(k_right_3+i-1)
     1             * (dbl_mb(k_singles+i-1) + dbl_mb(k_doubles+i-1))
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c
c                 den2 = < R( | exp( T_1 + T_2 ) >
c
                  den2 = den2 + factor * dbl_mb(k_denom+i-1)
     1             * (dbl_mb(k_singles+i-1) + dbl_mb(k_doubles+i-1))
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
            if (.not.MA_POP_STACK(l_denom)) 
     1        call errquit('cr_qcisd_t_U',16,MA_ERR)
            if (.not.MA_POP_STACK(l_right_3)) 
     1        call errquit('cr_qcisd_t_U',15,MA_ERR)
            if (.not.MA_POP_STACK(l_right_2)) 
     1        call errquit('cr_qcisd_t_U',14,MA_ERR)
            if (.not.MA_POP_STACK(l_right_1)) 
     1        call errquit('cr_qcisd_t_U',13,MA_ERR)
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('cr_qcisd_t_U',12,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('cr_qcisd_t_U',11,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
            endif
            count = count + 1
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      next = nxtask(-nprocs,1)
      call cr_ccsd_t_E(dbl_mb(k_denom),d_i1_3,
     1  k_t1_local,d_t2,k_i1_offset_3,k_t1_offset,k_t2_offset,
     2  l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      call cr_ccsd_t_N(dbl_mb(k_right_3),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
c
c     This is used to accumulate contribution from all nodes.
c     Would not a parallel reduction be better?
c
      call ga_zero(g_energy)
c
      call ga_acc(g_energy,1,1,1,1,e1_1,1,1.0d0)
      call ga_acc(g_energy,1,1,2,2,e1_2,1,1.0d0)
      call ga_acc(g_energy,2,2,1,1,e2_1,1,1.0d0)
      call ga_acc(g_energy,2,2,2,2,e2_2,1,1.0d0)
      call ga_acc(g_energy,3,3,1,1,e3_1,1,1.0d0)
      call ga_acc(g_energy,3,3,2,2,den2,1,1.0d0)
c
      call ga_sync()
c
      call ga_get(g_energy,1,1,1,1,e1_1,1)
      call ga_get(g_energy,1,1,2,2,e1_2,1)
      call ga_get(g_energy,2,2,1,1,e2_1,1)
      call ga_get(g_energy,2,2,2,2,e2_2,1)
      call ga_get(g_energy,3,3,1,1,e3_1,1)
      call ga_get(g_energy,3,3,2,2,den2,1)
c
      if (.not.ga_destroy(g_energy))
     1  call errquit('cr_qcisd_t_U: GA problem',1,GA_ERR)
c
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end




      SUBROUTINE cr_qcisd_t_right(a_c,d_a,k_a_offset,d_b,k_b_offset,
     2                            t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b)
C
C     i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1 * P( 9 ) * a ( p4 h1 )_t * b ( p5 p6 h2 h3 )_v
C
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
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
      DOUBLE PRECISION a_c(*)
      LOGICAL skipped
      DO p4b = noab+1,noab+nvab
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      skipped = .true.
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) skipped = .false.
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) skipped = .false.
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) skipped = .false.
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) skipped = .false.
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) skipped = .false.
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) skipped = .false.
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) skipped = .false.
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) skipped = .false.
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) skipped = .false.
      IF (.not.skipped) THEN
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
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsd_t_singles_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
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
     & ERRQUIT('ccsd_t_singles_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsd_t_singles_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p4b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p4b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsd_t_singles_1',3,MA_E
     &RR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsd_t_singles_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsd_t_singles_1',5,MA_ERR)
!
!     OLD CALL FOR V2
!
!      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),
!     &(h3b_2
!     & - 1 + (noab+nvab) * (h2b_2 - 1 + (noab+nvab) * (p6b_2 - 1 + (noab
!     &+nvab) * (p5b_2 - 1)))),h3b_2,h2b_2,p6b_2,p5b_2)
!
!     CALL FOR T2
!
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h3b_2
     & - 1 + noab * (h2b_2 - 1 + noab * (p6b_2 - noab - 1 + nvab * (p5b_
     &2 - noab - 1)))))
!
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h2b-1),int_mb(k_range+h3b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsd_t_singles_1',6,MA_E
     &RR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsd_t_singles_1',7
     &,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsd_t_singles_1',8
     &,MA_ERR)
      END IF
      END IF
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),6,4,3,5,2,1,1.0d0)
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),6,4,3,2,5,1,-1.0d0)
      END IF
      IF ((t_p4b .eq. p4b) .and. (t_p5b .eq. p5b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),6,4,3,2,1,5,1.0d0)
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),4,6,3,5,2,1,-1.0d0)
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),4,6,3,2,5,1,1.0d0)
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p4b) .and. (t_p6b .eq. p6b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),4,6,3,2,1,5,-1.0d0)
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h2b) .and. (t_h3b .eq. h
     &3b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),4,3,6,5,2,1,1.0d0)
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h1b) .and. (t_h3b .eq. h
     &3b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),4,3,6,2,5,1,-1.0d0)
      END IF
      IF ((t_p4b .eq. p5b) .and. (t_p5b .eq. p6b) .and. (t_p6b .eq. p4b)
     & .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h3b) .and. (t_h3b .eq. h
     &1b)) THEN
      CALL TCE_SORTACC_6(dbl_mb(k_c_sort),a_c,int_mb(k_range+h3b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1),int_
     &mb(k_range+h1b-1),int_mb(k_range+p4b-1),4,3,6,2,1,5,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsd_t_singles_1',9
     &,MA_ERR)
      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
