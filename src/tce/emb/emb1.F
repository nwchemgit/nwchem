      SUBROUTINE emb1(d_x1,k_x1_offset,d_x2,k_x2_offset,
     1                        d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     6                        d_c1,k_c1_offset,size_c1,
     7                        d_c2,k_c2_offset,size_c2,
     8                        excit,energy1,energy2,size_t1,size_x1)
C
C     $Id$
C 
c from this point on we assume that corresponding one- and two-
c body components of R(d_x1,d_x2) and T(d_t1,d_t2) operators,
c and corresponding excitation energy are available.  
c objects c1 and c2 are created in tce_energy right before 
c calling this procedure. The same structure of x1 and c1 and
c x2 and c2 is assumed (irrep_c = irrep_x) 
c
c 
c Simplified version of bra-states, full form of moments
c
c

      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1
      integer k_t1_offset
      integer d_t2
      integer k_t2_offset
c *** sizes ***
      integer size_c1,size_c2
      integer size_ex1,size_ex2
      integer size_d2
c *************
      integer d_x1
      integer k_x1_offset
      integer d_x2
      integer k_x2_offset
      double precision r0xx  ! r0 
      double precision dr0xx ! r0*r0
      double precision d1xx  ! <(singles)+singles>
      double precision d2xx  ! <(doubles)+doubles>
      double precision d12xx ! <(singles+doubles)+(singles+doubles)>
      double precision d1xxt,d2xxt
      double precision d1xxr,d2xxr 
      double precision d1xxtr,d2xxtr 
      double precision excit ! eomsd excitation energy
      logical lr0            ! (true) r0*M3-calculated 
      integer d_ex1
      integer k_ex1_offset
      integer d_ex2
      integer k_ex2_offset 
      integer d_c1
      integer k_c1_offset
      integer d_c2
      integer k_c2_offset
      integer d_d2
      integer l_d2_offset
      integer k_d2_offset
c ***
      integer d_f1
      integer k_f1_offset
      integer d_v2
      integer k_v2_offset
      integer d_e
      integer k_e_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_den,l_den
      integer k_right1,l_right1  !r0*M3(T1,T2)
      integer k_right2,l_right2  !M3(T1,T2,X1,X2)
      integer k_right3,l_right3
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer l_x1_local,k_x1_local
      integer size_t1,size_x1
c ---------------------------------------
c *** cr_ccsd_t_E ***
      integer k_den0,l_den0
c *******************
      integer size,i
      integer g_energy
      integer nxtask
      integer next
      integer nprocs
      integer count
c --- new intermediates ---
c  cr_ccsd_t_N1 or cr_ccsd_t_N
      integer d_i1_1,d_i1_2,d_i1_3
      integer k_i1_offset_1,k_i1_offset_2,k_i1_offset_3
      integer l_i1_offset_1,l_i1_offset_2,l_i1_offset_3
c  cr_ccsd_t_E
      integer d_i1_4,k_i1_offset_4,l_i1_offset_4
c  cr_ccsd_t_N2
      integer d_i2_1,d_i2_2,d_i2_3,d_i2_4,d_i2_5,d_i2_6
      integer k_i2_offset_1,k_i2_offset_2,k_i2_offset_3
      integer k_i2_offset_4,k_i2_offset_5,k_i2_offset_6
      integer l_i2_offset_1,l_i2_offset_2,l_i2_offset_3
      integer l_i2_offset_4,l_i2_offset_5,l_i2_offset_6
c  q3rexpt
      integer d_i3_1,d_i3_2
      integer k_i3_offset_1,k_i3_offset_2
      integer l_i3_offset_1,l_i3_offset_2
c  c2vt1_e
      integer d_ui1_1,d_ui1_2
      integer k_ui1_offset_1,k_ui1_offset_2
      integer l_ui1_offset_1,l_ui1_offset_2
c  c2f_e
      integer d_zi1_1
      integer k_zi1_offset_1
      integer l_zi1_offset_1 
c 
      double precision bdenex,bdengr
      double precision coeffx
c ----------------------------
      integer k_c2vs,k_c2ve,k_c2fs,k_c2fe
      integer l_c2vs,l_c2ve,l_c2fs,l_c2fe
      integer l_c1v,k_c1v
c ----------------------------
      double precision energy1,energy2
      double precision factor
      double precision den1,num1
      double precision den2,num2
      character*255 filename
      external nxtask
c *** debug ***
      double precision xxx
c *************
c
c Getting R0
      call  nr0(d_f1,d_e,d_t1,d_v2,d_x1,d_x2,k_f1_offset,
     &k_e_offset,k_t1_offset,k_v2_offset,k_x1_offset,k_x2_offset) 
      call reconcilefile(d_e,1)
      call get_block(d_e,r0xx,1,0)
      if(dabs(excit).gt.1.0d-7) then
         r0xx = r0xx/excit
      else
         write(6,1000)
      end if
      dr0xx = r0xx*r0xx
      lr0 = .true.
      if(dabs(r0xx).lt.1.0d-7) lr0 = .false.
c 
c Now on ga with handle d_e we store corresponding R0 value 
c
c
c
c Calculating one- and two-body overlaps
c
      if(lr0) then !symmetry of the reference
         call tce_zero(d_ex1,size_ex1)
         call tce_zero(d_ex2,size_ex2)
         call c1_c1(d_t1,d_ex1,k_t1_offset,k_ex1_offset) 
         call reconcilefile(d_ex1,1)
         d1xxt = 0.0d0 
         call get_block(d_ex1,d1xxt,1,0)
         call t2t12(d_c2,d_t1,d_t2,k_c2_offset,k_t1_offset,k_t2_offset)
         call c2_c2(d_c2,d_ex2,k_c2_offset,k_ex2_offset)
         call reconcilefile(d_ex2,1)
         d2xxt = 0.0d0
         call get_block(d_ex2,d2xxt,1,0)
c forming vector R1T1+R2, d_d2's irrep corresponds to irrep_x
c (in this case this is fully symmetric situation)
c on d_c2 we have (T2+1/2T1*T1)|Phi\rangle
            irrep_d=irrep_c
            call tce_zero(d_ex1,size_ex1)
            call tce_x2_offset(l_d2_offset,k_d2_offset,size_d2)
            call tce_filename('d2',filename)
            call createfile(filename,d_d2,size_d2)
            call tce_zero(d_d2,size_d2)
         call c2excit2(d_d2,d_t1,d_x1,d_x2,k_d2_offset,k_t1_offset,k_
     &   x1_offset,k_x2_offset)
c <C2+ D2>  here
      call c2_d2(d_c2,d_d2,d_ex1,k_c2_offset,k_d2_offset,k_ex1_offset)
         call reconcilefile(d_ex1,1)
         d2xxtr = 0.0d0
         call get_block(d_ex1,d2xxtr,1,0)
            call deletefile(d_d2)
            if (.not.ma_pop_stack(l_d2_offset))
     1        call errquit("tce_energy: MA problem",36,MA_ERR)
ccxxxx       call c2excit(d_c2,d_t1,d_t2,d_e,d_x1,d_x2,k_c2_offset,
ccxxxx     & k_t1_offset,k_t2_offset,k_e_offset,k_x1_offset,k_x2_offset)
ccxxxx       call c2_c2(d_c2,d_ex2,k_c2_offset,k_ex2_offset)
ccxxxx       call reconcilefile(d_ex2,1)
ccxxxx       d2xx = 0.0d0
ccxxxx       call get_block(d_ex2,d2xx,1,0)
c
         call tce_zero(d_ex1,size_ex1)
         call tce_zero(d_ex2,size_ex2)
         call x1_t1(d_ex1,d_t1,d_x1,k_ex1_offset,k_t1_offset,
     &   k_x1_offset)
         call reconcilefile(d_ex1,1)
         d1xxtr=0.0d0
         call get_block(d_ex1,d1xxtr,1,0)
c
         call tce_zero(d_ex1,size_ex1)
         call tce_zero(d_ex2,size_ex2)
         call tce_zero(d_c1,size_c1)
         call tce_zero(d_c2,size_c2)
         call c1excit2(d_c1,d_x1,k_c1_offset,k_x1_offset)
         call c2excit2(d_c2,d_t1,d_x1,d_x2,k_c2_offset,k_t1_offset,k_
     &   x1_offset,k_x2_offset)
         call c1_c1(d_c1,d_ex1,k_c1_offset,k_ex1_offset)
         call reconcilefile(d_ex1,1)
         d1xxr = 0.0d0 
         call get_block(d_ex1,d1xxr,1,0)
         call c2_c2(d_c2,d_ex2,k_c2_offset,k_ex2_offset)
         call reconcilefile(d_ex2,1)
         d2xxr = 0.0d0
         call get_block(d_ex2,d2xxr,1,0)
         d12xx = d1xxt*r0xx*r0xx+d2xxt*r0xx*r0xx+d1xxr+d2xxr+
     &           2.0d0*d1xxtr*r0xx+2.0d0*d2xxtr*r0xx

      else         !symmetry different form the symmetry of the reference
c        <(R1+R0*T1)^{\dagger} (R1+R0*T1)>
         call c1excit2(d_c1,d_x1,k_c1_offset,k_x1_offset)
         call c2excit2(d_c2,d_t1,d_x1,d_x2,k_c2_offset,k_t1_offset,k_
     &   x1_offset,k_x2_offset)
         call c1_c1(d_c1,d_ex1,k_c1_offset,k_ex1_offset)
         call reconcilefile(d_ex1,1)
         d1xx = 0.0d0 
         call get_block(d_ex1,d1xx,1,0)
c        <(R2+R1T1)^{\dagger} (R2+R1T1)>
         call c2_c2(d_c2,d_ex2,k_c2_offset,k_ex2_offset)
         call reconcilefile(d_ex2,1)
         d2xx = 0.0d0
         call get_block(d_ex2,d2xx,1,0)
         d12xx = d1xx+d2xx
      end if
c
c
c ---- initializing c1 and c2
c   
      if(lr0) then
         coeffx = 1.0d0
         call tce_zero(d_c1,size_c1)
         call tce_zero(d_c2,size_c2)
         call copyfile(d_t1,d_c1,size_c1)
         call dscalfile(r0xx,d_c1,size_c1)
         call daxpyfile(1,coeffx,d_x1,d_c1,size_c1)
         call t2t12(d_c2,d_t1,d_t2,k_c2_offset,k_t1_offset,k_t2_offset)
         call dscalfile(r0xx,d_c2,size_c2)
         call c2excit2(d_c2,d_t1,d_x1,d_x2,k_c2_offset,k_t1_offset,k_
     &   x1_offset,k_x2_offset)
      else
         call tce_zero(d_c1,size_c1)
         call tce_zero(d_c2,size_c2)
         call copyfile(d_x1,d_c1,size_c1)
         call c2excit2(d_c2,d_t1,d_x1,d_x2,k_c2_offset,k_t1_offset,k_
     &   x1_offset,k_x2_offset)
      end if
c - T1/X1 LOCALIZATION ----------
c    opening l_t1_local and l_x1_local
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        if (.not.MA_PUSH_GET(mt_dbl,size_x1,'x1_local',
     1      l_x1_local,k_x1_local))
     1      call errquit('x1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
        call ma_zero(dbl_mb(k_x1_local),size_x1)
c    copy d_t1 ==> l_t1_local
c    copy x1(ivec) ==> l_x1_local
cc        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
cc        call ga_get(d_x1,1,size_x1,1,1,dbl_mb(k_x1_local),1)
       call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
       call get_block(d_x1,dbl_mb(k_x1_local),size_x1,0)
c -------------------------------
c
c
c     Caution! k_right & k_den are not even allocated yet
c     but they won't be used.
c --------------- initialization -----------------------------------
      if(lr0) then
      call cr_ccsd_t_N(dbl_mb(k_right1),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      end if
c
c      call cr_eomccsd_t_N2(dbl_mb(k_right2),d_f1,d_i2_1,d_i2_2,d_i2_3,
c     &d_i2_4,d_i2_5,d_i2_6,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offset,
c     &k_i2_offset_1,k_i2_offset_2,k_i2_offset_3,k_i2_offset_4,
c     &k_i2_offset_5,k_i2_offset_6,k_t1_offset,k_t2_offset,k_v2_offset,
c     &k_x1_offset,k_x2_offset,l_i2_offset_1,l_i2_offset_2,l_i2_offset_3,
c     &l_i2_offset_4,l_i2_offset_5,l_i2_offset_6,t_h1b,t_h2b,t_h3b,t_p4b,
c     &t_p5b,t_p6b,1)
      call creomsd_t_n2_mem(dbl_mb(k_right2),d_f1,d_i2_1,d_i2_2,
     &d_i2_3,d_i2_4,k_t1_local,d_t2,d_v2,k_x1_local,
     &d_x2,k_f1_offset,k_i2_offset_1,
     &k_i2_offset_2,k_i2_offset_3,k_i2_offset_4,k_t1_offset,k_t2_offset,
     &k_v2_offset,k_x1_offset,k_x2_offset,l_i2_offset_1,l_i2_offset_2,
     &l_i2_offset_3,l_i2_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
      if(lr0) then !--
      call cr_ccsd_t_E(dbl_mb(k_den0),d_i1_4,
     1  k_t1_local,d_t2,k_i1_offset_4,k_t1_offset,k_t2_offset,
     2  l_i1_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      call q3rexpt2(dbl_mb(k_den),d_i3_1,k_t1_local,d_t2,
     &k_x1_local,d_x2,
     &k_i3_offset_1,k_t1_offset,k_t2_offset,k_x1_offset,k_x2_offset,
     &l_i3_offset_1,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      else !--
      call q3rexpt2(dbl_mb(k_den),d_i3_1,k_t1_local,d_t2,
     &k_x1_local,d_x2,
     &k_i3_offset_1,k_t1_offset,k_t2_offset,k_x1_offset,k_x2_offset,
     &l_i3_offset_1,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      end if !--
c
      call c2vt1_s(dbl_mb(k_c2vs),d_v2,d_c2,k_v2_offset,k_c2_offset,
     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
ccc      call c2vt1_e(dbl_mb(k_c2ve),d_ui1_1,d_ui1_2,d_t1,d_v2,d_c2,
ccc     &k_ui1_offset_1,k_ui1_offset_2,k_t1_offset,k_v2_offset,k_c2_offset,
ccc     &l_ui1_offset_1,l_ui1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
ccc     &1)
c
ccc      call c2f_s(dbl_mb(k_c2fs),d_f1,d_c2,k_f1_offset,k_c2_offset,
ccc     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
ccc      call c2f_e(dbl_mb(k_c2fe),d_zi1_1,d_t1,d_v2,d_c2,k_zi1_offset_1,
ccc     &k_t1_offset,k_v2_offset,k_c2_offset,l_zi1_offset_1,
ccc     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
      call c1v(dbl_mb(k_c1v),d_v2,d_c1,k_v2_offset,k_c1_offset,
     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c 
c ------------------------------------------------------------------
c k_c2vs,k_c2ve,k_c2fs,k_c2fe
c
c     Get the numerator
c
      num1 = 0.0d0
      den1 = 0.0d0
c
      num2 = 0.0d0
      den2 = 0.0d0
c
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('ccsd_t: GA problem',0,GA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
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
     5               int_mb(k_sym+t_h3b-1)))))).eq.irrep_x) then
c
c Symmetry control (above)
c 
            size = int_mb(k_range+t_p4b-1)
     1           * int_mb(k_range+t_p5b-1)
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
          if(lr0) then 
            if (.not.MA_PUSH_GET(mt_dbl,size,'moment 1,3',
     1        l_right1,k_right1)) call errquit('eomccsd_t',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'t1t2_t13',
     1        l_den0,k_den0)) call errquit('eomccsd_t',3,MA_ERR)
          end if
            if (.not.MA_PUSH_GET(mt_dbl,size,'moment 2,3',
     1        l_right2,k_right2)) call errquit('eomccsd_t',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'denominator',
     1        l_den,k_den)) call errquit('ccsd_t',3,MA_ERR)
c
            if (.not.MA_PUSH_GET(mt_dbl,size,'c2vs',
     1        l_c2vs,k_c2vs)) call errquit('c2vs',3,MA_ERR)
ccc            if (.not.MA_PUSH_GET(mt_dbl,size,'c2ve',
ccc     1        l_c2ve,k_c2ve)) call errquit('c2ve',3,MA_ERR)
ccc            if (.not.MA_PUSH_GET(mt_dbl,size,'c2fs',
ccc     1        l_c2fs,k_c2fs)) call errquit('c2fs',3,MA_ERR)
ccc            if (.not.MA_PUSH_GET(mt_dbl,size,'c2fe',
ccc     1        l_c2fe,k_c2fe)) call errquit('c2fe',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'c1v',
     1        l_c1v,k_c1v)) call errquit('c2fe',3,MA_ERR)
c
          if(lr0) then
            do i = 1, size
             dbl_mb(k_right1+i-1) = 0.0d0
             dbl_mb(k_den0+i-1) = 0.0d0
            enddo
          end if
            do i = 1, size
cccx             dbl_mb(k_right3+i-1) = 0.0d0
             dbl_mb(k_right2+i-1) = 0.0d0
             dbl_mb(k_den+i-1) = 0.0d0
             dbl_mb(k_c2vs+i-1) = 0.0d0
ccc             dbl_mb(k_c2ve+i-1) = 0.0d0
ccc             dbl_mb(k_c2fs+i-1) = 0.0d0
ccc             dbl_mb(k_c2fe+i-1) = 0.0d0 
             dbl_mb(k_c1v+i-1) = 0.0d0
            enddo
c
c Moments are calculated here
c 
c
      if(lr0) then   
      call cr_ccsd_t_N(dbl_mb(k_right1),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
      end if
c
c      call cr_eomccsd_t_N2(dbl_mb(k_right2),d_f1,d_i2_1,d_i2_2,d_i2_3,
c     &d_i2_4,d_i2_5,d_i2_6,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offset,
c     &k_i2_offset_1,k_i2_offset_2,k_i2_offset_3,k_i2_offset_4,
c     &k_i2_offset_5,k_i2_offset_6,k_t1_offset,k_t2_offset,k_v2_offset,
c     &k_x1_offset,k_x2_offset,l_i2_offset_1,l_i2_offset_2,l_i2_offset_3,
c     &l_i2_offset_4,l_i2_offset_5,l_i2_offset_6,t_h1b,t_h2b,t_h3b,t_p4b,
c     &t_p5b,t_p6b,2)
      call creomsd_t_n2_mem(dbl_mb(k_right2),d_f1,d_i2_1,d_i2_2,
     &d_i2_3,d_i2_4,k_t1_local,d_t2,d_v2,
     &k_x1_local,d_x2,k_f1_offset,k_i2_offset_1,
     &k_i2_offset_2,k_i2_offset_3,k_i2_offset_4,k_t1_offset,k_t2_offset,
     &k_v2_offset,k_x1_offset,k_x2_offset,l_i2_offset_1,l_i2_offset_2,
     &l_i2_offset_3,l_i2_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
      call c2vt1_s(dbl_mb(k_c2vs),d_v2,d_c2,k_v2_offset,k_c2_offset,
     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
ccc      call c2vt1_e(dbl_mb(k_c2ve),d_ui1_1,d_ui1_2,d_t1,d_v2,d_c2,
ccc     &k_ui1_offset_1,k_ui1_offset_2,k_t1_offset,k_v2_offset,k_c2_offset,
ccc     &l_ui1_offset_1,l_ui1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
ccc     &2)
c
ccc      call c2f_s(dbl_mb(k_c2fs),d_f1,d_c2,k_f1_offset,k_c2_offset,
ccc     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
ccc      call c2f_e(dbl_mb(k_c2fe),d_zi1_1,d_t1,d_v2,d_c2,k_zi1_offset_1,
ccc     &k_t1_offset,k_v2_offset,k_c2_offset,l_zi1_offset_1,
ccc     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
      call c1v(dbl_mb(k_c1v),d_v2,d_c1,k_v2_offset,k_c1_offset,
     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
c
c
c Q3(R0+R1+R2)exp(T1+T2)|Ref> calculated here
c
      if(lr0) then
      call cr_ccsd_t_E(dbl_mb(k_den0),d_i1_4,
     1  k_t1_local,d_t2,k_i1_offset_4,k_t1_offset,k_t2_offset,
     2  l_i1_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
      call q3rexpt2(dbl_mb(k_den),d_i3_1,k_t1_local,
     &d_t2,k_x1_local,d_x2,
     &k_i3_offset_1,k_t1_offset,k_t2_offset,k_x1_offset,k_x2_offset,
     &l_i3_offset_1,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
      else
      call q3rexpt2(dbl_mb(k_den),d_i3_1,k_t1_local,d_t2,
     &k_x1_local,d_x2,
     &k_i3_offset_1,k_t1_offset,k_t2_offset,k_x1_offset,k_x2_offset,
     &l_i3_offset_1,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
      end if
c
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
c
c
c
      if(lr0) then  !excited-state symm. = symm. of ref.
            i = 0
            do t_p4 = 1, int_mb(k_range+t_p4b-1)
             do t_p5 = 1, int_mb(k_range+t_p5b-1)
              do t_p6 = 1, int_mb(k_range+t_p6b-1)
               do t_h1 = 1, int_mb(k_range+t_h1b-1)
                do t_h2 = 1, int_mb(k_range+t_h2b-1)
                 do t_h3 = 1, int_mb(k_range+t_h3b-1)
                  i = i + 1
          bdenex=(-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     5         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1)
     9         +excit)
          bdengr=(-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     5         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c <(C2V1-bar + C1V2-bar){\dagger}*R3>
c                  num1 = num1 + factor * 
c     1            (r0xx*dbl_mb(k_right1+i-1)+dbl_mb(k_right2+i-1))*
c     2            (dbl_mb(k_c2vs+i-1)+dbl_mb(k_c2ve+i-1)+
c     2             dbl_mb(k_c2fs+i-1)+dbl_mb(k_c2fe+i-1)+
c     2             dbl_mb(k_c1v+i-1))
c     3      / bdenex
                  num1 = num1 + factor *
     1            (r0xx*dbl_mb(k_right1+i-1)+dbl_mb(k_right2+i-1))*
     2            (dbl_mb(k_c2vs+i-1)+dbl_mb(k_c1v+i-1))
     3      / bdenex
c Multiplied by 2
                  den1 = den1 + factor * 2.0 *
     1            (r0xx*dbl_mb(k_den0+i-1)+dbl_mb(k_den+i-1))*
     2            (r0xx*dbl_mb(k_right1+i-1)+dbl_mb(k_right2+i-1))
     3      / bdenex
c ----
                  den1 = den1 + factor * 
     1             (r0xx*dbl_mb(k_den0+i-1)+dbl_mb(k_den+i-1))*
     2             (r0xx*dbl_mb(k_den0+i-1)+dbl_mb(k_den+i-1))
c ---- Multiplied by 2
                  den1 = den1 + factor * 2.0 *
     1            (r0xx*dbl_mb(k_den0+i-1)+dbl_mb(k_den+i-1))*
     2            (r0xx*dbl_mb(k_right1+i-1))
     3      / bdengr
c ----
                  den1 = den1 + factor * 
     1            (r0xx*dbl_mb(k_right1+i-1)*
     2             r0xx*dbl_mb(k_right1+i-1))
     3      / (bdengr*bdengr)
c ----
                  den1 = den1 + factor *
     2            (r0xx*dbl_mb(k_right1+i-1)+dbl_mb(k_right2+i-1))*
     2            (r0xx*dbl_mb(k_right1+i-1)+dbl_mb(k_right2+i-1))
     3      / (bdenex*bdenex)
c ----
                  den1 = den1 + factor * 2.0 *
     2            (r0xx*dbl_mb(k_right1+i-1)+dbl_mb(k_right2+i-1))*
     3            (r0xx*dbl_mb(k_right1+i-1))
     4      / (bdenex*bdengr)
c
c
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo 
      else       !excited-state symm. different from the symm. of ref.
            i = 0
            do t_p4 = 1, int_mb(k_range+t_p4b-1)
             do t_p5 = 1, int_mb(k_range+t_p5b-1)
              do t_p6 = 1, int_mb(k_range+t_p6b-1)
               do t_h1 = 1, int_mb(k_range+t_h1b-1)
                do t_h2 = 1, int_mb(k_range+t_h2b-1)
                 do t_h3 = 1, int_mb(k_range+t_h3b-1)
                  i = i + 1
          bdenex=(-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     5         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1)
     9         +excit)
          bdengr=(-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     5         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c---
ccc                  num1 = num1 + factor * 
ccc     1            dbl_mb(k_right2+i-1)*
ccc     2            (dbl_mb(k_c2vs+i-1)+dbl_mb(k_c2ve+i-1)+
ccc     2             dbl_mb(k_c2fs+i-1)+dbl_mb(k_c2fe+i-1)+
ccc     2             dbl_mb(k_c1v+i-1))
ccc     2      / bdenex 
                  num1 = num1 + factor *
     1            dbl_mb(k_right2+i-1)*
     2            (dbl_mb(k_c2vs+i-1)+dbl_mb(k_c1v+i-1))
     2      / bdenex
c
c ----
                  den1 = den1 + factor * 2.0d0 *
     1            dbl_mb(k_den+i-1)*dbl_mb(k_right2+i-1)
     2      / bdenex 
c ----
                  den1 = den1 + factor * dbl_mb(k_den+i-1)*
     1                          dbl_mb(k_den+i-1)
c ----
                  den1 = den1 + factor * dbl_mb(k_right2+i-1) *
     1                   dbl_mb(k_right2+i-1)/(bdenex*bdenex) 
c
c
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo 
      end if
c
c
c
            if (.not.MA_POP_STACK(l_c1v))
     1        call errquit('c1v',6,MA_ERR)
ccc            if (.not.MA_POP_STACK(l_c2fe))
ccc     1        call errquit('c2fe',6,MA_ERR)
ccc            if (.not.MA_POP_STACK(l_c2fs))
ccc     1        call errquit('c2fs',6,MA_ERR)
ccc            if (.not.MA_POP_STACK(l_c2ve))
ccc     1        call errquit('c2ve',6,MA_ERR)
            if (.not.MA_POP_STACK(l_c2vs))
     1        call errquit('c2vs',6,MA_ERR)
c
            if (.not.MA_POP_STACK(l_den)) 
     1        call errquit('eomccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_right2)) 
     1        call errquit('eomccsd_t',6,MA_ERR)
          if(lr0) then
            if (.not.MA_POP_STACK(l_den0)) 
     1        call errquit('eomccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_right1)) 
     1        call errquit('eomccsd_t',6,MA_ERR)
          end if
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
c --- toggle = 3 ---
c
      call c1v(dbl_mb(k_c1v),d_v2,d_c1,k_v2_offset,k_c1_offset,
     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
c
ccc      call c2f_e(dbl_mb(k_c2fe),d_zi1_1,d_t1,d_v2,d_c2,k_zi1_offset_1,
ccc     &k_t1_offset,k_v2_offset,k_c2_offset,l_zi1_offset_1,
ccc     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
c
ccc      call c2f_s(dbl_mb(k_c2fs),d_f1,d_c2,k_f1_offset,k_c2_offset,
ccc     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
c
ccc      call c2vt1_e(dbl_mb(k_c2ve),d_ui1_1,d_ui1_2,d_t1,d_v2,d_c2,
ccc     &k_ui1_offset_1,k_ui1_offset_2,k_t1_offset,k_v2_offset,k_c2_offset,
ccc     &l_ui1_offset_1,l_ui1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
ccc     &3)
c
      call c2vt1_s(dbl_mb(k_c2vs),d_v2,d_c2,k_v2_offset,k_c2_offset,
     &t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
c
      if(lr0) then
      call q3rexpt2(dbl_mb(k_den),d_i3_1,k_t1_local,d_t2,
     &k_x1_local,d_x2,
     &k_i3_offset_1,k_t1_offset,k_t2_offset,k_x1_offset,k_x2_offset,
     &l_i3_offset_1,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      call cr_ccsd_t_E(dbl_mb(k_den0),d_i1_4,
     1  k_t1_local,d_t2,k_i1_offset_4,k_t1_offset,k_t2_offset,
     2  l_i1_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      else
      call q3rexpt2(dbl_mb(k_den),d_i3_1,k_t1_local,d_t2,
     &k_x1_local,d_x2,
     &k_i3_offset_1,k_t1_offset,k_t2_offset,k_x1_offset,k_x2_offset,
     &l_i3_offset_1,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      end if
c
c      call cr_eomccsd_t_N2(dbl_mb(k_right2),d_f1,d_i2_1,d_i2_2,d_i2_3,
c     &d_i2_4,d_i2_5,d_i2_6,d_t1,d_t2,d_v2,d_x1,d_x2,k_f1_offset,
c     &k_i2_offset_1,k_i2_offset_2,k_i2_offset_3,k_i2_offset_4,
c     &k_i2_offset_5,k_i2_offset_6,k_t1_offset,k_t2_offset,k_v2_offset,
c     &k_x1_offset,k_x2_offset,l_i2_offset_1,l_i2_offset_2,l_i2_offset_3,
c     &l_i2_offset_4,l_i2_offset_5,l_i2_offset_6,t_h1b,t_h2b,t_h3b,t_p4b,
c     &t_p5b,t_p6b,3)
      call creomsd_t_n2_mem(dbl_mb(k_right2),d_f1,d_i2_1,d_i2_2,
     &d_i2_3,d_i2_4,k_t1_local,d_t2,d_v2,
     &k_x1_local,d_x2,k_f1_offset,k_i2_offset_1,
     &k_i2_offset_2,k_i2_offset_3,k_i2_offset_4,k_t1_offset,k_t2_offset,
     &k_v2_offset,k_x1_offset,k_x2_offset,l_i2_offset_1,l_i2_offset_2,
     &l_i2_offset_3,l_i2_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
c
      if(lr0) then
      call cr_ccsd_t_N(dbl_mb(k_right1),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      end if
c ------------------
c
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,num1,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,num1,1)
c
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den1,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den1,1)
c
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,num2,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,num2,1)
c
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den2,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den2,1)
c
cccx      call ga_get(g_energy,1,1,1,1,den2,1)
      if (.not.ga_destroy(g_energy))
     1  call errquit('ccsd_t: GA problem',1,GA_ERR)
c      den1 = den1 + 
c      den2 = den2 +
       energy1 = num1/(dr0xx+d12xx+den1)
cccx       energy2 = num2/(dr0xx+d12xx+den2)
c *** debug *********************************
c       write(6,*)'## from cr-eom ----- num2',num2
c       write(6,*)'## from cr-eom ----- den2',dr0xx+d12xx+den2
c      write(6,*)'from cr_eomsd_t - ALL DONE'
c      write(6,*)'dr0xx',dr0xx
c      write(6,*)'d12xx',d12xx
c      write(6,*)'den1',den1
c      write(6,*)'num1--- ',num1
c      write(6,*)'den1--- ',den1
c      write(6,*)'den --- ',dr0xx+d12xx+den1
c      call util_flush(6)
c *******************************************
c format statements
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_x1_local))
     &      call errquit('l_x1_local',4,MA_ERR)
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
 1000    format('corresponding excitation energy = 0')
c
      return
      end
c
