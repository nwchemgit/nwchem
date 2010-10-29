      SUBROUTINE cr_eomccsd_t_gr(d_t1,k_t1_offset,d_t2,k_t2_offset,
     2                        d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                        d_e,k_e_offset,
     4                        d_ex1,k_ex1_offset,size_ex1,
     5                        d_ex2,k_ex2_offset,size_ex2,
     7                        d_c2,k_c2_offset,size_c2,
     8                        energy1)
C
C     $Id$
C 
c 
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
c      integer size_d2
c *************
c      integer d_x1
c      integer k_x1_offset
c      integer d_x2
c      integer k_x2_offset
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
      integer d_c2
      integer k_c2_offset
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
      integer k_right1,l_right1  !r0*M3(T1,T2)
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
c  cr_ccsd_t_N
      integer d_i1_1,d_i1_2
      integer k_i1_offset_1,k_i1_offset_2
      integer l_i1_offset_1,l_i1_offset_2
c  cr_ccsd_t_E
      integer d_i1_4,k_i1_offset_4,l_i1_offset_4
c  cr_ccsd_t_N2
c      integer d_i2_1,d_i2_2,d_i2_3,d_i2_4,d_i2_5,d_i2_6
c      integer k_i2_offset_1,k_i2_offset_2,k_i2_offset_3
c      integer k_i2_offset_4,k_i2_offset_5,k_i2_offset_6
c      integer l_i2_offset_1,l_i2_offset_2,l_i2_offset_3
c      integer l_i2_offset_4,l_i2_offset_5,l_i2_offset_6
c  q3rexpt
c      integer d_i3_1,d_i3_2
c      integer k_i3_offset_1,k_i3_offset_2
c      integer l_i3_offset_1,l_i3_offset_2
c ----------------------------
      double precision energy1
      double precision factor
      double precision den1,num1
      character*255 filename
      external nxtask
c *** debug ***
      double precision xxx
c *************
c
c Getting R0 for the ground-state R0=1
c
      r0xx=1.0d0
      dr0xx = r0xx*r0xx
      lr0 = .true.
c *** debug ***
c      write(6,*)'gr step1'
c      call util_flush(6)
c *************
c 
c Calculating one- and two-body overlaps
c
      if(lr0) then !symmetry of the reference
         call tce_zero(d_ex1,size_ex1)
         call tce_zero(d_ex2,size_ex2)
         call c1_c1(d_t1,d_ex1,k_t1_offset,k_ex1_offset) 
         call reconcilefile(d_ex1,1)
c *** debug ***
c      write(6,*)'gr step1'
c      call util_flush(6)
c *************
         d1xxt = 0.0d0 
         call get_block(d_ex1,d1xxt,1,0)
         call t2t12(d_c2,d_t1,d_t2,k_c2_offset,k_t1_offset,k_t2_offset)
         call c2_c2(d_c2,d_ex2,k_c2_offset,k_ex2_offset)
         call reconcilefile(d_ex2,1)
         d2xxt = 0.0d0
         call get_block(d_ex2,d2xxt,1,0)
         d12xx = d1xxt*r0xx*r0xx+d2xxt*r0xx*r0xx
c *** debug ***
c      write(6,*)'gr step2'
c      call util_flush(6)
c *************
      end if
c
c     Caution! k_right & k_den are not even allocated yet
c     but they won't be used.
c --------------- initialization -----------------------------------
      if(lr0) then
      call cr_ccsd_t_N(dbl_mb(k_right1),d_f1,d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      end if
c
c *** debug ***
c      write(6,*)'gr step3'
c      call util_flush(6)
c *************
      if(lr0) then
      call cr_ccsd_t_E(dbl_mb(k_den0),d_i1_4,
     1  d_t1,d_t2,k_i1_offset_4,k_t1_offset,k_t2_offset,
     2  l_i1_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      end if 
c *** debug ***
c      write(6,*)'gr step4'
c      call util_flush(6)
c *************
c ------------------------------------------------------------------
c
c
c     Get the numerator
c
      num1 = 0.0d0
      den1 = 0.0d0
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
          if(lr0) then
            do i = 1, size
             dbl_mb(k_right1+i-1) = 0.0d0
             dbl_mb(k_den0+i-1) = 0.0d0
            enddo
          end if

c
c Moments are calculated here
c 
c
      if(lr0) then   
      call cr_ccsd_t_N(dbl_mb(k_right1),d_f1,d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
      end if
c
c Q3 exp(T1+T2)|Ref> calculated here
c
      if(lr0) then
      call cr_ccsd_t_E(dbl_mb(k_den0),d_i1_4,
     1  d_t1,d_t2,k_i1_offset_4,k_t1_offset,k_t2_offset,
     2  l_i1_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
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
c <(R3^{\dagger})*M3>
                  num1 = num1 + factor * 
     1            (r0xx*dbl_mb(k_right1+i-1))*
     2            (r0xx*dbl_mb(k_right1+i-1))
     3      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     5         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c <[R0*(T1T2+1/6T1^3)]^{\dagger}*M3>
                  num1 = num1 + factor * 
     1            (r0xx*dbl_mb(k_den0+i-1))*
     2            (r0xx*dbl_mb(k_right1+i-1))
c
                  den1 = den1 + factor * 
     1            (r0xx*dbl_mb(k_den0+i-1))*
     2            (r0xx*dbl_mb(k_right1+i-1))
     3      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     5         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                  den1 = den1 + factor * 
     1             (r0xx*dbl_mb(k_den0+i-1))*
     2             (r0xx*dbl_mb(k_den0+i-1))
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
      if(lr0) then
      call cr_ccsd_t_E(dbl_mb(k_den0),d_i1_4,
     1  d_t1,d_t2,k_i1_offset_4,k_t1_offset,k_t2_offset,
     2  l_i1_offset_4,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      end if
c
      if(lr0) then
      call cr_ccsd_t_N(dbl_mb(k_right1),d_f1,d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
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
      if (.not.ga_destroy(g_energy))
     1  call errquit('ccsd_t: GA problem',1,GA_ERR)
       energy1 = num1/(dr0xx+d12xx+den1)
c format statements
 1000    format('corresponding excitation energy = 0')
c
      return
      end
c
