      SUBROUTINE tce_dens_hh_print(d_a,k_a_offset,r1_thresh,irrep)
C     dens_hh ( h6 h5 ) - fortran77 convention
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER h5b,h5
      INTEGER h6b,h6
      integer l_b,k_b,size_b
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_h5b,irrep_h6b
      character*5 spinnames(2),spin_h5b,spin_h6b
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0) !nodezero
      if (nodezero) then
         write(LuOut,9000) '1PDM h-h block largest contributions'
      DO h5b = 1,noab
      DO h6b = 1,noab
      IF (int_mb(k_spin+h5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+h5b-1),int_mb(k_sym+h6b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
c
      size_b = int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_b,'d1hh',l_b,k_b))
     &       call errquit('tce_print_hh_1dm: MA problem',0,MA_ERR)
       call get_hash_block(d_a,dbl_mb(k_b),size_b,
     &       int_mb(k_a_offset),h6b-1+noab*(h5b-1))
      i = 0
      do h5=1, int_mb(k_range+h5b-1)
      do h6=1, int_mb(k_range+h6b-1)
        i=i+1
                    if (dabs(dbl_mb(k_b+i-1)).gt.r1_thresh) then
                      call sym_irrepname(geom,int_mb(k_sym+h5b-1)+1,
     1                                   irrep_h5b)
                      call sym_irrepname(geom,int_mb(k_sym+h6b-1)+1,
     1                                   irrep_h6b)
                      spin_h5b = spinnames(int_mb(k_spin+h5b-1))
                      spin_h6b = spinnames(int_mb(k_spin+h6b-1))
                      if (nodezero) write(LuOut,9010)
     1                (int_mb(k_mo_index+int_mb(k_offset+h5b-1)+h5-1)+
     1                 2-int_mb(k_spin+h5b-1))/2,
     1                irrep_h5b,spin_h5b,
     2                (int_mb(k_mo_index+int_mb(k_offset+h6b-1)+h6-1)+
     1                 2-int_mb(k_spin+h6b-1))/2,
     2                irrep_h6b,spin_h6b,
     2                dbl_mb(k_b+i-1)
                    endif       
      enddo
      enddo
c
      if (.not.ma_pop_stack(l_b))
     1   call errquit('1DM_hh:MA_problem',1,MA_ERR)
c
      END IF
      END IF
      END IF
      END DO
      END DO
c
      write(LuOut,9001)
c
      end if  !nodezero
c
 9000 format(/,1x,a)
 9001 format(/)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
c
      return 
      end
c
      SUBROUTINE tce_dens_hp_print(d_a,k_a_offset,r1_thresh,irrep)
C     dens_pp ( p6 h5 ) - fortran77 convention
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER h5b,h5
      INTEGER p6b,p6
      integer l_b,k_b,size_b
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_h5b,irrep_p6b
      character*5 spinnames(2),spin_h5b,spin_p6b
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0) !nodezero
      if (nodezero) then
         write(LuOut,9000) '1PDM h-p block largest contributions'
c
      DO h5b = 1,noab
      DO p6b = noab+1,noab+nvab
      IF (int_mb(k_spin+h5b-1) .eq. int_mb(k_spin+p6b-1)) THEN
      IF (ieor(int_mb(k_sym+h5b-1),int_mb(k_sym+p6b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+h5b-1)+int_mb(k_spin+p6b-1
     &).ne.4)) THEN
c
      size_b = int_mb(k_range+h5b-1) * int_mb(k_range+p6b-1)
       if (.not.ma_push_get(mt_dbl,size_b,'d1hh',l_b,k_b))
     &       call errquit('tce_print_pp_1dm: MA problem',0,MA_ERR)
       call get_hash_block(d_a,dbl_mb(k_b),size_b,
     &       int_mb(k_a_offset),(p6b-noab-1)+nvab*(h5b-1))
      i = 0
      do h5=1, int_mb(k_range+h5b-1)
      do p6=1, int_mb(k_range+p6b-1)
        i=i+1
                    if (dabs(dbl_mb(k_b+i-1)).gt.r1_thresh) then
                      call sym_irrepname(geom,int_mb(k_sym+h5b-1)+1,
     1                                   irrep_h5b)
                      call sym_irrepname(geom,int_mb(k_sym+p6b-1)+1,
     1                                   irrep_p6b)
                      spin_h5b = spinnames(int_mb(k_spin+h5b-1))
                      spin_p6b = spinnames(int_mb(k_spin+p6b-1))
                      if (nodezero) write(LuOut,9010)
     1                (int_mb(k_mo_index+int_mb(k_offset+h5b-1)+h5-1)+
     1                 2-int_mb(k_spin+h5b-1))/2,
     1                irrep_h5b,spin_h5b,
     2                (int_mb(k_mo_index+int_mb(k_offset+p6b-1)+p6-1)+
     1                 2-int_mb(k_spin+p6b-1))/2,
     2                irrep_p6b,spin_p6b,
     2                dbl_mb(k_b+i-1)
                    endif
      enddo
      enddo
c
      if (.not.ma_pop_stack(l_b))
     1   call errquit('1DM_hh:MA_problem',1,MA_ERR)
c
      END IF
      END IF
      END IF
      END DO
      END DO
c
      write(LuOut,9001)
c
      end if  !nodezero
c
 9000 format(/,1x,a)
 9001 format(/)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
c
      return
      end
c
      SUBROUTINE tce_dens_ph_print(d_a,k_a_offset,r1_thresh,irrep)
C     dens_hh ( h6 p5 ) - fortran77 convention
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER p5b,p5
      INTEGER h6b,h6
      integer l_b,k_b,size_b
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_p5b,irrep_h6b
      character*5 spinnames(2),spin_p5b,spin_h6b
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0) !nodezero
      if (nodezero) then
         write(LuOut,9000) '1PDM p-h block largest contributions'
c
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+h6b-1
     &).ne.4)) THEN
c
      size_b = int_mb(k_range+p5b-1) * int_mb(k_range+h6b-1)
       if (.not.ma_push_get(mt_dbl,size_b,'d1hh',l_b,k_b))
     &       call errquit('tce_print_hh_1dm: MA problem',0,MA_ERR)
       call get_hash_block(d_a,dbl_mb(k_b),size_b,
     &       int_mb(k_a_offset),h6b-1+noab*(p5b-noab-1))
      i = 0
      do p5=1, int_mb(k_range+p5b-1)
      do h6=1, int_mb(k_range+h6b-1)
        i=i+1
                    if (dabs(dbl_mb(k_b+i-1)).gt.r1_thresh) then
                      call sym_irrepname(geom,int_mb(k_sym+p5b-1)+1,
     1                                   irrep_p5b)
                      call sym_irrepname(geom,int_mb(k_sym+h6b-1)+1,
     1                                   irrep_h6b)
                      spin_p5b = spinnames(int_mb(k_spin+p5b-1))
                      spin_h6b = spinnames(int_mb(k_spin+h6b-1))
                      if (nodezero) write(LuOut,9010)
     1                (int_mb(k_mo_index+int_mb(k_offset+p5b-1)+p5-1)+
     1                 2-int_mb(k_spin+p5b-1))/2,
     1                irrep_p5b,spin_p5b,
     2                (int_mb(k_mo_index+int_mb(k_offset+h6b-1)+h6-1)+
     1                 2-int_mb(k_spin+h6b-1))/2,
     2                irrep_h6b,spin_h6b,
     2                dbl_mb(k_b+i-1)
                    endif
      enddo
      enddo
c
      if (.not.ma_pop_stack(l_b))
     1   call errquit('1DM_ph:MA_problem',1,MA_ERR)
c
      END IF
      END IF
      END IF
      END DO
      END DO
c
      write(LuOut,9001)
c
      end if  !nodezero
c
 9000 format(/,1x,a)
 9001 format(/)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
c
      return
      end
c
c
      SUBROUTINE tce_dens_pp_print(d_a,k_a_offset,r1_thresh,irrep)
C     dens_pp ( p6 p5 ) - fortran77 convention
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER p5b,p5
      INTEGER p6b,p6
      integer l_b,k_b,size_b
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_p5b,irrep_p6b
      character*5 spinnames(2),spin_p5b,spin_p6b
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0) !nodezero
      if (nodezero) then
         write(LuOut,9000) '1PDM p-p block largest contributions'
c
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+p6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+p6b-1)) .eq. irrep) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &).ne.4)) THEN
c
      size_b = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1)
       if (.not.ma_push_get(mt_dbl,size_b,'d1hh',l_b,k_b))
     &       call errquit('tce_print_pp_1dm: MA problem',0,MA_ERR)
       call get_hash_block(d_a,dbl_mb(k_b),size_b,
     &       int_mb(k_a_offset),(p6b-noab-1)+nvab*(p5b-noab-1))
      i = 0
      do p5=1, int_mb(k_range+p5b-1)
      do p6=1, int_mb(k_range+p6b-1)
         i=i+1
                    if (dabs(dbl_mb(k_b+i-1)).gt.r1_thresh) then
                      call sym_irrepname(geom,int_mb(k_sym+p5b-1)+1,
     1                                   irrep_p5b)
                      call sym_irrepname(geom,int_mb(k_sym+p6b-1)+1,
     1                                   irrep_p6b)
                      spin_p5b = spinnames(int_mb(k_spin+p5b-1))
                      spin_p6b = spinnames(int_mb(k_spin+p6b-1))
                      if (nodezero) write(LuOut,9010)
     1                (int_mb(k_mo_index+int_mb(k_offset+p5b-1)+p5-1)+
     1                 2-int_mb(k_spin+p5b-1))/2,
     1                irrep_p5b,spin_p5b,
     2                (int_mb(k_mo_index+int_mb(k_offset+p6b-1)+p6-1)+
     1                 2-int_mb(k_spin+p6b-1))/2,
     2                irrep_p6b,spin_p6b,
     2                dbl_mb(k_b+i-1)
                    endif
      enddo
      enddo
c
      if (.not.ma_pop_stack(l_b))
     1   call errquit('1DM_hh:MA_problem',1,MA_ERR)
c
      END IF
      END IF
      END IF
      END DO
      END DO
c 
      write(LuOut,9001)
c
      end if  !nodezero
c
 9000 format(/,1x,a)
 9001 format(/)
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
c
      return
      end

