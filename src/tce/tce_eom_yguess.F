      subroutine tce_eom_yguess(needy1,needy2,needy3,needy4,
     1  size_y1,size_y2,size_y3,size_y4,
     2  k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
c
c $Id$
c
c     Generates initial trial vectors that are all singly excited 
c     determinants within the same irrep_y
c
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "errquit.fh"
#include "util.fh"
      logical needy1
      logical needy2
      logical needy3
      logical needy4
      integer size_y1
      integer size_y2
      integer size_y3
      integer size_y4
      integer k_y1_offset
      integer k_y2_offset
      integer k_y3_offset
      integer k_y4_offset
      logical nodezero
      integer p1b
      integer h3b
      integer p1,h3
      integer l_y1,k_y1
      integer l_diff,k_diff
      integer size
      integer ivec, jvec
c      double precision maxdiff
      double precision nextmaxdiff
      integer i,ibumper
      character*255 filename
c
c     Determine threshold
c
      nodezero = (ga_nodeid().eq.0)
      ibumper=200000
ccx      if (.not.ma_push_get(mt_dbl,maxtrials,'diff',l_diff,k_diff))
ccx     1  call errquit('tce_eom_yguess: MA problem',0,MA_ERR)
ccx      call dfill(maxtrials,1.0d99,dbl_mb(k_diff),1)
      if (.not.ma_push_get(mt_dbl,ibumper,'diff',l_diff,k_diff))
     1  call errquit('tce_eom_xguess: MA problem',0,MA_ERR)
      call dfill(ibumper,1.0d99,dbl_mb(k_diff),1)
      ivec = 0
      do p1b = noab+1,noab+nvab
        do h3b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h3b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)).eq.
     1      irrep_y) then
          do p1 = 1,int_mb(k_range+p1b-1)
            do h3 = 1,int_mb(k_range+h3b-1)
              ivec = ivec + 1
              dbl_mb(k_diff+ivec-1) = 
     1        (+dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2         -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1))
            enddo
          enddo
          endif
          endif
          endif
        enddo
      enddo
      nroots_reduced = min(ivec, nroots)
      maxdiff = 0.0d0
  100 continue
      nextmaxdiff = 1.0d99
      jvec = 0
      do i = 1,ivec
        if (dbl_mb(k_diff+i-1).lt.maxdiff) jvec = jvec + 1
        if (dbl_mb(k_diff+i-1).ge.maxdiff) then
          if (dbl_mb(k_diff+i-1).lt.nextmaxdiff)
     1      nextmaxdiff = dbl_mb(k_diff+i-1)
        endif
      enddo
      if (jvec.ge.nroots_reduced) goto 110
      maxdiff = nextmaxdiff + 0.001d0
      goto 100
  110 continue
      if (.not.ma_pop_stack(l_diff))
     1  call errquit('tce_eom_yguess: MA problem',1,MA_ERR)
c
c     Make trial Y1
c
      nytrials = 0
      do h3b = 1,noab
        do p1b = noab+1,noab+nvab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h3b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)).eq.
     1      irrep_y) then
          size = int_mb(k_range+p1b-1) * int_mb(k_range+h3b-1)
          if (.not.ma_push_get(mt_dbl,size,'y1',l_y1,k_y1))
     1      call errquit('tce_eom_yguess: MA problem',2,MA_ERR)
          ivec = 0
          do h3 = 1,int_mb(k_range+h3b-1)
            do p1 = 1,int_mb(k_range+p1b-1)
              ivec = ivec + 1
              if (dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     1           -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     2           .le.maxdiff) then
              nytrials = nytrials + 1
              call dfill(size,0.0d0,dbl_mb(k_y1),1)
              dbl_mb(k_y1+ivec-1) = 1.0d0
              if (nytrials.gt.maxtrials) call errquit
     1          ('tce_eom_yguess: problem too large',nytrials,
     2          CALC_ERR)
              if (needy1) then
                call tce_filenameindexed(nytrials,'y1',filename)
                call createfile(filename,y1(nytrials),size_y1)
                call put_hash_block(y1(nytrials),dbl_mb(k_y1),size,
     1            int_mb(k_y1_offset),((h3b-1)*nvab+p1b-noab-1))
                call gatodra(y1(nytrials))
                call gatoeaf(y1(nytrials))
                y1_exist(nytrials) = .true.
              endif
              if (needy2) then
                call tce_filenameindexed(nytrials,'y2',filename)
                call createfile(filename,y2(nytrials),size_y2)
                call gatodra(y2(nytrials))
                call gatoeaf(y2(nytrials))
                y2_exist(nytrials) = .true.
              endif
              if (needy3) then
                call tce_filenameindexed(nytrials,'y3',filename)
                call createfile(filename,y3(nytrials),size_y3)
                call gatodra(y3(nytrials))
                call gatoeaf(y3(nytrials))
                y3_exist(nytrials) = .true.
              endif
              if (needy4) then
                call tce_filenameindexed(nytrials,'y4',filename)
                call createfile(filename,y4(nytrials),size_y4)
                call gatodra(y4(nytrials))
                call gatoeaf(y4(nytrials))
                y4_exist(nytrials) = .true.
              endif
              endif
            enddo
          enddo
          if (.not.ma_pop_stack(l_y1))
     1      call errquit('tce_eom_yguess: MA problem',3,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      if (nytrials.lt.nroots_reduced)
     1  call errquit('there is a bug in the program',0,UNKNOWN_ERR)
      if (nodezero) then
        write(LuOut,9000) nytrials
        if (nroots.gt.nroots_reduced) then
          write(LuOut,9010) nroots, nroots_reduced
        endif
      endif
      return
 9000 format(/,1x,'No. of initial left vectors ',i4)
 9010 format(  1x,'No. of roots reduced from    ',i4,' to ',i4)
      end
