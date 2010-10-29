      subroutine cis_xguess(size_x1,k_x1_offset,iopen_f)
c
c $Id$
c
c     Generates initial trial vectors that are all singly excited 
c     determinants within the same irrep_x
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
c      logical guessx2
c      parameter (guessx2 = .true.)
c      logical needx1
c      logical needx2
c      logical needx3
c      logical needx4
      integer size_x1
c      integer size_x2
c      integer size_x3
c      integer size_x4
      integer k_x1_offset
c      integer k_x2_offset
c      integer k_x3_offset
c      integer k_x4_offset
      logical nodezero
      integer p1b,p2b
      integer h3b,h4b
      integer p1,p2,h3,h4
      integer l_x1,k_x1
      integer l_x2,k_x2
      integer l_diff,k_diff
      integer size
      integer ivec, jvec,iopen_f
c     double precision maxdiff
      double precision nextmaxdiff
      integer i,ibumper
      character*255 filename
c
c     Determine threshold
c
      nodezero = (ga_nodeid().eq.0)
      ibumper=200000
      if (.not.ma_push_get(mt_dbl,ibumper,'diff',l_diff,k_diff))
     1  call errquit('tce_eom_xguess: MA problem',0,MA_ERR)
      call dfill(maxtrials,1.0d99,dbl_mb(k_diff),1)
      ivec = 0
      do p1b = noab+1,noab+nvab
        do h3b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h3b-1).ne.4)) then
              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)).eq.
     1          irrep_x) then
                do p1 = 1,int_mb(k_range+p1b-1)
                  do h3 = 1,int_mb(k_range+h3b-1)
                    ivec = ivec + 1
                   dbl_mb(k_diff+ivec-1) = 
     1             (+dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1))
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
     1  call errquit('cis_xguess: MA problem',1,MA_ERR)
c
c     Make trial X1
c
ccx      maxdiff = 0.800d0 !uracil
c
      nxtrials = 0
      do p1b = noab+1,noab+nvab
        do h3b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h3b-1).ne.4)) then
              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)).eq.
     1          irrep_x) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+h3b-1)
                if (.not.ma_push_get(mt_dbl,size,'x1',l_x1,k_x1))
     1          call errquit('cis_xguess: MA problem',2,MA_ERR)
                ivec = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do h3 = 1,int_mb(k_range+h3b-1)
                    ivec = ivec + 1
                    if (dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     1                 -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     2                 .le.maxdiff) then
                      nxtrials = nxtrials + 1
                      call dfill(size,0.0d0,dbl_mb(k_x1),1)
                      dbl_mb(k_x1+ivec-1) = 1.0d0
                      if (nxtrials.gt.maxtrials) call errquit
     1                  ('cis_xguess: problem too large',nxtrials,
     2                  CALC_ERR)
        if(iopen_f.eq.1) then
                      call tce_filenameindexed(nxtrials,'x1',filename)
                      call createfile(filename,x1(nxtrials),size_x1)
                      call put_hash_block(x1(nxtrials),dbl_mb(k_x1),
     1                     size,int_mb(k_x1_offset),
     2                     ((p1b-noab-1)*noab+h3b-1))
                      call gatodra(x1(nxtrials))
                      call gatoeaf(x1(nxtrials))
                      x1_exist(nxtrials) = .true.
         end if
                    endif
                  enddo
                enddo
                if (.not.ma_pop_stack(l_x1))
     1          call errquit('cis_xguess: MA problem',3,MA_ERR)
              endif
            endif
          endif
        enddo
      enddo
c
      if (nxtrials.lt.nroots_reduced)
     1  call errquit('there is a bug in the program',0,UNKNOWN_ERR)
      if (nodezero) then
        write(LuOut,9000) nxtrials
        if (nroots.gt.nroots_reduced) then
          write(LuOut,9010) nroots, nroots_reduced
        endif
      endif
      return
 9000 format(/,1x,'No. of initial right vectors ',i4)
 9010 format(  1x,'No. of roots reduced from    ',i4,' to ',i4)
      end
