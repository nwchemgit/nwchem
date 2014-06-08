ckbn  Kiran Bhaskaran-Nair modified from tce_eom_xguess.F
      subroutine tce_eom_ipxguess(rtdb,needx1,needx2,needx3,needx4,
     1  size_x1,size_x2,size_x3,size_x4,
     2  k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
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
#include "rtdb.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "errquit.fh"
#include "util.fh"
      logical guessx2
      parameter (guessx2 = .true.)
      logical needx1
      logical needx2
      logical needx3
      logical needx4
      integer rtdb
      integer size_x1
      integer size_x2
      integer size_x3
      integer size_x4
      integer k_x1_offset
      integer k_x2_offset
      integer k_x3_offset
      integer k_x4_offset
      logical nodezero
      integer p1b,p2b
      integer h3b,h4b
      integer p1,p2,h3,h4
      integer l_x1,k_x1
      integer l_x2,k_x2
      integer l_diff,k_diff
      integer size
      integer ivec, jvec
c      double precision maxdiff
      double precision nextmaxdiff,maxeorb
      integer i,ibumper
      logical defmeo
      character*255 filename
c
c     Determine threshold
c
      nodezero = (ga_nodeid().eq.0)
      ibumper=200000
cx      if (.not.ma_push_get(mt_dbl,maxtrials,'diff',l_diff,k_diff))
cx     1  call errquit('tce_eom_xguess: MA problem',0,MA_ERR)
cx      call dfill(maxtrials,1.0d99,dbl_mb(k_diff),1)
      if (.not.ma_push_get(mt_dbl,ibumper,'diff',l_diff,k_diff))
     1  call errquit('tce_eom_xguess: MA problem',0,MA_ERR)
      call dfill(ibumper,1.0d99,dbl_mb(k_diff),1)
      ivec = 0
ckbn      do p1b = noab+1,noab+nvab
      do p1b = 1,1
        do h3b = 1,noab
ckbn          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) then
          if (ip_unused_spin  .eq. int_mb(k_spin+h3b-1)) then
ckbn          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
ckbn     1      +int_mb(k_spin+h3b-1).ne.4)) then
          if ((.not.restricted).or.(ip_unused_spin 
     1      +int_mb(k_spin+h3b-1).ne.4)) then
ckbn          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)).eq.
ckbn     1      irrep_x) then
          if (ieor(ip_unused_sym ,int_mb(k_sym+h3b-1)).eq.
     1      irrep_x) then
ckbn          do p1 = 1,int_mb(k_range+p1b-1)
            do h3 = 1,int_mb(k_range+h3b-1)
              ivec = ivec + 1
              dbl_mb(k_diff+ivec-1) = 
ckbn     1        (+dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
ckbn     2         -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1))
     2         -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
            enddo
ckbn          enddo
          endif
          endif
          endif
        enddo
      enddo
      if (model.eq.'cis') then
        nroots_reduced = ivec
      else
        nroots_reduced = min(ivec, nroots)
      endif
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
     1  call errquit('tce_eom_xguess: MA problem',1,MA_ERR)
c
       maxeorb=0.0d0
       defmeo=.false.
       if (.not.rtdb_get(rtdb,'tce:maxeorb',mt_dbl,1,maxeorb)) then
         defmeo=.false.
       else
         defmeo=.true.
         maxdiff=maxeorb 
       end if
c
c
c     Make trial X1
c
      nxtrials = 0
ckbn      do p1b = noab+1,noab+nvab
      do p1b = 1,1
        do h3b = 1,noab
ckbn          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) then
          if (ip_unused_spin .eq. int_mb(k_spin+h3b-1)) then
ckbn          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
ckbn     1      +int_mb(k_spin+h3b-1).ne.4)) then
          if ((.not.restricted).or.(ip_unused_spin 
     1      +int_mb(k_spin+h3b-1).ne.4)) then
ckbn          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)).eq.
ckbn     1      irrep_x) then
          if (ieor(ip_unused_sym ,int_mb(k_sym+h3b-1)).eq.
     1      irrep_x) then
ckbn          size = int_mb(k_range+p1b-1) * int_mb(k_range+h3b-1)
          size = 1 * int_mb(k_range+h3b-1)
          if (.not.ma_push_get(mt_dbl,size,'x1',l_x1,k_x1))
     1      call errquit('tce_eom_xguess: MA problem',2,MA_ERR)
          ivec = 0
ckbn          do p1 = 1,int_mb(k_range+p1b-1)
          do p1 = 1,1
            do h3 = 1,int_mb(k_range+h3b-1)
              ivec = ivec + 1
ckbn              if (dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
              if ( !dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     1           -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     2           .le.maxdiff) then
              nxtrials = nxtrials + 1
              call dfill(size,0.0d0,dbl_mb(k_x1),1)
              dbl_mb(k_x1+ivec-1) = 1.0d0
              if (nxtrials.gt.maxtrials) call errquit
     1          ('tce_eom_xguess: problem too large',nxtrials,
     2          CALC_ERR)
              if (needx1) then
                call tce_filenameindexed(nxtrials,'x1',filename)
                call createfile(filename,x1(nxtrials),size_x1)
                call put_hash_block(x1(nxtrials),dbl_mb(k_x1),size,
     1            int_mb(k_x1_offset),((p1b-noab-1)*noab+h3b-1))
                call gatodra(x1(nxtrials))
                call gatoeaf(x1(nxtrials))
                x1_exist(nxtrials) = .true.
              endif
              if (needx2) then
                call tce_filenameindexed(nxtrials,'x2',filename)
                call createfile(filename,x2(nxtrials),size_x2)
                call gatodra(x2(nxtrials))
                call gatoeaf(x2(nxtrials))
                x2_exist(nxtrials) = .true.
              endif
              if (needx3) then
                call tce_filenameindexed(nxtrials,'x3',filename)
                call createfile(filename,x3(nxtrials),size_x3)
                call gatodra(x3(nxtrials))
                call gatoeaf(x3(nxtrials))
                x3_exist(nxtrials) = .true.
              endif
              if (needx4) then
                call tce_filenameindexed(nxtrials,'x4',filename)
                call createfile(filename,x4(nxtrials),size_x4)
                call gatodra(x4(nxtrials))
                call gatoeaf(x4(nxtrials))
                x4_exist(nxtrials) = .true.
              endif
              endif
            enddo
          enddo
          if (.not.ma_pop_stack(l_x1))
     1      call errquit('tce_eom_xguess: MA problem',3,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
c
      if(defmeo) go to 300
c
c     Make trial X2
c
      maxdiff=0.30d0
      if (guessx2.and.needx2) then
      do p1b = noab+1,noab+nvab
ckbn        do p2b = p1b,noab+nvab
        do p2b = 1,1
          do h3b = 1,noab
            do h4b = h3b,noab
ckbn              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
ckbn     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
              if (int_mb(k_spin+p1b-1)+ ip_unused_spin 
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
ckbn              if ((.not.restricted).or.
ckbn     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
ckbn     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+ ip_unused_spin +
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
ckbn              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
ckbn     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
ckbn     2          .eq. irrep_x) then
              if (ieor(int_mb(k_sym+p1b-1),ieor( ip_unused_sym ,
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep_x) then
ckbn              size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
ckbn     1             * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
              size = int_mb(k_range+p1b-1) * 1
     1             * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
              if (.not.ma_push_get(mt_dbl,size,'x2',l_x2,k_x2))
     1          call errquit('tce_eom_xguess: MA problem',0,MA_ERR)
              do p1 = 1,int_mb(k_range+p1b-1)
ckbn                do p2 = 1,int_mb(k_range+p2b-1)
                do p2 = 1,1
                  do h3 = 1,int_mb(k_range+h3b-1)
                    do h4 = 1,int_mb(k_range+h4b-1)
                if (+dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
ckbn     1              +dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     3              -dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     4              .le. maxdiff*2.0d0) then
ckbn                if (((p1b.ne.p2b).or.(p1.lt.p2)).and.
ckbn     1              ((h3b.ne.h4b).or.(h3.lt.h4))) then
                if ( 
     1              ((h3b.ne.h4b).or.(h3.lt.h4))) then
                  nxtrials = nxtrials + 1
                  call dfill(size,0.0d0,dbl_mb(k_x2),1)
ckbn                  if ((p1b.eq.p2b).and.(h3b.eq.h4b)) then 
ckbn                  dbl_mb(k_x2+(((p1-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p2-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
ckbn     3              = 0.5d0
ckbn                  dbl_mb(k_x2+(((p2-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p1-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
ckbn     3              = -0.5d0
ckbn                  dbl_mb(k_x2+(((p1-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p2-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h4-1))*int_mb(k_range+h4b-1)+h3-1)
ckbn     3              = -0.5d0
ckbn                  dbl_mb(k_x2+(((p2-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p1-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h4-1))*int_mb(k_range+h4b-1)+h3-1)
ckbn     3              = 0.5d0
ckbn                  else if (p1b.eq.p2b) then
ckbn                  dbl_mb(k_x2+(((p1-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p2-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
ckbn     3              = dsqrt(0.5d0)
ckbn                  dbl_mb(k_x2+(((p2-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p1-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
ckbn     3              = -dsqrt(0.5d0)
ckbn                  else if (h3b.eq.h4b) then
ckbn                  dbl_mb(k_x2+(((p1-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p2-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
ckbn     3              = dsqrt(0.5d0)
ckbn                  dbl_mb(k_x2+(((p1-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p2-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h4-1))*int_mb(k_range+h4b-1)+h3-1)
ckbn     3              = -dsqrt(0.5d0)
ckbn                  else
ckbn                  dbl_mb(k_x2+(((p1-1)*int_mb(k_range+p2b-1)
ckbn     1                         +(p2-1))*int_mb(k_range+h3b-1)
ckbn     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
ckbn     3              = 1.0d0
ckbn                  endif
                  if ((h3b.eq.h4b)) then 
                  dbl_mb(k_x2+(((p1-1)*1
     1                         +(p2-1))*int_mb(k_range+h3b-1)
     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
     3              = 0.5d0
                  dbl_mb(k_x2+(((p2-1)*1
     1                         +(p1-1))*int_mb(k_range+h3b-1)
     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
     3              = -0.5d0
                  dbl_mb(k_x2+(((p1-1)*1
     1                         +(p2-1))*int_mb(k_range+h3b-1)
     2                         +(h4-1))*int_mb(k_range+h4b-1)+h3-1)
     3              = -0.5d0
                  dbl_mb(k_x2+(((p2-1)*1
     1                         +(p1-1))*int_mb(k_range+h3b-1)
     2                         +(h4-1))*int_mb(k_range+h4b-1)+h3-1)
     3              = 0.5d0
                  else if (p1b.eq.p2b) then
                  dbl_mb(k_x2+(((p1-1)*1
     1                         +(p2-1))*int_mb(k_range+h3b-1)
     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
     3              = dsqrt(0.5d0)
                  dbl_mb(k_x2+(((p2-1)*1
     1                         +(p1-1))*int_mb(k_range+h3b-1)
     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
     3              = -dsqrt(0.5d0)
                  else if (h3b.eq.h4b) then
                  dbl_mb(k_x2+(((p1-1)*1
     1                         +(p2-1))*int_mb(k_range+h3b-1)
     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
     3              = dsqrt(0.5d0)
                  dbl_mb(k_x2+(((p1-1)*1
     1                         +(p2-1))*int_mb(k_range+h3b-1)
     2                         +(h4-1))*int_mb(k_range+h4b-1)+h3-1)
     3              = -dsqrt(0.5d0)
                  else
                  dbl_mb(k_x2+(((p1-1)*1
     1                         +(p2-1))*int_mb(k_range+h3b-1)
     2                         +(h3-1))*int_mb(k_range+h4b-1)+h4-1)
     3              = 1.0d0
                  endif
                  if (nxtrials.gt.maxtrials) call errquit
     1              ('tce_eom_xguess: problem too large',nxtrials,
     2              CALC_ERR)
                  if (needx1) then
                    call tce_filenameindexed(nxtrials,'x1',filename)
                    call createfile(filename,x1(nxtrials),size_x1)
                    call gatodra(x1(nxtrials))
                    call gatoeaf(x1(nxtrials))
                    x1_exist(nxtrials) = .true.
                  endif
                  if (needx2) then
                    call tce_filenameindexed(nxtrials,'x2',filename)
                    call createfile(filename,x2(nxtrials),size_x2)
                    call put_hash_block(x2(nxtrials),dbl_mb(k_x2),size,
     1                int_mb(k_x2_offset),((((p1b-noab-1)*nvab
     2                +p2b-noab-1)*noab+h3b-1)*noab+h4b-1))
                    call gatodra(x2(nxtrials))
                    call gatoeaf(x2(nxtrials))
                    x2_exist(nxtrials) = .true.
                  endif
                  if (needx3) then
                    call tce_filenameindexed(nxtrials,'x3',filename)
                    call createfile(filename,x3(nxtrials),size_x3)
                    call gatodra(x3(nxtrials))
                    call gatoeaf(x3(nxtrials))
                    x3_exist(nxtrials) = .true.
                  endif
                  if (needx4) then
                    call tce_filenameindexed(nxtrials,'x4',filename)
                    call createfile(filename,x4(nxtrials),size_x4)
                    call gatodra(x4(nxtrials))
                    call gatoeaf(x4(nxtrials))
                    x4_exist(nxtrials) = .true.
                  endif
                endif
                endif
                    enddo
                  enddo
                enddo
              enddo
              if (.not.ma_pop_stack(l_x2))
     1          call errquit('tce_eom_xguess: MA problem',3,MA_ERR)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
c
 300  continue
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
