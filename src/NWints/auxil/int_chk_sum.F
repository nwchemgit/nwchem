      subroutine int_chk_sum(basisin,print_ints)
c $Id: int_chk_sum.F,v 1.3 1994-05-05 00:03:20 d3e129 Exp $
      implicit none
#include "mafdecls.h"
#include "bas.fh"
c::passed
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral and label values 
c::local
      logical debug
      logical MAstatus
      integer max1bf, maxg, mscr1, mscr2
c
      integer i,j,k,l, numleft, numtot, lend
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for labels buffers
      integer i_klab, h_klab, i_llab, h_llab   ! handles and index for labels buffers
c
      integer h_stats, i_stats
c
      integer ncont
      integer plevel
c
      data debug /.false./
      data plevel /0/
c 
c
c if debugging then up print level
c
      if (debug) plevel = 1
c
c get number of contractions (shells) ncont
c
      if (.not.bas_numcont(basisin,ncont))
     &       call errquit('int_chk_sum: numcont not obtained',ncont)
c
      call int_init(1,basisin)
      call int_mem(max1bf,maxg,mscr1,mscr2)
      memb = max(max1bf,maxg)
      mems = max(mscr1,mscr2)
c
c...  pray that ma has been initialized still need verification flag for
c...       consistency.
c
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &       call errquit('int_chk_sum: allocation memb failed',memb)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &       call errquit('int_chk_sum: allocation mems failed',mems)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,'i label buffer',h_ilab,i_ilab))
     &       call errquit('int_chk_sum: allocation ilab failed',memb)
      if (.not.MA_Push_Get(MT_Int,memb,'j label buffer',h_jlab,i_jlab))
     &       call errquit('int_chk_sum: allocation jlab failed',memb)
      if (.not.MA_Push_Get(MT_Int,memb,'k label buffer',h_klab,i_klab))
     &       call errquit('int_chk_sum: allocation klab failed',memb)
      if (.not.MA_Push_Get(MT_Int,memb,'l label buffer',h_llab,i_llab))
     &       call errquit('int_chk_sum: allocation llab failed',memb)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &       h_stats,i_stats))
     &       call errquit('int_chk_sum: allocation mems failed',(5*3))
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
      call ifill(memb,0,    int_mb(i_klab) ,1)
      call ifill(memb,0,    int_mb(i_llab) ,1)
      call dfill(memb,0.0,  dbl_mb(imemb)  ,1)
      call dfill(mems,0.0,  dbl_mb(imems)  ,1)
      call dfill((5*3),0.0, dbl_mb(i_stats),1)
c
c loop over shells and compute S sums
c      
      numtot = 0
      do 00100 i = 1,ncont
        do 00200 j = 1,i
          call int_1eov(basisin,basisin,i,j,
     &           mems,dbl_mb(imems),memb,dbl_mb(imemb))
          call int_lgen1e(basisin,basisin,i,j,1.0d-08,
     &           int_mb(i_ilab),int_mb(i_jlab),memb,
     &           dbl_mb(imemb),mems,dbl_mb(imems),numleft)
          if(print_ints)
     &      call int_pgen1e('overlap',basisin,basisin,i,j,
     &           int_mb(i_ilab),int_mb(i_jlab),numleft,
     &           dbl_mb(imemb),plevel)
          call stsum1(numleft,numtot,
     &           dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &           dbl_mb(i_stats))
00200   continue
00100 continue
      call pstat_sum('overlap',dbl_mb(i_stats),numtot)
c
c loop over shells and compute T sums
c      
      call dfill((5*3),0.0, dbl_mb(i_stats),1)
      numtot = 0
      do 00101 i = 1,ncont
        do 00201 j = 1,i
          call int_1eke(basisin,basisin,i,j,
     &           mems,dbl_mb(imems),memb,dbl_mb(imemb))
          call int_lgen1e(basisin,basisin,i,j,1.0d-08,
     &           int_mb(i_ilab),int_mb(i_jlab),memb,
     &           dbl_mb(imemb),mems,dbl_mb(imems),numleft)
          if(print_ints)
     &      call int_pgen1e('kinetic energy',basisin,basisin,i,j,
     &           int_mb(i_ilab),int_mb(i_jlab),numleft,
     &           dbl_mb(imemb),plevel)
          call stsum1(numleft,numtot,
     &           dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &           dbl_mb(i_stats))
00201   continue
00101 continue
      call pstat_sum('kinetic energy',dbl_mb(i_stats),numtot)
c
c loop over shells and compute V sums
c      
      call dfill((5*3),0.0, dbl_mb(i_stats),1)
      numtot = 0
      do 00102 i = 1,ncont
        do 00202 j = 1,i
          call int_1epe(basisin,basisin,i,j,
     &           mems,dbl_mb(imems),memb,dbl_mb(imemb))
          call int_lgen1e(basisin,basisin,i,j,1.0d-08,
     &           int_mb(i_ilab),int_mb(i_jlab),memb,
     &           dbl_mb(imemb),mems,dbl_mb(imems),numleft)
          if(print_ints)
     &      call int_pgen1e('pot. energy',basisin,basisin,i,j,
     &           int_mb(i_ilab),int_mb(i_jlab),numleft,
     &           dbl_mb(imemb),plevel)
          call stsum1(numleft,numtot,
     &           dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &           dbl_mb(i_stats))
00202   continue
00102 continue
      call pstat_sum('pot. energy',dbl_mb(i_stats),numtot)
c
c loop over shells and compute h1 sums
c      
      call dfill((5*3),0.0, dbl_mb(i_stats),1)
      numtot = 0
      do 00103 i = 1,ncont
        do 00203 j = 1,i
          call int_1eh1(basisin,basisin,i,j,
     &           mems,dbl_mb(imems),memb,dbl_mb(imemb))
          call int_lgen1e(basisin,basisin,i,j,1.0d-08,
     &           int_mb(i_ilab),int_mb(i_jlab),memb,
     &           dbl_mb(imemb),mems,dbl_mb(imems),numleft)
          if(print_ints)
     &      call int_pgen1e('h1',basisin,basisin,i,j,
     &           int_mb(i_ilab),int_mb(i_jlab),numleft,
     &           dbl_mb(imemb),plevel)
          call stsum1(numleft,numtot,
     &           dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &           dbl_mb(i_stats))
00203   continue
00103 continue
      call pstat_sum('h1',dbl_mb(i_stats),numtot)
c
c do canonical shell list for 2e integrals
c
      numtot = 0
      call dfill((5*3),0.0, dbl_mb(i_stats),1)
      do 01000 i = 1,ncont
        do 02000 j = 1, i
          do 03000 k = 1, i
            lend = k
            if (i.eq.k) lend = j
            do 04000 l = 1, lend
              call int_2e4c(basisin,i,j,basisin,k,l,
     &               mems,dbl_mb(imems),memb,dbl_mb(imemb))
              call int_l2e4c(basisin,i,j,basisin,k,l,
     &               1.0d-08,.true.,memb,dbl_mb(imemb),numleft,
     &               int_mb(i_ilab),int_mb(i_jlab),
     &               int_mb(i_klab),int_mb(i_llab),
     &               mems,dbl_mb(imems))
              call stsum2(numleft,numtot,
     &               dbl_mb(imemb),
     &               int_mb(i_ilab),int_mb(i_jlab),
     &               int_mb(i_klab),int_mb(i_llab),
     &               dbl_mb(i_stats),print_ints)
04000       continue
03000     continue
02000   continue
01000 continue
      call pstat_sum('two electron',dbl_mb(i_stats),numtot)      
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_llab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_klab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &       call errquit('int_chk_sum: error poping memory stack',0)
      end
      subroutine pstat_sum(msg,stats,numtot)
      implicit none
c::passed
      character*(*) msg
      integer numtot
      double precision stats(5,3)
c::local
      write(6,10000)
      write(6,*)' summary statistics for ',msg,' integrals'
      write(6,*)' total number of integrals processed: ',numtot
      write(6,10001)
      write(6,10002)stats
      write(6,10001)
10000 format(/,75('-'))
10001 format(75('-'))
10002 format(5(1pd12.5,3x))
      end      
      subroutine stsum1(num,numtot,Gen,ilab,jlab,stats)
      implicit none
c::passed
      integer num                   ! number of 1e moieties
c
c...   numtot must be zeroed prior to first call to stsum1 ...
c
      integer numtot                ! accumulater
c
      integer ilab(num),jlab(num)   ! labels for moieties
      double precision Gen(num)     ! 1e moieties
      double precision stats(5,3)   ! stats array (zeroed elsewhere)
c::local 
      integer ii,i,j
      double precision x,x2,denom,s1,d1
c
      numtot = numtot + num
      do 00100 ii = 1,num
        i     = ilab(ii)
        j     = jlab(ii)
        s1    = i+j
        d1    = abs(i-j)
        x     = Gen(ii)
        x2    = x*x
        denom = 1.0d00 + x2
c
        stats(1,1) = stats(1,1) + x
        stats(2,1) = stats(2,1) + abs(x)
        stats(3,1) = stats(3,1) + x2
        stats(4,1) = stats(4,1) + x/denom
        stats(5,1) = stats(5,1) + x2/denom
c
        stats(1,2) = stats(1,2) + x        * s1
        stats(2,2) = stats(2,2) + abs(x)   * s1
        stats(3,2) = stats(3,2) + x2       * s1
        stats(4,2) = stats(4,2) + x/denom  * s1
        stats(5,2) = stats(5,2) + x2/denom * s1
c
        stats(1,3) = stats(1,3) + x        * d1
        stats(2,3) = stats(2,3) + abs(x)   * d1
        stats(3,3) = stats(3,3) + x2       * d1
        stats(4,3) = stats(4,3) + x/denom  * d1
        stats(5,3) = stats(5,3) + x2/denom * d1
c
00100 continue
      end
      subroutine stsum2(num,numtot,Eri,ilab,jlab,klab,llab,
     &       stats,printem)
      implicit none
c::passed
      integer num                   ! number of 2e moieties
c
c...   numtot must be zeroed prior to first call to stsum1 ...
c
      integer numtot                ! accumulater
c
      integer ilab(num),jlab(num)   ! labels for moieties
      integer klab(num),llab(num)   ! labels for moieties
      double precision Eri(num)     ! 2e moieties
      double precision stats(5,3)   ! stats array (zeroed elsewhere)
      logical printem               ! print values as you go??
c::local 
      integer ii,i,j,k,l,l4
      double precision x,x2,denom,s1,d1
c
      integer isym2,isym4
      ISYM2(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
      ISYM4(I,J,K,L)=MAX(ISYM2(I,J),ISYM2(K,L))*
     &               (MAX(ISYM2(I,J),ISYM2(K,L))-1)/2+
     &               MIN(ISYM2(I,J),ISYM2(K,L))
c
      numtot = numtot + num
      do 00100 ii = 1,num
        i     = ilab(ii)
        j     = jlab(ii)
        k     = klab(ii)
        l     = llab(ii)
        s1    = i+j+k+l
        d1    = abs(abs(i-j)-abs(k-l))
        x     = Eri(ii)
        if (printem) then
          l4 = isym4(i,j,k,l)
          write(6,10001)l4,i,j,k,l,x
10001    format(i8,4(i4,1x),5x,'eri = ',1pd20.10)
       endif
        x2    = x*x
        denom = 1.0d00 + x2
c
        stats(1,1) = stats(1,1) + x
        stats(2,1) = stats(2,1) + abs(x)
        stats(3,1) = stats(3,1) + x2
        stats(4,1) = stats(4,1) + x/denom
        stats(5,1) = stats(5,1) + x2/denom
c
        stats(1,2) = stats(1,2) + x        * s1
        stats(2,2) = stats(2,2) + abs(x)   * s1
        stats(3,2) = stats(3,2) + x2       * s1
        stats(4,2) = stats(4,2) + x/denom  * s1
        stats(5,2) = stats(5,2) + x2/denom * s1
c
        stats(1,3) = stats(1,3) + x        * d1
        stats(2,3) = stats(2,3) + abs(x)   * d1
        stats(3,3) = stats(3,3) + x2       * d1
        stats(4,3) = stats(4,3) + x/denom  * d1
        stats(5,3) = stats(5,3) + x2/denom * d1
c
00100 continue
      end

