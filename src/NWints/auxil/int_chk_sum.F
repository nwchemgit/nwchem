      subroutine int_chk_sum(rtdb,basisin,print_ints)
c $Id: int_chk_sum.F,v 1.19 1996-01-05 02:00:12 d3e129 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral and label values 
c::local
      logical debug
      logical MAstatus
      integer max1bf, maxg, mscr1, mscr2, mscrh1, maxh1
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
      logical print69
      logical dumm69
c
      integer ncont
      integer plevel
c
      integer me, nproc, idoit
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
     &    call errquit('int_chk_sum: numcont not obtained',ncont)
c
      call int_mem(max1bf,maxg,mscr1,mscr2)
      memb = max(max1bf,maxg)
      mems = max(mscr1,mscr2)
      call int_mem_h1(maxh1, mscrh1)
      memb = max(memb,maxh1)
      mems = max(mems,mscrh1)
      call int_mem_2e3c(maxg,mscr2)
      memb = max(memb,maxg)
      mems = max(mems,mscr2)
      call int_mem_2e2c(maxg,mscr2)
      memb = max(memb,maxg)
      mems = max(mems,mscr2)
c
c...  pray that ma has been initialized still need verification flag for
c...       consistency.
c
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit('int_chk_sum: allocation memb failed',memb)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit('int_chk_sum: allocation mems failed',mems)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,'i label buffer',h_ilab,i_ilab))
     &    call errquit('int_chk_sum: allocation ilab failed',memb)
      if (.not.MA_Push_Get(MT_Int,memb,'j label buffer',h_jlab,i_jlab))
     &    call errquit('int_chk_sum: allocation jlab failed',memb)
      if (.not.MA_Push_Get(MT_Int,memb,'k label buffer',h_klab,i_klab))
     &    call errquit('int_chk_sum: allocation klab failed',memb)
      if (.not.MA_Push_Get(MT_Int,memb,'l label buffer',h_llab,i_llab))
     &    call errquit('int_chk_sum: allocation llab failed',memb)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit('int_chk_sum: allocation mems failed',(5*3))
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
      call ifill(memb,0,    int_mb(i_klab) ,1)
      call ifill(memb,0,    int_mb(i_llab) ,1)
      call dfill(memb,0.0d0,  dbl_mb(imemb)  ,1)
      call dfill(mems,0.0d0,  dbl_mb(imems)  ,1)
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
c
c initialize parallelization flags
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      call ga_sync()
      if ((me.eq.0).and.(nproc.gt.1).and.print_ints)
     &    write(6,*)' only integrals on node 0 will be printed. '
c
c loop over shells and compute S sums
c      
      print69 = .false.
      if (rtdb_get(rtdb,'int:ov69',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      idoit = (me-1)
      numtot = 0
      do 00100 i = 1,ncont
        do 00200 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
*            call int_1eov(basisin,i,basisin,j,
*     &          mems,dbl_mb(imems),memb,dbl_mb(imemb))
*            call int_lgen1e(basisin,i,basisin,j,1.0d-08,
*     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call int_l1eov(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            if(print_ints.and.me.eq.0)
     &          call int_pgen1e('overlap',basisin,i,basisin,j,
     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
     &          dbl_mb(imemb),plevel)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),(print69.and.(me.eq.0)),
     &          'overlap')
          endif
00200   continue
00100 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('overlap',dbl_mb(i_stats),numtot)
      call ga_sync()
c
c
c loop over shells and compute T sums
c      
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
      print69 = .false.
      if (rtdb_get(rtdb,'int:ke69',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      idoit = (me-1)
      numtot = 0
      do 00101 i = 1,ncont
        do 00201 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
*            call int_1eke(basisin,i,basisin,j,
*     &          mems,dbl_mb(imems),memb,dbl_mb(imemb))
*            call int_lgen1e(basisin,i,basisin,j,1.0d-08,
*     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call int_l1eke(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            if(print_ints.and.me.eq.0)
     &          call int_pgen1e('kinetic energy',basisin,i,basisin,j,
     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
     &          dbl_mb(imemb),plevel)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),(print69.and.(me.eq.0)),
     &          'kinetic')
          endif
00201   continue
00101 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('kinetic energy',dbl_mb(i_stats),numtot)
      call ga_sync()
c
c loop over shells and compute V sums
c      
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
      print69 = .false.
      if (rtdb_get(rtdb,'int:pe69',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      idoit = (me-1)
      numtot = 0
      do 00102 i = 1,ncont
        do 00202 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
*            call int_1epe(basisin,i,basisin,j,
*     &          mems,dbl_mb(imems),memb,dbl_mb(imemb))
*            call int_lgen1e(basisin,i,basisin,j,1.0d-08,
*     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call int_l1epe(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            if(print_ints.and.me.eq.0)
     &          call int_pgen1e('pot. energy',basisin,i,basisin,j,
     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
     &          dbl_mb(imemb),plevel)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),(print69.and.(me.eq.0))
     &          ,'potential')
          endif
00202   continue
00102 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('pot. energy',dbl_mb(i_stats),numtot)
      call ga_sync()
c
c loop over shells and compute h1 sums
c      
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
      idoit = (me-1)
      print69 = .false.
      if (rtdb_get(rtdb,'int:h169',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      numtot = 0
      do 00103 i = 1,ncont
        do 00203 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
*            call int_1eh1(basisin,i,basisin,j,
*     &          mems,dbl_mb(imems),memb,dbl_mb(imemb))
*            call int_lgen1e(basisin,i,basisin,j,1.0d-08,
*     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call int_l1eh1(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            if(print_ints.and.me.eq.0)
     &          call int_pgen1e('h1',basisin,i,basisin,j,
     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
     &          dbl_mb(imemb),plevel)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),(print69.and.(me.eq.0)),
     &          ' H1 (T+V)')
          endif
00203   continue
00103 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('h1',dbl_mb(i_stats),numtot)
      call ga_sync()
c
c do canonical shell list for 2e integrals (4 center)
c
      print69 = .false.
      if (rtdb_get(rtdb,'int:4c69',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      idoit = (me-1)
      numtot = 0
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
      do 01000 i = 1,ncont
        if((me.eq.0).and.(mod(i,10).eq.0)) then
          write(6,*)' 2e4c outer loop counter ',i,' of ',ncont
          call util_flush(6)
        endif
        do 02000 j = 1, i
          do 03000 k = 1, i
            lend = k
            if (i.eq.k) lend = j
            do 04000 l = 1, lend
              idoit = idoit + 1
              if(mod(idoit,nproc).eq.0) then
*                call int_2e4c(basisin,i,j,basisin,k,l,
*     &              mems,dbl_mb(imems),memb,dbl_mb(imemb))
                call int_l2e4c(basisin,i,j,basisin,k,l,
c     &              1.0d-08,.true.,memb,dbl_mb(imemb),numleft,
     &              1.0d-08,.false.,memb,dbl_mb(imemb),numleft,
     &              int_mb(i_ilab),int_mb(i_jlab),
     &              int_mb(i_klab),int_mb(i_llab),
     &              mems,dbl_mb(imems))
                call stsum2(numleft,numtot,
     &              dbl_mb(imemb),
     &              int_mb(i_ilab),int_mb(i_jlab),
     &              int_mb(i_klab),int_mb(i_llab),
     &              dbl_mb(i_stats),
     &              (print_ints.and.me.eq.0),(print69.and.me.eq.0))
              endif
04000       continue
03000     continue
02000   continue
01000 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('two electron',dbl_mb(i_stats),numtot)      
      call ga_sync()
c
c do canonical shell list for 2e integrals (3 center)
c
      print69 = .false.
      if (rtdb_get(rtdb,'int:3c69',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      idoit = (me-1)
      numtot = 0
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
      do 01001 i = 1,ncont
        if((me.eq.0).and.(mod(i,10).eq.0)) then
          write(6,*)' 2e3c outer loop counter ',i,' of ',ncont
          call util_flush(6)
        endif
        do 02001 j = 1, i
          do 03001 k = 1, j
            idoit = idoit + 1
            if(mod(idoit,nproc).eq.0) then
*              call int_2e3c(basisin,i,basisin,j,k,
*     &              mems,dbl_mb(imems),memb,dbl_mb(imemb))
              call int_l2e3c(basisin,i,basisin,j,k,
     &            1.0d-08,.true.,memb,dbl_mb(imemb),numleft,
     &            int_mb(i_ilab),int_mb(i_jlab),
     &            int_mb(i_klab),
     &            mems,dbl_mb(imems))
              call stsum2_3c(numleft,numtot,
     &            dbl_mb(imemb),
     &            int_mb(i_ilab),int_mb(i_jlab),
     &            int_mb(i_klab),
     &            dbl_mb(i_stats),(print_ints.and.me.eq.0),
     &            (print69.and.(me.eq.0)))
            endif
03001     continue
02001   continue
01001 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum
     &    ('3 center two electron',dbl_mb(i_stats),numtot)      
      call ga_sync()
c
c do canonical shell list for 2e integrals (2 center)
c
      print69 = .false.
      if (rtdb_get(rtdb,'int:2c69',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
      idoit = (me-1)
      numtot = 0
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
      do 01002 i = 1,ncont
        if((me.eq.0).and.(mod(i,10).eq.0)) then
          write(6,*)' 2e2c outer loop counter ',i,' of ',ncont
          call util_flush(6)
        endif
        do 02002 j = 1, i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
*            call int_2e2c(basisin,i,basisin,j,
*     &          mems,dbl_mb(imems),memb,dbl_mb(imemb))
            call int_l2e2c(basisin,i,basisin,j,
     &          1.0d-08,memb,dbl_mb(imemb),numleft,
     &          int_mb(i_ilab),int_mb(i_jlab),
     &          mems,dbl_mb(imems))
            call stsum2_2c(numleft,numtot,
     &          dbl_mb(imemb),
     &          int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),(print_ints.and.me.eq.0),
     &            (print69.and.(me.eq.0)))
          endif
02002   continue
01002 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum
     &    ('2 center two electron',dbl_mb(i_stats),numtot)      
      call ga_sync()
      call ga_sync()
      call ga_sync()
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
      call ga_sync()
      call ga_sync()
      end
      subroutine prt_stat_sum(msg,stats,numtot)
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
      call util_flush(6)
10000 format(/,90('-'))
10001 format(90('-'))
10002 format(5(1pd15.8,3x))
      end      
      subroutine stsum1(num,numtot,Gen,ilab,jlab,stats,print69,tag69)
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
      logical print69
      character*(*) tag69
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
        if(print69)
     &      write(69,'(1x,i5,i5,d20.10,3x,a)') i, j, x, tag69
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
      subroutine stsum_inve(num,Eri,ilab,jlab,klab,llab)
      implicit none
c::passed
      integer num                   ! number of 2e moieties
      integer ilab(num),jlab(num)   ! labels for moieties
      integer klab(num),llab(num)   ! labels for moieties
      double precision Eri(num)     ! 2e moieties
c
      integer ii, i, j, k, l
      double precision x
      logical coul, exch
c
      double precision ene
      integer icall, nused, nout
      save icall, nused, nout
      save ene
      data icall /0/
      data nused /0/
      data nout /0/
c
      if (icall.eq.0) ene = 0.0d00
      icall = icall + 1
c
      do ii = 1,num
        i = ilab(ii)
        j = jlab(ii)
        k = klab(ii)
        l = llab(ii)
        x = Eri(ii)
        coul = (i.eq.j).and.(k.eq.l)
        exch = (i.eq.k).and.(j.eq.l)
        if (coul.and.exch) then
          ene = ene + 2*x
          nused = nused + 1
        elseif (coul.or.exch) then
          ene = ene + x
          nused = nused + 1
        else
          nout = nout + 1
        endif
      enddo
c
      write(6,*)' stsum_inve: icall  =',icall
      write(6,*)' stsum_inve: nused  =',nused
      write(6,*)' stsum_inve: nout   =',nout
      write(6,*)' stsum_inve: n_sum  =',(nused+nout)
      write(6,*)' stsum_inve: energy =',ene
c
      end
      subroutine stsum2(num,numtot,Eri,ilab,jlab,klab,llab,
     &       stats,printem,print69)
      implicit none
c::passed
      integer num                   ! number of 2e moieties
c
c...   numtot must be zeroed prior to first call to stsum2 ...
c
      integer numtot                ! accumulater
c
      integer ilab(num),jlab(num)   ! labels for moieties
      integer klab(num),llab(num)   ! labels for moieties
      double precision Eri(num)     ! 2e moieties
      double precision stats(5,3)   ! stats array (zeroed elsewhere)
      logical printem               ! print values as you go??
      logical print69               ! print values to unit 69??
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
10001     format(i8,4(i4,1x),5x,'eri = ',1pd20.10)
        endif
        if (print69) then
          l4 = isym4(i,j,k,l)
          write(69,10002)l4,i,j,k,l,x,'4 center eri'
10002     format(i8,4(i4,1x),1x,f30.20,1x,a)
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
      subroutine stsum2_3c(num,numtot,Eri,ilab,jlab,klab,
     &       stats,printem,print69)
      implicit none
c::passed
      integer num                   ! number of 2e moieties
c
c...   numtot must be zeroed prior to first call to stsum_3c
c
      integer numtot                ! accumulater
c
      integer ilab(num),jlab(num)   ! labels for moieties
      integer klab(num)             ! labels for moieties
      double precision Eri(num)     ! 2e moieties
      double precision stats(5,3)   ! stats array (zeroed elsewhere)
      logical printem               ! print values as you go??
      logical print69
c::local 
      integer ii,i,j,k
      double precision x,x2,denom,s1,d1
c
      numtot = numtot + num
      do 00100 ii = 1,num
        i     = ilab(ii)
        j     = jlab(ii)
        k     = klab(ii)
        s1    = i+j+k
        d1    = abs(abs(i-j)-abs(k))
        x     = Eri(ii)
        if (printem) then
          write(6,10001)i,j,k,x
10001     format(1x,3(i4,1x),5x,'eri(3c) = ',1pd20.10)
        endif
        if (print69) then
          write(69,10002)i,j,k,x,'3 center eri'
10002     format(3(i4,1x),1x,f30.20,1x,a)
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
      subroutine stsum2_2c(num,numtot,Eri,ilab,jlab,
     &       stats,printem,print69)
      implicit none
c::passed
      integer num                   ! number of 2e moieties
c
c...   numtot must be zeroed prior to first call to stsum2_2c
c
      integer numtot                ! accumulater
c
      integer ilab(num),jlab(num)   ! labels for moieties
      double precision Eri(num)     ! 2e moieties
      double precision stats(5,3)   ! stats array (zeroed elsewhere)
      logical printem               ! print values as you go??
      logical print69
c::local 
      integer ii,i,j
      double precision x,x2,denom,s1,d1
c
      numtot = numtot + num
      do 00100 ii = 1,num
        i     = ilab(ii)
        j     = jlab(ii)
        s1    = i+j
        d1    = abs(abs(i)-abs(j))
        x     = Eri(ii)
        if (printem) then
          write(6,10001)i,j,x
10001    format(1x,2(i4,1x),5x,'eri(2c) = ',1pd20.10)
       endif
        if (print69) then
          write(69,10002)i,j,x,'2 center eri'
10002     format(2(i4,1x),1x,f30.20,1x,a)
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

