c $Id$
*
*rak:      subroutine int_chk_sum(rtdb,basisin,print_ints)
*rak:      implicit none
*rak:#include "mafdecls.fh"
*rak:#include "errquit.fh"
*rak:#include "global.fh"
*rak:#include "bas.fh"
*rak:#include "msgids.fh"
*rak:#include "rtdb.fh"
*rak:#include "context.fh"
*rak:c::passed
*rak:      integer rtdb       ! [input] rtdb handle
*rak:      integer basisin    ! [input] basis set handle 
*rak:      logical print_ints ! [input] also print integral and label values 
*rak:c::local
*rak:      logical debug
*rak:      logical MAstatus
*rak:      integer max1bf, maxg, mscr1, mscr2, mscrh1, maxh1
*rak:c
*rak:      integer i,j,k,l, numleft, numtot, lend
*rak:      integer memb      ! memory buffer for integrals (1es and eri)
*rak:      integer mems      ! scratch space needed for integrals
*rak:      integer hmemb, hmems ! memory handles
*rak:      integer imemb, imems ! memory indexes
*rak:c
*rak:      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for labels buffers
*rak:      integer i_klab, h_klab, i_llab, h_llab   ! handles and index for labels buffers
*rak:c
*rak:      integer h_stats, i_stats
*rak:c
*rak:      logical print69
*rak:      logical dumm69
*rak:c
*rak:      integer ncont
*rak:      integer plevel
*rak:c
*rak:      integer me, nproc, idoit
*rak:c      
*rak:      data debug /.false./
*rak:      data plevel /0/
*rak:c 
*rak:c
*rak:c if debugging then up print level
*rak:c
*rak:      if (debug) plevel = 1
*rak:c
*rak:c get number of contractions (shells) ncont
*rak:c
*rak:      if (.not.bas_numcont(basisin,ncont))
*rak:     &    call errquit('int_chk_sum: numcont not obtained',ncont)
*rak:c
*rak:      call int_mem(max1bf,maxg,mscr1,mscr2)
*rak:      memb = max(max1bf,maxg)
*rak:      mems = max(mscr1,mscr2)
*rak:      call int_mem_h1(maxh1, mscrh1)
*rak:      memb = max(memb,maxh1)
*rak:      mems = max(mems,mscrh1)
*rak:      call int_mem_2e3c(maxg,mscr2)
*rak:      memb = max(memb,maxg)
*rak:      mems = max(mems,mscr2)
*rak:      call int_mem_2e2c(maxg,mscr2)
*rak:      memb = max(memb,maxg)
*rak:      mems = max(mems,mscr2)
*rak:c
*rak:c...  pray that ma has been initialized still need verification flag for
*rak:c...       consistency.
*rak:c
*rak:c
*rak:c.... eri/1e buffer
*rak:      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
*rak:     &    call errquit('int_chk_sum: allocation memb failed',memb)
*rak:c.... scratch for integrals
*rak:      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
*rak:     &    call errquit('int_chk_sum: allocation mems failed',mems)
*rak:c...  label buffers
*rak:      if (.not.MA_Push_Get(MT_Int,memb,'i label buffer',h_ilab,i_ilab))
*rak:     &    call errquit('int_chk_sum: allocation ilab failed',memb)
*rak:      if (.not.MA_Push_Get(MT_Int,memb,'j label buffer',h_jlab,i_jlab))
*rak:     &    call errquit('int_chk_sum: allocation jlab failed',memb)
*rak:      if (.not.MA_Push_Get(MT_Int,memb,'k label buffer',h_klab,i_klab))
*rak:     &    call errquit('int_chk_sum: allocation klab failed',memb)
*rak:      if (.not.MA_Push_Get(MT_Int,memb,'l label buffer',h_llab,i_llab))
*rak:     &    call errquit('int_chk_sum: allocation llab failed',memb)
*rak:c...  summary stats 1e and 2e.  i.e., allocate bigger
*rak:      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
*rak:     &    h_stats,i_stats))
*rak:     &    call errquit('int_chk_sum: allocation mems failed',(5*3))
*rak:c
*rak:c zero memory 
*rak:      call ifill(memb,0,    int_mb(i_ilab) ,1)
*rak:      call ifill(memb,0,    int_mb(i_jlab) ,1)
*rak:      call ifill(memb,0,    int_mb(i_klab) ,1)
*rak:      call ifill(memb,0,    int_mb(i_llab) ,1)
*rak:      call dfill(memb,0.0d0,  dbl_mb(imemb)  ,1)
*rak:      call dfill(mems,0.0d0,  dbl_mb(imems)  ,1)
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:c
*rak:c initialize parallelization flags
*rak:c
*rak:      me = ga_nodeid()
*rak:      nproc = ga_nnodes()
*rak:      call ga_sync()
*rak:      if ((me.eq.0).and.(nproc.gt.1).and.print_ints)
*rak:     &    write(6,*)' only integrals on node 0 will be printed. '
*rak:c
*rak:c loop over shells and compute S sums
*rak:c      
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:ov69',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      idoit = (me-1)
*rak:      numtot = 0
*rak:      do 00100 i = 1,ncont
*rak:        do 00200 j = 1,i
*rak:          idoit = idoit + 1
*rak:          if(mod(idoit,nproc).eq.0) then
*rak:*            call int_1eov(basisin,i,basisin,j,
*rak:*     &          mems,dbl_mb(imems),memb,dbl_mb(imemb))
*rak:*            call int_lgen1e(basisin,i,basisin,j,1.0d-08,
*rak:*     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*rak:*     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
*rak:            call int_l1eov(basisin,i,basisin,j,1.0d-08,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*rak:     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
*rak:            if(print_ints.and.me.eq.0)
*rak:     &          call int_pgen1e('overlap',basisin,i,basisin,j,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
*rak:     &          dbl_mb(imemb),plevel)
*rak:            call stsum1(numleft,numtot,
*rak:     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
*rak:     &          dbl_mb(i_stats),(print69.and.(me.eq.0)),
*rak:     &          'overlap')
*rak:          endif
*rak:00200   continue
*rak:00100 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum('overlap',dbl_mb(i_stats),numtot)
*rak:      call ga_sync()
*rak:c
*rak:c
*rak:c loop over shells and compute T sums
*rak:c      
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:ke69',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      idoit = (me-1)
*rak:      numtot = 0
*rak:      do 00101 i = 1,ncont
*rak:        do 00201 j = 1,i
*rak:          idoit = idoit + 1
*rak:          if(mod(idoit,nproc).eq.0) then
*rak:            call int_l1eke(basisin,i,basisin,j,1.0d-08,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*rak:     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
*rak:            if(print_ints.and.me.eq.0)
*rak:     &          call int_pgen1e('kinetic energy',basisin,i,basisin,j,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
*rak:     &          dbl_mb(imemb),plevel)
*rak:            call stsum1(numleft,numtot,
*rak:     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
*rak:     &          dbl_mb(i_stats),(print69.and.(me.eq.0)),
*rak:     &          'kinetic')
*rak:          endif
*rak:00201   continue
*rak:00101 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum('kinetic energy',dbl_mb(i_stats),numtot)
*rak:      call ga_sync()
*rak:c
*rak:c loop over shells and compute V sums
*rak:c      
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:pe69',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      idoit = (me-1)
*rak:      numtot = 0
*rak:      do 00102 i = 1,ncont
*rak:        do 00202 j = 1,i
*rak:          idoit = idoit + 1
*rak:          if(mod(idoit,nproc).eq.0) then
*rak:            call int_l1epe(basisin,i,basisin,j,1.0d-08,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*rak:     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
*rak:            if(print_ints.and.me.eq.0)
*rak:     &          call int_pgen1e('pot. energy',basisin,i,basisin,j,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
*rak:     &          dbl_mb(imemb),plevel)
*rak:            call stsum1(numleft,numtot,
*rak:     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
*rak:     &          dbl_mb(i_stats),(print69.and.(me.eq.0))
*rak:     &          ,'potential')
*rak:          endif
*rak:00202   continue
*rak:00102 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum('pot. energy',dbl_mb(i_stats),numtot)
*rak:      call ga_sync()
*rak:c
*rak:c loop over shells and compute h1 sums
*rak:c      
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:      idoit = (me-1)
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:h169',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      numtot = 0
*rak:      do 00103 i = 1,ncont
*rak:        do 00203 j = 1,i
*rak:          idoit = idoit + 1
*rak:          if(mod(idoit,nproc).eq.0) then
*rak:            call int_l1eh1(basisin,i,basisin,j,1.0d-08,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),memb,
*rak:     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
*rak:            if(print_ints.and.me.eq.0)
*rak:     &          call int_pgen1e('h1',basisin,i,basisin,j,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),numleft,
*rak:     &          dbl_mb(imemb),plevel)
*rak:            call stsum1(numleft,numtot,
*rak:     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
*rak:     &          dbl_mb(i_stats),(print69.and.(me.eq.0)),
*rak:     &          ' H1 (T+V)')
*rak:          endif
*rak:00203   continue
*rak:00103 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum('h1',dbl_mb(i_stats),numtot)
*rak:      call ga_sync()
*rak:c
*rak:c do canonical shell list for 2e integrals (4 center)
*rak:c
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:4c69',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      idoit = (me-1)
*rak:      numtot = 0
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:      do 01000 i = 1,ncont
*rak:        if((me.eq.0).and.(mod(i,10).eq.0)) then
*rak:          write(6,*)' 2e4c outer loop counter ',i,' of ',ncont
*rak:          call util_flush(6)
*rak:        endif
*rak:        do 02000 j = 1, i
*rak:          do 03000 k = 1, i
*rak:            lend = k
*rak:            if (i.eq.k) lend = j
*rak:            do 04000 l = 1, lend
*rak:              idoit = idoit + 1
*rak:              if(mod(idoit,nproc).eq.0) then
*rak:                call int_l2e4c(basisin,i,j,basisin,k,l,
*rak:c     &              1.0d-08,.true.,memb,dbl_mb(imemb),numleft,
*rak:     &              1.0d-08,.false.,memb,dbl_mb(imemb),numleft,
*rak:     &              int_mb(i_ilab),int_mb(i_jlab),
*rak:     &              int_mb(i_klab),int_mb(i_llab),
*rak:     &              mems,dbl_mb(imems))
*rak:                call stsum2(numleft,numtot,
*rak:     &              dbl_mb(imemb),
*rak:     &              int_mb(i_ilab),int_mb(i_jlab),
*rak:     &              int_mb(i_klab),int_mb(i_llab),
*rak:     &              dbl_mb(i_stats),
*rak:     &              (print_ints.and.me.eq.0),(print69.and.me.eq.0))
*rak:              endif
*rak:04000       continue
*rak:03000     continue
*rak:02000   continue
*rak:01000 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum('two electron',dbl_mb(i_stats),numtot)      
*rak:      call ga_sync()
*rak:c
*rak:c do canonical shell list for 2e integrals (3 center)
*rak:c
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:3c69',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      idoit = (me-1)
*rak:      numtot = 0
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:      do 01001 i = 1,ncont
*rak:        if((me.eq.0).and.(mod(i,10).eq.0)) then
*rak:          write(6,*)' 2e3c outer loop counter ',i,' of ',ncont
*rak:          call util_flush(6)
*rak:        endif
*rak:        do 02001 j = 1, i
*rak:          do 03001 k = 1, j
*rak:            idoit = idoit + 1
*rak:            if(mod(idoit,nproc).eq.0) then
*rak:              call int_l2e3c(basisin,i,basisin,j,k,
*rak:     &            1.0d-08,.true.,memb,dbl_mb(imemb),numleft,
*rak:     &            int_mb(i_ilab),int_mb(i_jlab),
*rak:     &            int_mb(i_klab),
*rak:     &            mems,dbl_mb(imems))
*rak:              call stsum2_3c(numleft,numtot,
*rak:     &            dbl_mb(imemb),
*rak:     &            int_mb(i_ilab),int_mb(i_jlab),
*rak:     &            int_mb(i_klab),
*rak:     &            dbl_mb(i_stats),(print_ints.and.me.eq.0),
*rak:     &            (print69.and.(me.eq.0)))
*rak:            endif
*rak:03001     continue
*rak:02001   continue
*rak:01001 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum
*rak:     &    ('3 center two electron',dbl_mb(i_stats),numtot)      
*rak:      call ga_sync()
*rak:c
*rak:c do canonical shell list for 2e integrals (2 center)
*rak:c
*rak:      print69 = .false.
*rak:      if (rtdb_get(rtdb,'int:2c69',MT_LOG,1,dumm69)) then
*rak:        print69 = dumm69
*rak:      endif
*rak:      idoit = (me-1)
*rak:      numtot = 0
*rak:      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
*rak:      do 01002 i = 1,ncont
*rak:        if((me.eq.0).and.(mod(i,10).eq.0)) then
*rak:          write(6,*)' 2e2c outer loop counter ',i,' of ',ncont
*rak:          call util_flush(6)
*rak:        endif
*rak:        do 02002 j = 1, i
*rak:          idoit = idoit + 1
*rak:          if(mod(idoit,nproc).eq.0) then
*rak:            call int_l2e2c(basisin,i,basisin,j,
*rak:     &          1.0d-08,memb,dbl_mb(imemb),numleft,
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),
*rak:     &          mems,dbl_mb(imems))
*rak:            call stsum2_2c(numleft,numtot,
*rak:     &          dbl_mb(imemb),
*rak:     &          int_mb(i_ilab),int_mb(i_jlab),
*rak:     &          dbl_mb(i_stats),(print_ints.and.me.eq.0),
*rak:     &            (print69.and.(me.eq.0)))
*rak:          endif
*rak:02002   continue
*rak:01002 continue
*rak:      call ga_igop(Msg_IChkI,numtot,1,'+')
*rak:      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
*rak:      if (me.eq.0)
*rak:     &    call prt_stat_sum
*rak:     &    ('2 center two electron',dbl_mb(i_stats),numtot)      
*rak:      call ga_sync()
*rak:      call ga_sync()
*rak:      call ga_sync()
*rak:c
*rak:c free memory in order
*rak:c
*rak:      MAstatus = MA_Pop_Stack(h_stats)
*rak:      MAstatus = MAstatus .and. MA_Pop_Stack(h_llab)
*rak:      MAstatus = MAstatus .and. MA_Pop_Stack(h_klab)
*rak:      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
*rak:      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
*rak:      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
*rak:      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
*rak:      if (.not.MAstatus)
*rak:     &       call errquit('int_chk_sum: error poping memory stack',0)
*rak:      call ga_sync()
*rak:      call ga_sync()
*rak:      end
      subroutine prt_stat_sum(msg,stats,numtot)
      implicit none
#include "inp.fh"
#include "stdio.fh"
c::passed
      character*(*) msg
      integer numtot
      double precision stats(5,3)
c::local
      integer lmsg
c
      lmsg = inp_strlen(msg)
      write(luout,10000)
      write(luout,*)' summary statistics for ',
     &    msg(1:lmsg),' integrals'
      write(luout,*)' total number of integrals processed: ',
     &    numtot
      write(luout,10001)
      write(luout,10002)stats
      write(luout,10001)
      call util_flush(luout)
10000 format(/,90('-'))
10001 format(90('-'))
10002 format(5(1pd15.8,3x))
      end      
!$HP$ OPTIMIZE OFF
      subroutine stsum1(num,numtot,Gen,ilab,jlab,stats,
     &    printem,print69,tag69)
      implicit none
#include "stdio.fh"
#include "inp.fh"
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
      logical printem               ! print to stdout
      logical print69               ! print to file 69 
*...................................!      (opened by calling routine)
      character*(*) tag69
      integer ltag69
c::local 
      integer ii,i,j
      double precision x,x2,denom,s1,d1
c
      ltag69 = inp_strlen(tag69)
      numtot = numtot + num
      do 00100 ii = 1,num
        i     = ilab(ii)
        j     = jlab(ii)
        s1    = i+j
        d1    = abs(i-j)
        x     = Gen(ii)
        x2    = x*x
        if(printem) then
          write(luout,'(1x,i5,i5,d20.10,3x,a)') i, j, x,
     &        tag69(1:ltag69)
        endif
        if(print69) then
          write(69,'(1x,i5,i5,d20.10,3x,a)') i, j, x, tag69(1:ltag69)
        endif
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
      if (num.ne.-1) then
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
      else
        write(6,*)' stsum_inve: icall  =',icall
        write(6,*)' stsum_inve: nused  =',nused
        write(6,*)' stsum_inve: nout   =',nout
        write(6,*)' stsum_inve: n_sum  =',(nused+nout)
        write(6,*)' stsum_inve: energy =',ene
      endif
c
      end
      subroutine stsum2(num,numtot,Eri,ilab,jlab,klab,llab,
     &       stats,printem,print69,tagout)
      implicit none
#include "stdio.fh"
#include "inp.fh"
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
      logical printem               ! print values as you go?
      logical print69               ! print values to unit 69?
      character*(*) tagout
c::local 
      integer ltagout
      integer ii,i,j,k,l,l4
      double precision x,x2,denom,s1,d1
c
      integer isym2,isym4
      ISYM2(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
      ISYM4(I,J,K,L)=MAX(ISYM2(I,J),ISYM2(K,L))*
     &               (MAX(ISYM2(I,J),ISYM2(K,L))-1)/2+
     &               MIN(ISYM2(I,J),ISYM2(K,L))
c
      ltagout = inp_strlen(tagout)
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
          write(6,10001)l4,i,j,k,l,x,tagout(1:ltagout)
10001     format(i8,4(i4,1x),5x,1pd20.10,1x,a)
        endif
        if (print69) then
          l4 = isym4(i,j,k,l)
          write(69,10002)l4,i,j,k,l,x,tagout(1:ltagout)
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
      subroutine stsum3c(num,numtot,Eri,ilab,jlab,klab,
     &       stats,printem,print69,tagout)
      implicit none
#include "stdio.fh"
#include "inp.fh"
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
      logical printem               ! print values as you go?
      logical print69
      character*(*) tagout
c::local 
      integer ltagout
*rak:      integer l31, l32
      integer ii,i,j,k
      double precision x,x2,denom,s1,d1
c
*rak:      integer isym2,isym3
*rak:      ISYM2(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*rak:      ISYM3(I,J,K)=MAX(ISYM2(I,J),K)*
*rak:     &               (MAX(ISYM2(I,J),K)-1)/2+
*rak:     &               MIN(ISYM2(I,J),K)
c
      ltagout = inp_strlen(tagout)
      numtot = numtot + num
      do 00100 ii = 1,num
        i     = ilab(ii)
        j     = jlab(ii)
        k     = klab(ii)
        s1    = i+j+k
        d1    = abs(abs(i-j)-abs(k))
        x     = Eri(ii)
        if (printem) then
*rak:          l31 = isym3(i,j,k)
*rak:          l32 = isym3(j,k,i)
*rak:          write(6,10001)l31,l32,i,j,k,x,tagout(1:ltagout)
*rak:10001     format(1x,2(i8,1x),3(i4,1x),5x,1pd20.10,1x,a)
          write(6,10001)i,j,k,x,tagout(1:ltagout)
10001     format(1x,3(i4,1x),5x,1pd20.10,1x,a)
        endif
        if (print69) then
*rak:          l31 = isym3(i,j,k)
*rak:          l32 = isym3(j,k,i)
*rak:          write(69,10002)l31,l32,i,j,k,x,tagout(1:ltagout)
*rak:10002     format(2(i8,1x),3(i4,1x),1x,f30.20,1x,a)
          write(69,10002)i,j,k,x,tagout(1:ltagout)
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
!$HP$ OPTIMIZE ON
*rak:      subroutine stsum2_2c(num,numtot,Eri,ilab,jlab,
*rak:     &       stats,printem,print69)
*rak:      implicit none
*rak:c::passed
*rak:      integer num                   ! number of 2e moieties
*rak:c
*rak:c...   numtot must be zeroed prior to first call to stsum2_2c
*rak:c
*rak:      integer numtot                ! accumulater
*rak:c
*rak:      integer ilab(num),jlab(num)   ! labels for moieties
*rak:      double precision Eri(num)     ! 2e moieties
*rak:      double precision stats(5,3)   ! stats array (zeroed elsewhere)
*rak:      logical printem               ! print values as you go?
*rak:      logical print69
*rak:c::local 
*rak:      integer ii,i,j
*rak:      double precision x,x2,denom,s1,d1
*rak:c
*rak:      numtot = numtot + num
*rak:      do 00100 ii = 1,num
*rak:        i     = ilab(ii)
*rak:        j     = jlab(ii)
*rak:        s1    = i+j
*rak:        d1    = abs(i-j)
*rak:        x     = Eri(ii)
*rak:        if (printem) then
*rak:          write(6,10001)i,j,x
*rak:10001    format(1x,2(i4,1x),5x,'eri(2c) = ',1pd20.10)
*rak:       endif
*rak:        if (print69) then
*rak:          write(69,10002)i,j,x,'2 center eri'
*rak:10002     format(2(i4,1x),1x,f30.20,1x,a)
*rak:        endif
*rak:        x2    = x*x
*rak:        denom = 1.0d00 + x2
*rak:c
*rak:        stats(1,1) = stats(1,1) + x
*rak:        stats(2,1) = stats(2,1) + abs(x)
*rak:        stats(3,1) = stats(3,1) + x2
*rak:        stats(4,1) = stats(4,1) + x/denom
*rak:        stats(5,1) = stats(5,1) + x2/denom
*rak:c
*rak:        stats(1,2) = stats(1,2) + x        * s1
*rak:        stats(2,2) = stats(2,2) + abs(x)   * s1
*rak:        stats(3,2) = stats(3,2) + x2       * s1
*rak:        stats(4,2) = stats(4,2) + x/denom  * s1
*rak:        stats(5,2) = stats(5,2) + x2/denom * s1
*rak:c
*rak:        stats(1,3) = stats(1,3) + x        * d1
*rak:        stats(2,3) = stats(2,3) + abs(x)   * d1
*rak:        stats(3,3) = stats(3,3) + x2       * d1
*rak:        stats(4,3) = stats(4,3) + x/denom  * d1
*rak:        stats(5,3) = stats(5,3) + x2/denom * d1
*rak:c
*rak:00100 continue
*rak:      end
*
      subroutine int_chk_sum_ov(rtdb,basisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer max1bf, maxg, mscr1, mscr2
c
      integer i, j, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers
c
      integer h_stats, i_stats
c
      character*255 name69
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
     &    call errquit('int_chk_sum_ov: numcont not obtained',ncont,
     &           INT_ERR)
c
      call int_mem(max1bf,maxg,mscr1,mscr2)
      memb = max1bf
      mems = mscr1
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit('int_chk_sum_ov: allocation memb failed',memb,
     &            MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit('int_chk_sum_ov: allocation mems failed',mems,
     &           MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit('int_chk_sum_ov: allocation ilab failed',memb,
     &           MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit('int_chk_sum_ov: allocation jlab failed',memb,
     &          MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_ov: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
      call dfill(memb,0.0d0,  dbl_mb(imemb)  ,1)
      call dfill(mems,0.0d0,  dbl_mb(imems)  ,1)
      call dfill((5*3),0.0d0, dbl_mb(i_stats),1)
c
c initialize parallelization flags
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      if ((me.eq.0).and.(nproc.gt.1).and.print_ints)
     &    write(6,*)' only integrals on node 0 will be printed. '
c
      call ga_sync()
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:ovprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.ov',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute S sums
c      
      idoit = (me-1)
      numtot = 0
      do 00100 i = 1,ncont
        do 00200 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
            call int_l1eov(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),
     &          (print_ints.and.(me.eq.0)),
     &          print69,
     &          'overlap')
          endif
00200   continue
00100 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('overlap',dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_ov: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end
      subroutine int_chk_sum_ke(rtdb,basisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer max1bf, maxg, mscr1, mscr2
c
      integer i, j, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers
c
      integer h_stats, i_stats
c
      character*255 name69
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
     &    call errquit('int_chk_sum_ke: numcont not obtained',ncont,
     &           INT_ERR)
c
      call int_mem(max1bf,maxg,mscr1,mscr2)
      memb = max1bf
      mems = mscr1
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit('int_chk_sum_ke: allocation memb failed',memb,
     &           MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit('int_chk_sum_ke: allocation mems failed',mems,
     &             MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit('int_chk_sum_ke: allocation ilab failed',memb,
     &            MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit('int_chk_sum_ke: allocation jlab failed',memb,
     &            MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_ke: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:keprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.ke',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute T sums
c      
      idoit = (me-1)
      numtot = 0
      do 00100 i = 1,ncont
        do 00200 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
            call int_l1eke(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),
     &          (print_ints.and.(me.eq.0)),
     &          print69,
     &          'kinetic')
          endif
00200   continue
00100 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('kinetic',dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_ke: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end
      subroutine int_chk_sum_pe(rtdb,basisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer max1bf, maxg, mscr1, mscr2
c
      integer i, j, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers
c
      integer h_stats, i_stats
c
      character*255 name69
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
     &    call errquit('int_chk_sum_pe: numcont not obtained',ncont,
     &             INT_ERR)
c
      call int_mem(max1bf,maxg,mscr1,mscr2)
      memb = max1bf
      mems = mscr1
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit('int_chk_sum_pe: allocation memb failed',memb,
     &           MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit('int_chk_sum_pe: allocation mems failed',mems,
     &            MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit('int_chk_sum_pe: allocation ilab failed',memb,
     &            MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit('int_chk_sum_pe: allocation jlab failed',memb,
     &           MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_pe: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:peprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.pe',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute V sums
c      
      idoit = (me-1)
      numtot = 0
      do 00100 i = 1,ncont
        do 00200 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
            call int_l1epe(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),
     &          (print_ints.and.(me.eq.0)),
     &          print69,
     &          'potential')
          endif
00200   continue
00100 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('potential',dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_pe: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end
      subroutine int_chk_sum_h1(rtdb,basisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer maxh1, mscrh1
c
      integer i, j, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers
c
      integer h_stats, i_stats
c
      character*255 name69
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
     &    call errquit('int_chk_sum_h1: numcont not obtained',ncont,
     &          INT_ERR)
c
      call int_mem_h1(maxh1,mscrh1)
      memb = maxh1
      mems = mscrh1
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit('int_chk_sum_h1: allocation memb failed',memb,
     &           MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit('int_chk_sum_h1: allocation mems failed',mems,
     &            MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit('int_chk_sum_h1: allocation ilab failed',memb,
     &             MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit('int_chk_sum_h1: allocation jlab failed',memb,
     &             MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_h1: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:h1print',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.h1',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute H1 sums
c      
      idoit = (me-1)
      numtot = 0
      do 00100 i = 1,ncont
        do 00200 j = 1,i
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
            call int_l1eh1(basisin,i,basisin,j,1.0d-08,
     &          int_mb(i_ilab),int_mb(i_jlab),memb,
     &          dbl_mb(imemb),mems,dbl_mb(imems),numleft)
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),int_mb(i_ilab),int_mb(i_jlab),
     &          dbl_mb(i_stats),
     &          (print_ints.and.(me.eq.0)),
     &          print69,
     &          'H1 (T + V)')
          endif
00200   continue
00100 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('H1 (T + V)',dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_h1: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end
      subroutine int_chk_sum_3ov(rtdb,basisin,ribasisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
#include "stdio.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      integer ribasisin  ! [input] basis set handle for fitting basis
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer max3ov, mscr3ov
c
      integer i, j, k, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers
      integer i_klab, h_klab
c
      integer h_stats, i_stats
c
      character*11 tagextra
      character*40 tag69
c
      character*255 name69
      logical print69
      logical dumm69
c
      integer ncont, rincont
      integer plevel
      integer ribasisuse
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
     &    call errquit('int_chk_sum_3ov: numcont not obtained',ncont,
     &            INT_ERR)
c
c  if (ribasisin = -1) then no ri basis is loaded so use 
c      regular ao basis
c
      if (ribasisin.eq.-1) then
        ribasisuse = basisin
        tagextra = '<ao,ao,ao> '
*        write(luout,*)
*     &      '3 center overlap integrals computed <ao,ao,ao>'
      else
        ribasisuse = ribasisin
        tagextra = '<ao,ao,fit>'
*        write(luout,*)
*     &      '3 center overlap integrals computed <ao,ao,fit>'
      endif
      tag69 = '3 center overlap'//' '//tagextra
c
      if (.not.bas_numcont(ribasisuse,rincont))
     &    call errquit
     &    ('int_chk_sum_3ov: rincont not obtained',rincont, INT_ERR)
c
      call int_mem_3ov(max3ov,mscr3ov)
      memb = max3ov
      mems = mscr3ov
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit('int_chk_sum_3ov: allocation memb failed',memb,
     &                MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit('int_chk_sum_3ov: allocation mems failed',mems,
     &           MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit('int_chk_sum_3ov: allocation ilab failed',memb,
     &               MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit('int_chk_sum_3ov: allocation jlab failed',memb, 
     &               MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'k label buffer',h_klab,i_klab))
     &    call errquit('int_chk_sum_3ov: allocation klab failed',memb,
     &               MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_3ov: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
      call ifill(memb,0,    int_mb(i_klab) ,1)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:3ovprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.3ov',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute 3 center overlap sums
c      
      idoit = (me-1)
      numtot = 0
      do i = 1,ncont
        do j = 1,i
          do k = 1,rincont
            idoit = idoit + 1
            if(mod(idoit,nproc).eq.0) then
              call int_l1e3ov(
     &            basisin,i,
     &            basisin,j,
     &            ribasisuse,k,
     &            1.0d-08,
     &            memb,dbl_mb(imemb),
     &            int_mb(i_ilab),
     &            int_mb(i_jlab),
     &            int_mb(i_klab),
     &            numleft,mems,dbl_mb(imems))
              call stsum3c(numleft,numtot,
     &            dbl_mb(imemb),
     &            int_mb(i_ilab),
     &            int_mb(i_jlab),
     &            int_mb(i_klab),
     &            dbl_mb(i_stats),
     &            (print_ints.and.(me.eq.0)),
     &            print69,
     &            tag69)
            endif
          enddo
        enddo
      enddo
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum(tag69,dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_klab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_3ov: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end
      subroutine int_chk_sum_2e4c(rtdb,basisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
#include "stdio.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      logical print_ints ! [input] also print integral and label 
*........................!         values to stdout.
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
      integer i_ilab, h_ilab ! handles and index for labels buffers
      integer i_jlab, h_jlab ! handles and index for labels buffers
      integer i_klab, h_klab ! handles and index for labels buffers
      integer i_llab, h_llab ! handles and index for labels buffers
c
      integer h_stats, i_stats
c
      character*255 name69
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
     &    call errquit('int_chk_sum_2e4c: numcont not obtained',ncont, 
     &           INT_ERR)
c
      call int_mem(max1bf,maxg,mscr1,mscr2)
      memb = maxg
      mems = mscr2
c
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation memb failed',memb, MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation mems failed',mems, MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,'i label buffer',
     &    h_ilab,i_ilab))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation ilab failed',memb, MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,'j label buffer',
     &    h_jlab,i_jlab))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation jlab failed',memb, MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,'k label buffer',
     &    h_klab,i_klab))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation klab failed',memb, MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,'l label buffer',
     &    h_llab,i_llab))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation llab failed',memb, MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_2e4c: allocation stats failed',(5*3), MEM_ERR)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:2e4cprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.2e4c',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c do canonical shell list for 2e integrals (4 center)
c
      idoit = (me-1)
      numtot = 0
      do 01000 i = 1,ncont
        if((me.eq.0).and.(mod(i,10).eq.0)) then
          write(luout,*)' 2e4c outer loop counter ',i,' of ',ncont
          call util_flush(luout)
        endif
        do 02000 j = 1, i
          do 03000 k = 1, i
            lend = k
            if (i.eq.k) lend = j
            do 04000 l = 1, lend
              idoit = idoit + 1
              if(mod(idoit,nproc).eq.0) then
                call int_l2e4c(basisin,i,j,basisin,k,l,
     &              1.0d-08,.false.,memb,dbl_mb(imemb),numleft,
     &              int_mb(i_ilab),int_mb(i_jlab),
     &              int_mb(i_klab),int_mb(i_llab),
     &              mems,dbl_mb(imems))
                call stsum2(numleft,numtot,
     &              dbl_mb(imemb),
     &              int_mb(i_ilab),int_mb(i_jlab),
     &              int_mb(i_klab),int_mb(i_llab),
     &              dbl_mb(i_stats),
     &              (print_ints.and.(me.eq.0)),
     &              print69,'two electron (2e4c)')
              endif
04000       continue
03000     continue
02000   continue
01000 continue
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum('two electron (2e4c)',
     &    dbl_mb(i_stats),numtot)      
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
     &       call errquit
     &    ('int_chk_sum_2e4c: error poping memory stack',0, MEM_ERR)
* close file      
      close(unit=69,status='keep')
      call ga_sync()
      end
      subroutine int_chk_sum_2e3c(rtdb,basisin,ribasisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
#include "stdio.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      integer ribasisin  ! [input] basis set handle for fitting basis
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer max2e3c, mscr2e3c
c
      integer i, j, k, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers
      integer i_klab, h_klab
c
      integer h_stats, i_stats
c
      character*11 tagextra
      character*40 tag69
c
      character*255 name69
      logical print69
      logical dumm69
c
      integer ncont, rincont
      integer plevel
      integer ribasisuse
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
     &    call errquit('int_chk_sum_2e3c: numcont not obtained',ncont,
     &          MEM_ERR)
c
c  if (ribasisin = -1) then no ri basis is loaded so use 
c      regular ao basis
c
      if (ribasisin.eq.-1) then
        ribasisuse = basisin
        tagextra = '<ao|ao,ao> '
*        write(luout,*)
*     &      '2 electron 3 center ERI integrals computed <ao,ao,ao>'
      else
        ribasisuse = ribasisin
        tagextra = '<fit|ao,ao>'
*        write(luout,*)
*     &      '2 electron 3 center ERI integrals computed <fit,ao,ao>'
      endif
      tag69 = '2e-3 center ERI'//' '//tagextra
c
      if (.not.bas_numcont(ribasisuse,rincont))
     &    call errquit
     &    ('int_chk_sum_2e3c: rincont not obtained',rincont, MEM_ERR)
c
      call int_mem_2e3c(max2e3c,mscr2e3c)
      memb = max2e3c
      mems = mscr2e3c
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit
     &    ('int_chk_sum_2e3c: allocation memb failed',memb, MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit
     &    ('int_chk_sum_2e3c: allocation mems failed',mems, MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit
     &    ('int_chk_sum_2e3c: allocation ilab failed',memb, MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit
     &    ('int_chk_sum_2e3c: allocation jlab failed',memb, MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'k label buffer',h_klab,i_klab))
     &    call errquit
     &    ('int_chk_sum_2e3c: allocation klab failed',memb, MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_2e3c: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
      call ifill(memb,0,    int_mb(i_klab) ,1)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:2e3cprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.2e3c',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute 3 center overlap sums
c      
      idoit = (me-1)
      numtot = 0
      do i = 1,rincont
        do j = 1,ncont
          do k = 1,j
            idoit = idoit + 1
            if(mod(idoit,nproc).eq.0) then
              call int_l2e3c(
     &            ribasisuse,i,
     &            basisin,j,k,
     &            1.0d-08,.false.,
     &            memb,dbl_mb(imemb),
     &            numleft,
     &            int_mb(i_ilab),
     &            int_mb(i_jlab),
     &            int_mb(i_klab),
     &            mems,dbl_mb(imems))
              call stsum3c(numleft,numtot,
     &            dbl_mb(imemb),
     &            int_mb(i_ilab),
     &            int_mb(i_jlab),
     &            int_mb(i_klab),
     &            dbl_mb(i_stats),
     &            (print_ints.and.(me.eq.0)),
     &            print69,
     &            tag69)
            endif
          enddo
        enddo
      enddo
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum(tag69,dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_klab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_2e3c: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end
      subroutine int_chk_sum_2e2c(rtdb,basisin,ribasisin,print_ints)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "context.fh"
#include "stdio.fh"
c::passed
      integer rtdb       ! [input] rtdb handle
      integer basisin    ! [input] basis set handle 
      integer ribasisin  ! [input] basis set handle for fitting basis
      logical print_ints ! [input] also print integral/label values 
c::local
      logical debug
      logical MAstatus
      integer max2e2c, mscr2e2c
c
      integer i, j, numleft, numtot
      integer memb      ! memory buffer for integrals (1es and eri)
      integer mems      ! scratch space needed for integrals
      integer hmemb, hmems ! memory handles
      integer imemb, imems ! memory indexes
c
      integer i_ilab, h_ilab, i_jlab, h_jlab   ! handles and index for
*                                              ! labels and buffers

      integer h_stats, i_stats
c
      character*8 tagextra
      character*40 tag69
c
      character*255 name69
      logical print69
      logical dumm69
c
      integer ncont, rincont
      integer plevel
      integer ribasisuse
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
     &    call errquit('int_chk_sum_2e2c: numcont not obtained',ncont,
     &           INT_ERR)
c
c  if (ribasisin = -1) then no ri basis is loaded so use 
c      regular ao basis
c
      if (ribasisin.eq.-1) then
        ribasisuse = basisin
c                   12345678901
        tagextra = '<ao|ao> '
*        write(luout,*)
*     &      '2 electron 2 center ERI integrals computed <ao|ao>'
      else
        ribasisuse = ribasisin
        tagextra = '<fit|ao>'
*        write(luout,*)
*     &      '2 electron 2 center ERI integrals computed <fit|ao>'
      endif
      tag69 = '2e-2 center ERI'//' '//tagextra
c
      if (.not.bas_numcont(ribasisuse,rincont))
     &    call errquit
     &    ('int_chk_sum_2e2c: rincont not obtained',rincont, INT_ERR)
c
      call int_mem_2e2c(max2e2c,mscr2e2c)
      memb = max2e2c
      mems = mscr2e2c
c
c.... eri/1e buffer
      if (.not.MA_Push_Get(MT_Dbl,memb,'int buffer',hmemb,imemb))
     &    call errquit
     &    ('int_chk_sum_2e2c: allocation memb failed',memb, MEM_ERR)
c.... scratch for integrals
      if (.not.MA_Push_Get(MT_Dbl,mems,'int scr buffer',hmems,imems))
     &    call errquit
     &    ('int_chk_sum_2e2c: allocation mems failed',mems, MEM_ERR)
c...  label buffers
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'i label buffer',h_ilab,i_ilab))
     &    call errquit
     &    ('int_chk_sum_2e2c: allocation ilab failed',memb, MEM_ERR)
      if (.not.MA_Push_Get(MT_Int,memb,
     &    'j label buffer',h_jlab,i_jlab))
     &    call errquit
     &    ('int_chk_sum_2e2c: allocation jlab failed',memb, MEM_ERR)
c...  summary stats 1e and 2e.  i.e., allocate bigger
      if (.not.MA_Push_Get(MT_Dbl,(5*3),'summary stats',
     &    h_stats,i_stats))
     &    call errquit
     &    ('int_chk_sum_2e2c: allocation stats failed',(5*3), MEM_ERR)
c
c zero memory 
      call ifill(memb,0,    int_mb(i_ilab) ,1)
      call ifill(memb,0,    int_mb(i_jlab) ,1)
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
      print69 = .false.
      if (rtdb_get(rtdb,'intcsum:2e2cprint',MT_LOG,1,dumm69)) then
        print69 = dumm69
      endif
*
*-- open files for each node
*
      if (print69) then
        call util_file_name('csum.2e2c',.false.,.true.,name69)
        open(unit=69,file=name69,
     &      form='formatted',access='sequential',
     &      status='unknown')
      endif
c
c loop over shells and compute 3 center overlap sums
c      
      idoit = (me-1)
      numtot = 0
      do i = 1,rincont
        do j = 1,ncont
          idoit = idoit + 1
          if(mod(idoit,nproc).eq.0) then
            call int_l2e2c(
     &          ribasisuse,i,
     &          basisin,j,
     &          1.0d-08,
     &          memb,dbl_mb(imemb),
     &          numleft,
     &          int_mb(i_ilab),
     &          int_mb(i_jlab),
     &          mems,dbl_mb(imems))
            call stsum1(numleft,numtot,
     &          dbl_mb(imemb),
     &          int_mb(i_ilab),
     &          int_mb(i_jlab),
     &          dbl_mb(i_stats),
     &          (print_ints.and.(me.eq.0)),
     &          print69,
     &          tag69)
          endif
        enddo
      enddo
      call ga_igop(Msg_IChkI,numtot,1,'+')
      call ga_dgop(Msg_IChkD,dbl_mb(i_stats),(5*3),'+')
      if (me.eq.0)
     &    call prt_stat_sum(tag69,dbl_mb(i_stats),numtot)
c
c free memory in order
c
      MAstatus = MA_Pop_Stack(h_stats)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_jlab)
      MAstatus = MAstatus .and. MA_Pop_Stack(h_ilab)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmems)
      MAstatus = MAstatus .and. MA_Pop_Stack(hmemb)
      if (.not.MAstatus)
     &    call errquit
     &    ('int_chk_sum_2e2c: error poping memory stack',0, MEM_ERR)
c
* close file      
      close(unit=69,status='keep')
c      
      call ga_sync()
      end

