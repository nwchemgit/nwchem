      subroutine int_init(nbas, bases)
c $Id: int_init.F,v 1.11 1994-08-17 08:40:52 d3e129 Exp $
c
c initializes integral code to data structers for a integral computation
c
      implicit none
#include "bas.fh"
#include "apiP.fh"
c::passed
      integer nbas        ! number of basis sets to be used
      integer bases(nbas) ! basis set handles
c::local
      integer ibas, ang2use, angm
      logical status
c
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'int_init')
00100 continue

      if (.not.status) then
        write(6,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(6,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('int_init: basis handles hosed ',nbas)
      endif
*      write(6,*)' int_init: basis set handles valid '
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('int_init: angm error',angm)
        ang2use = max(ang2use,angm)
00300 continue
c.      call defNxyz(ang2use)
c.. for derivatives
      call defNxyz(ang2use+1)
c
c... generate memory requirements and store in structures in apiP.fh
c
#if defined(INT_TIMING)
      call int_tim_init()
#endif
      call exact_mem(bases,nbas)
      call sp_init(nbas,bases)
      call int_acc_std()
      init_int = 1
      end
      block data api_data
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
      implicit none
#include "apiP.fh"
c     
      data mem_2e4c  /0/
      data init_int  /0/
      data mem_1e    /0/
      data isz_2e4c  /0/
      data isz_1e    /0/
      data mem_1e_h1 /0/
      data isz_2e3c  /0/
      data isz_2e2c  /0/
      data mem_2e3c  /0/
      data mem_2e2c  /0/
      data isz_3ov   /0/
      data mem_3ov   /0/
      data init_intd /0/
      data DCexp     /0.0000000000D00/
      data DCcoeff   /1.0000000000D00/
      data val_int_acc /0.0d00/
      end
      subroutine sp_init(nbas,bases)
      implicit none
#include "apispP.fh"
#include "bas.fh"
#include "mafdecls.fh"
c::passed
      integer nbas
      integer bases(nbas)
c::local
      integer ncont,icont,ibas
      integer mytype,nprim,ngen,sphorcart
      integer myang, testang
c
      no_sp_conts = .true.
      do 00100 ibas = 1,nbas
        if (.not.bas_numcont(bases(ibas),ncont))
     &         call errquit('sp_init: bas_numcont failed',911)
        do 00200 icont = 1,ncont
          if(.not.bas_continfo
     &           (bases(ibas),icont,mytype,nprim,ngen,sphorcart))
     &           call errquit('sp_init: bas_continfo failed',911)

          if (mytype.lt.0) then
            no_sp_conts = .false.
            goto 00201
          endif
00200   continue
00100 continue
00201 continue
c
      if (no_sp_conts) return
c... arbitrary set max number of blocks to 2 for sp's if spd's set to 3
      max_sp_blksz = 2
c
      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block A',
     &       hszspa,kszspa))
     &       call errquit('sp_init:ma_alloc_get A failed',911)
      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block B',
     &       hszspb,kszspb))
     &       call errquit('sp_init:ma_alloc_get B failed',911)
      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block C',
     &       hszspc,kszspc))
     &       call errquit('sp_init:ma_alloc_get C failed',911)
      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block D',
     &       hszspd,kszspd))
     &       call errquit('sp_init:ma_alloc_get D failed',911)
c
      myang = -56565
      do 00300 ibas = 1,nbas
        if (.not.bas_high_angular(bases(ibas),testang))
     &         call errquit('sp_init:bas_high_angular failed',911)
        myang = max(myang,testang)
00300 continue
      myang = (myang+1)*(myang+2)/2
* max of (sp,sp|sp,sp) and (sp,X,X,X) and permutations
      sp_scat_sz = max((4*4*4*4),(4*4*4*myang))
      sp_scat_sz = max(sp_scat_sz,(4*4*myang*myang))
      sp_scat_sz = max(sp_scat_sz,(4*myang*myang*myang))
      sp_scat_sz = max(sp_scat_sz,(myang*myang*myang*myang))
c
      if (.not.MA_alloc_get(mt_int,sp_scat_sz,'sp scatter array',
     &       hspscat,kspscat))
     &       call errquit('sp_init:ma_alloc_get scat failed',911)
      end
#if defined(INT_TIMING)
      subroutine int_tim_init()
      implicit none
#include "int_tim.fh"
#include "pstat.fh"
c
      if(.not.Pstat_init(num_int_tim,0,' '))
     &    call errquit('int_tim_init: pstat_init failed?',911)
c
      if(.not.Pstat_allocate('hfmkr_1',Pstat_Qstat,0,0,hand_int(1)))
     &    call errquit('int_tim_init:pstat_allocate_1 failed?',911)
      if(.not.Pstat_allocate('hfmkr_2',Pstat_Qstat,0,0,hand_int(2)))
     &    call errquit('int_tim_init:pstat_allocate_2 failed?',911)
      if(.not.Pstat_allocate('hfmkr_3',Pstat_Qstat,0,0,hand_int(3)))
     &    call errquit('int_tim_init:pstat_allocate_3 failed?',911)
      if(.not.Pstat_allocate('hfmkr_4',Pstat_Qstat,0,0,hand_int(4)))
     &    call errquit('int_tim_init:pstat_allocate_4 failed?',911)
      if(.not.Pstat_allocate('hfmkr_5',Pstat_Qstat,0,0,hand_int(5)))
     &    call errquit('int_tim_init:pstat_allocate_5 failed?',911)
      if(.not.Pstat_allocate('hfmkr_6',Pstat_Qstat,0,0,hand_int(6)))
     &    call errquit('int_tim_init:pstat_allocate_6 failed?',911)
      if(.not.Pstat_allocate('hfmkr_7',Pstat_Qstat,0,0,hand_int(7)))
     &    call errquit('int_tim_init:pstat_allocate_7 failed?',911)
      if(.not.Pstat_allocate('hfmkr_8',Pstat_Qstat,0,0,hand_int(8)))
     &    call errquit('int_tim_init:pstat_allocate_8 failed?',911)
      if(.not.Pstat_allocate('hfmkr_9',Pstat_Qstat,0,0,hand_int(9)))
     &    call errquit('int_tim_init:pstat_allocate_9 failed?',911)
      if(.not.Pstat_allocate('hfmkr',Pstat_Qstat,0,0,hand_int(10)))
     &    call errquit('int_tim_init:pstat_allocate_10 failed?',911)
      if(.not.Pstat_allocate('hfmkr all',Pstat_Qstat,0,0,hand_int(11)))
     &    call errquit('int_tim_init:pstat_allocate_11 failed?',911)
      if(.not.Pstat_allocate('gamma',Pstat_Qstat,0,0,hand_int(12)))
     &    call errquit('int_tim_init:pstat_allocate_12 failed?',911)
      if(.not.Pstat_allocate('hf2',Pstat_Qstat,0,0,hand_int(13)))
     &    call errquit('int_tim_init:pstat_allocate_13 failed?',911)
      if(.not.Pstat_allocate('hfset',Pstat_Qstat,0,0,hand_int(14)))
     &    call errquit('int_tim_init:pstat_allocate_14 failed?',911)
      if(.not.Pstat_allocate('hfmke',Pstat_Qstat,0,0,hand_int(15)))
     &    call errquit('int_tim_init:pstat_allocate_15 failed?',911)
      if(.not.Pstat_allocate('hf2mkr',Pstat_Qstat,0,0,hand_int(16)))
     &    call errquit('int_tim_init:pstat_allocate_16 failed?',911)
      if(.not.Pstat_allocate('hferi',Pstat_Qstat,0,0,hand_int(17)))
     &    call errquit('int_tim_init:pstat_allocate_17 failed?',911)
c
      end
      subroutine int_tim_terminate()
      implicit none
#include "int_tim.fh"
#include "pstat.fh"
c
      call pstat_print_all()
c
      if(.not.Pstat_terminate())
     &    call errquit('int_tim_terminate:pstat_terminate failed?',911)
c
      end
#endif
