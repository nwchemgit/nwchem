      subroutine int_init(rtdb, nbas, bases)
c $Id: int_init.F,v 1.24 1996-07-11 22:21:15 d3e129 Exp $
c
c initializes integral code to data structers for a integral computation
c
      implicit none
#include "bas.fh"
#include "apiP.fh"
#include "candoP.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c::passed
      integer rtdb        ! run time data base handle
      integer nbas        ! number of basis sets to be used
      integer bases(nbas) ! basis set handles
c::local
      integer txs_mem_min ! memory from texas
      integer ibas, ang2use, angm
      logical status
c
c      block data api_data
c
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
c     
      mem_2e4c  = 0
      init_int  = 0
      mem_1e    = 0
      isz_2e4c  = 0
      isz_1e    = 0
      mem_1e_h1 = 0
      isz_2e3c  = 0
      isz_2e2c  = 0
      mem_2e3c  = 0
      mem_2e2c  = 0
      isz_3ov   = 0
      mem_3ov   = 0
      iszb_2e4c = 0
      memb_2e4c = 0
      init_intd = 0
      DCexp     = 0.0D00
      DCcoeff   = 1.0D00
      val_int_acc = 0.0d00
c
      intd_memthresh = 0
      numd_tot       = 0
      numd_okay      = 0
      numd_red       = 0
c
      if(init_int.eq.1) then
        write(6,*)' warning nested int_inits'
        write(6,*)' int_init already called '
        call util_flush(6)
      endif
c
c initialize cando information from rtdb
c
      user_cando_sp   = .false.
      user_cando_nw   = .false.
      user_cando_txs  = .false.
      def_cando_sp    = .false.
      def_cando_nw    = .false.
      def_cando_txs   = .false.
c
      if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,status)) then
        user_cando_sp = .true.
        def_cando_sp  = status
        if (ga_nodeid().eq.0) then
          write(6,*)
     &        ' int_init: cando_sp set to always be ',def_cando_sp
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
        user_cando_nw = .true.
        def_cando_nw  = status
        if (ga_nodeid().eq.0) then
          write(6,*)
     &        ' int_init: cando_nw set to always be ',def_cando_nw
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
        user_cando_txs = .true.
        def_cando_txs  = status
        if (ga_nodeid().eq.0) then
          write(6,*)
     &        ' int_init: cando_txs set to always be ',def_cando_txs
          call util_flush(6)
        endif
      endif
* sanity checking: e.g., you only want to turn off a particular integral
* code never always turn it on.
*
      if (def_cando_sp.or.def_cando_nw.or.def_cando_txs) then
        if (ga_nodeid().eq.0) then
          write(6,*)' you are trying to turn an integral code on?? '
          write(6,*)' sp  ', def_cando_sp
          write(6,*)' nw  ', def_cando_nw
          write(6,*)' txs ', def_cando_txs
          call util_flush(6)
        endif
        call errquit
     &      ('int_init: logic error with user cando settings',911)
      endif
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
      call defNxyz(ang2use)
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exact_mem(bases,nbas)
      call sp_init(nbas,bases)
      call init70               ! To generate tables etc.
      call int_acc_std()
      call texas_init(rtdb,nbas,bases,txs_mem_min)
      init_int = 1
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


