c $Id: int_init.F,v 1.38 1998-07-29 19:01:19 d3e129 Exp $
*
c:tex-%API Initialization and Termination Routines
c:tex-\subsection{int\_init}
c:tex-This is the main initialization routine for integrals.
c:tex-Default memory requirements, accuracy thresholds, and other 
c:tex-initializations for all base integral codes are set here. 
c:tex-This routine will read (from the rtdb) any integral 
c:tex-settings changed by the user.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_init(rtdb, nbas, bases)
c:tex-\end{verbatim}
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
#include "int_nbf.fh"
c::functions
      logical  spcart_init
      external spcart_init
      logical  int_ecp_init
      external int_ecp_init
      logical  texas_check_basis_ok
      external texas_check_basis_ok
c::passed
c:tex-\begin{verbatim}
      integer rtdb        ! [input] run time data base handle
      integer nbas        ! [input] number of basis sets to be used
      integer bases(nbas) ! [input] basis set handles
c:tex-\end{verbatim}
c::local
      integer txs_mem_min ! memory from texas
      integer ibas, ang2use, angm, type
      logical status
      integer nqmax_texas  ! maximum number of quartets in texas blocking interface
      parameter (nqmax_texas = 10000)
c
c      block data api_data
c
c
c Block data structure to initialize the common block variables in the
c  internal basis set object data structures
c
c     
      call int_mem_zero()
c
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
c initialize type-> nbf maps
c
      int_nbf_x(-1) = 4
      int_nbf_s(-1) = 4
      do type = 0,int_nbf_max_ang
        int_nbf_x(type) = (type+1)*(type+2)/2
        int_nbf_s(type) = 2*type+1
      enddo
     
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
*
      if (.not.user_cando_txs) then
        if (.not.texas_check_basis_ok(nbas,bases)) then
          user_cando_txs = .true.
          def_cando_txs = .false.
          if (ga_nodeid().eq.0) then
            write(6,*)
     &          ' int_init: internal texas instability possible cando_',
     &          'txs set to always be ',def_cando_txs
            call util_flush(6)
          endif
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
c check for both sp and gc shells
c
      call int_bothsp_gc_check(bases,nbas,'int_init')
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('int_init: angm error',angm)
        ang2use = max(ang2use,angm)
00300 continue
*
* test for higher than h functions  0123456
      if (ang2use.ge.6) call errquit
     &    ('only basis sets with s through h functions are allowed',
     &    911)
*
      call defNxyz(ang2use)
c
c initialize spcart stuff 
c
      if (.not.(spcart_init(ang2use,.true.,.false.))) then
        call errquit('int_init: spcart_init failed',911)
      endif
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exact_mem(rtdb,bases,nbas)
      call sp_init(nbas,bases)
      call init70               ! To generate tables etc.
      call int_acc_std()
* def u=f d=f -> f.and.!f -> f -> e = t
* no txs u=t d=f -> t.and.!f -> t -> e = f
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_init(rtdb,nbas,bases,nqmax_texas,txs_mem_min,
     *                  'scfd_int')
      endif
c
c See if any basis has an attached ECP
c
      any_ecp = .false.
      ecp_bsh = 0
      do ibas = 1,nbas
        if (bas_get_ecp_handle(bases(ibas),ecp_bsh)) then
          any_ecp = .true.
          goto 00001
        endif
      enddo
00001 continue
      if (any_ecp) then
        if (.not.ecp_check_handle(ecp_bsh,'int_init')) call errquit
     &        ('int_init: ecp handle is invalid fatal error',911)
      endif
*
c See if any basis has an attached SO potential
      any_so = .false.
      so_bsh = 0
      do ibas = 1,nbas
        if (bas_get_so_handle(bases(ibas),so_bsh)) then
          any_so = .true.
          goto 00002
        endif
      enddo
00002 continue
      if (any_so) then
        if (.not.so_check_handle(so_bsh,'int_init')) call errquit
     &        ('int_init: so handle is invalid fatal error',911)
      endif
      if (any_so.or.any_ecp) then
        if (.not.int_ecp_init(ecp_bsh,so_bsh,0)) call errquit
     &        ('int_init: int_ecp_init failed ',911)
      endif
      init_int = 1
      end
      subroutine sp_init(nbas,bases)
      implicit none
*_not_used?:#include "apispP.fh"
*_not_used?:#include "bas.fh"
*_not_used?:#include "mafdecls.fh"
*_not_used?:c::passed
      integer nbas
      integer bases(nbas)
*_not_used?:c::local
*_not_used?:      integer ncont,icont,ibas
*_not_used?:      integer mytype,nprim,ngen,sphorcart
*_not_used?:      integer myang, testang
*_not_used?:c
*_not_used?:      no_sp_conts = .true.
*_not_used?:      do 00100 ibas = 1,nbas
*_not_used?:        if (.not.bas_numcont(bases(ibas),ncont))
*_not_used?:     &         call errquit('sp_init: bas_numcont failed',911)
*_not_used?:        do 00200 icont = 1,ncont
*_not_used?:          if(.not.bas_continfo
*_not_used?:     &           (bases(ibas),icont,mytype,nprim,ngen,sphorcart))
*_not_used?:     &           call errquit('sp_init: bas_continfo failed',911)
*_not_used?:
*_not_used?:          if (mytype.lt.0) then
*_not_used?:            no_sp_conts = .false.
*_not_used?:            goto 00201
*_not_used?:          endif
*_not_used?:00200   continue
*_not_used?:00100 continue
*_not_used?:00201 continue
*_not_used?:c
*_not_used?:      if (no_sp_conts) return
*_not_used?:c... arbitrary set max number of blocks to 2 for sp's if spd's set to 3
*_not_used?:      max_sp_blksz = 2
*_not_used?:c
*_not_used?:      write(6,*)':::: allocating sp blocks of memory'
*_not_used?:      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block A',
*_not_used?:     &       hszspa,kszspa))
*_not_used?:     &       call errquit('sp_init:ma_alloc_get A failed',911)
*_not_used?:      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block B',
*_not_used?:     &       hszspb,kszspb))
*_not_used?:     &       call errquit('sp_init:ma_alloc_get B failed',911)
*_not_used?:      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block C',
*_not_used?:     &       hszspc,kszspc))
*_not_used?:     &       call errquit('sp_init:ma_alloc_get C failed',911)
*_not_used?:      if (.not.MA_alloc_get(mt_int,max_sp_blksz,'sp size block D',
*_not_used?:     &       hszspd,kszspd))
*_not_used?:     &       call errquit('sp_init:ma_alloc_get D failed',911)
*_not_used?:      write(6,*)':::done: allocating sp blocks of memory'
*_not_used?:c
*_not_used?:      myang = -56565
*_not_used?:      do 00300 ibas = 1,nbas
*_not_used?:        if (.not.bas_high_angular(bases(ibas),testang))
*_not_used?:     &         call errquit('sp_init:bas_high_angular failed',911)
*_not_used?:        myang = max(myang,testang)
*_not_used?:00300 continue
*_not_used?:      myang = (myang+1)*(myang+2)/2
*_not_used?:* max of (sp,sp|sp,sp) and (sp,X,X,X) and permutations
*_not_used?:      sp_scat_sz = max((4*4*4*4),(4*4*4*myang))
*_not_used?:      sp_scat_sz = max(sp_scat_sz,(4*4*myang*myang))
*_not_used?:      sp_scat_sz = max(sp_scat_sz,(4*myang*myang*myang))
*_not_used?:      sp_scat_sz = max(sp_scat_sz,(myang*myang*myang*myang))
*_not_used?:c
*_not_used?:      if (.not.MA_alloc_get(mt_int,sp_scat_sz,'sp scatter array',
*_not_used?:     &       hspscat,kspscat))
*_not_used?:     &       call errquit('sp_init:ma_alloc_get scat failed',911)
*_not_used?:      write(6,*)':::done: allocating sp blocks of memory'
      end
      logical function texas_check_basis_ok(nbas,bases)
      implicit none
#include "bas.fh"
c::passed
      integer nbas             ! [input] number of basis sets
      integer bases(nbas)      ! [input] basis set handles
*return:: texas_check_basis_ok ! [output] true texas is okay to use
c::local
      integer ibas       ! specific basis index
      integer ishell     ! specific shell index
      integer my_highang ! high angular momentum specified in basis set
      integer my_nshell  ! number of contracted shells in basis set
      integer my_type    ! type of a given contraction
      integer my_nprim   ! number of primitives in a given contraction
      integer my_ngen    ! number of contractions in a given contraction
      integer my_spc            ! spherical/cartesean fag for a given contraction
      double precision exponent ! largest exponent for basis set
c
      texas_check_basis_ok = .true.  
c
      do ibas = 1,nbas
        if (.not.bas_high_angular(bases(ibas),my_highang))
     &      call errquit
     &      ('texas_check_basis_ok: bas_high_angular failed',911)
        if (my_highang.ge.2) then
          if (.not.bas_numcont(bases(ibas),my_nshell))
     &        call errquit
     &        ('texas_check_basis_ok: bas_numcont failed',911)
          do ishell = 1, my_nshell
            if (.not.bas_continfo
     &          (bases(ibas),ishell,my_type,my_nprim,my_ngen,my_spc))
     &          call errquit
     &          ('texas_check_basis_ok: bas_continfo failed',911)
            if (my_type.ge.2) then
              if (.not.
     &            bas_cont_large_exponent(bases(ibas),ishell,exponent))
     &            call errquit
     &            ('texas_check_basis_ok: bas_cont_lrg_expnnt failed',
     &            911)
              if (exponent.gt.10000.0d00) then
                texas_check_basis_ok = .false.
                return
              endif
            endif
          enddo
        endif
      enddo
      end
      
