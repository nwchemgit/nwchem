c $Id: intd_init.F,v 1.22 1999-10-26 22:27:29 mg201 Exp $
*
c:tex-%API Initialization and Termination Routines
c:tex-\subsection{intd\_init}
c:tex-This is the main initialization routine for integral derivatives.
c:tex-Default memory requirements, accuracy thresholds, and other 
c:tex-initializations for all base integral codes are set here. 
c:tex-This routine will read (from the rtdb) any integral 
c:tex-settings changed by the user.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intd_init(rtdb,nbas,bases)
c:tex-\end{verbatim}
c
c Initializes integral derivative code 
c
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "apiP.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
#include "util.fh"
c::functions
      logical  spcart_init
      external spcart_init
      logical  int_ecp_init
      external int_ecp_init
      logical  int_rel_init
      external int_rel_init
      logical  texas_check_basis_ok
      external texas_check_basis_ok
c::passed
c:tex-\begin{verbatim}
      integer rtdb        ! [input] run time data base handle
      integer nbas        ! [input] number of basis sets to be used
      integer bases(nbas) ! [input] basis set handles
c:tex-\end{verbatim}
c::local
      integer ibas, ang2use, angm
      logical status, oprint
      integer intd_memtmp
      integer nqmax_texas  ! maximum number of quartets in texas blocking interface
      parameter (nqmax_texas = 10000)
      integer txs_mem_min
      integer type
      logical cando_txs_deriv
c
c
c     print info/warnings unless print set to none. errors always print.
c
      oprint = util_print('information',print_low)

*      call int_init(nbas,bases)
c
      call int_mem_zero()
c
      DCexp     = 0.0D00
      DCcoeff   = 1.0D00
      val_int_acc = 0.0d00
c
      if(init_intd.eq.1 .and. oprint) then
        write(luout,*)' warning nested intd_inits'
        write(luout,*)' intd_init already called '
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
cTEMPORARY WORKAROUND check for derivative flag intd:cando_txs
c
      cando_txs_deriv = .true.
      if (rtdb_get(rtdb,'intd:cando_txs',MT_LOG,1,status))
     $     cando_txs_deriv = status
      if (.not.cando_txs_deriv) then
         call int_app_set_no_texas(rtdb)
      endif
c
c initialize cando information from rtdb
c
      user_cando_sp   = .false.
      user_cando_nw   = .false.
      user_cando_txs  = .false.
      user_cando_hnd  = .false.
      def_cando_sp    = .false.
      def_cando_nw    = .false.
      def_cando_txs   = .false.
      def_cando_hnd   = .false.
c
      if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,status)) then
        user_cando_sp = .true.
        def_cando_sp  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intd_init: cando_sp set to always be ',def_cando_sp
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
        user_cando_nw = .true.
        def_cando_nw  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intd_init: cando_nw set to always be ',def_cando_nw
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
        user_cando_txs = .true.
        def_cando_txs  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intd_init: cando_txs set to always be ',def_cando_txs
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,status)) then
        user_cando_hnd = .true.
        def_cando_hnd  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intd_init: cando_hnd set to always be ',def_cando_hnd
          call util_flush(6)
        endif
      endif
*
      if (.not.user_cando_txs) then
        if (.not.texas_check_basis_ok(nbas,bases)) then
          user_cando_txs = .true.
          def_cando_txs = .false.
          if (ga_nodeid().eq.0 .and. oprint) then
            write(luout,*)
     &          ' intd_init: internal texas instability possible cando',
     &          '_txs set to always be ',def_cando_txs
            call util_flush(6)
          endif
        endif
      endif
* sanity checking: e.g., you only want to turn off a particular integral
* code never always turn it on.
*
      if (def_cando_sp.or.def_cando_nw.or.def_cando_txs.or.
     &    def_cando_hnd) then
        if (ga_nodeid().eq.0) then
          write(luout,*)' you are trying to turn an integral code on? '
          write(luout,*)' sp  ', def_cando_sp
          write(luout,*)' nw  ', def_cando_nw
          write(luout,*)' txs ', def_cando_txs
          write(luout,*)' hnd ', def_cando_hnd
          call util_flush(6)
        endif
        call errquit
     &      ('intd_init: logic error with user cando settings',911)
      endif
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'intd_init')
00100 continue
      if (.not.status) then
        write(luout,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(luout,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('intd_init: basis handles hosed ',nbas)
      endif
*      write(luout,*)' intd_init: basis set handles valid '
c
c check for both sp and gc shells
c
      call int_bothsp_gc_check(bases,nbas,'intd_init')
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('intd_init: angm error',angm)
        ang2use = max(ang2use,angm)
00300 continue
*
* test for higher than h functions  0123456
      if (ang2use.ge.6) call errquit
     &    ('only basis sets with s through h functions are allowed',
     &    911)
*
c.. for first derivatives add 1
      call defNxyz(ang2use+1)
c
c initialize spcart stuff 
c
      if (.not.(spcart_init((ang2use+1),.true.,.false.))) then
        call errquit('intd_init: spcart_init failed',911)
      endif
c
c Initialize relativistic basis set stuff
c
      any_rel = bas_rel_exists(rtdb)
      if (any_rel) then
        if (.not. int_rel_init(rtdb,nbas,bases)) call errquit
     &      ('int_init: int_rel_init failed ',911)
      end if
c
c.. read in approximate memory from rtdb
c.  
c.. parameter is default value used in hf2d
c.. input memory in words should be scaled by 1/12 for same range
      intd_memthresh = intd_memp
      if (rtdb_get(rtdb,'intd:approxmem',MT_INT,1,intd_memtmp)) then
        if(ga_nodeid().eq.0)then
          write(luout,'(/a,i10,a/)')
     &        ' approximate memory for derivative integrals set to:',
     &        intd_memtmp
        endif
        intd_memthresh = intd_memtmp/12
      endif
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exactd_mem(rtdb,bases,nbas)
      call sp_der_init(nbas,bases)
      call int_acc_std()
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
         call texas_init(rtdb,nbas,bases,nqmax_texas,txs_mem_min,
     &        'der1_int')
      endif
      call hnd_initd()
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
        if (.not.ecp_check_handle(ecp_bsh,'intd_init')) call errquit
     &        ('intd_init: ecp handle is invalid fatal error',911)
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
        if (.not.so_check_handle(so_bsh,'intd_init')) call errquit
     &        ('intd_init: so handle is invalid fatal error',911)
      endif
      if (any_so.or.any_ecp) then
        if (.not.int_ecp_init(ecp_bsh,so_bsh,1)) call errquit
     &        ('intd_init: int_ecp_init failed ',911)
      endif
      init_int = 1
      init_intd = 1
      end
      subroutine sp_der_init(nbas,bases)
      implicit none
#include "stdio.fh"
#include "apispP.fh"
#include "bas.fh"
#include "mafdecls.fh"
c
c:: functions
      logical api_is_ecp_basis
      external api_is_ecp_basis
c::passed
      integer nbas
      integer bases(nbas)
c::local
      integer ncont,icont,ibas
      integer mytype,nprim,ngen,sphorcart
c
      no_sp_conts = .true.
      do 00100 ibas = 1,nbas
        if (.not.api_is_ecp_basis(bases(ibas))) then
          if (.not.bas_numcont(bases(ibas),ncont))
     &        call errquit('sp_der_init: bas_numcont failed',911)
          do 00200 icont = 1,ncont
            if(.not.bas_continfo
     &          (bases(ibas),icont,mytype,nprim,ngen,sphorcart))
     &          call errquit('sp_der_init: bas_continfo failed',911)
            
            if (mytype.lt.0) then
              no_sp_conts = .false.
              goto 00201
            endif
00200     continue
        endif
00100 continue
00201 continue
c
      if (no_sp_conts) return
c
      write(luout,*)' sp functions cannot be used in the integral ',
     &       'derivative code yet'
      call errquit('sp_der_init',911)
      end

      
