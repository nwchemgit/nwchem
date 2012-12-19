c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Initialize the integral 2nd derivatives drivers
C>
C> This is the main initialization routine for integral second 
C> derivatives.
C> Default memory requirements, accuracy thresholds, and other 
C> initializations for all base integral codes are set here. 
C> This routine will read (from the rtdb) any integral 
C> settings changed by the user.
C>
c:tex-%API Initialization and Termination Routines
c:tex-\subsection{intdd\_init}
c:tex-This is the main initialization routine for integral second 
c:tex-derivatives.
c:tex-Default memory requirements, accuracy thresholds, and other 
c:tex-initializations for all base integral codes are set here. 
c:tex-This routine will read (from the rtdb) any integral 
c:tex-settings changed by the user.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intdd_init(rtdb,nbas,bases)
c:tex-\end{verbatim}
c
c Initializes integral second derivative code 
c
      implicit none
#include "stdio.fh"
#include "errquit.fh"
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
      logical  texas_check_basis_ok
      external texas_check_basis_ok
c::passed
c:tex-\begin{verbatim}
      integer rtdb        !< [Input] run time data base handle
      integer nbas        !< [Input] number of basis sets to be used
      integer bases(nbas) !< [Input] basis set handles
c:tex-\end{verbatim}
c::local
      integer ibas, ang2use, angm
      logical status, oprint
      integer intd_memtmp
      integer nqmax_texas  ! maximum number of quartets in texas blocking interface
      parameter (nqmax_texas = 10000)
      integer txs_mem_min
      integer type
      integer nbf2use, nbf2use_test, maxgdd
      logical cando_txs_deriv
c
c
c     print info/warnings unless print set to none. errors always print.
c
      oprint = util_print('information',print_low)
c
      call int_mem_zero()
c
      DCexp     = 0.0D00
      DCcoeff   = 1.0D00
      val_int_acc = 0.0d00
c
      if(init_intdd.eq.1 .and. oprint) then
        write(luout,*)' warning nested intdd_inits'
        write(luout,*)' intdd_init already called '
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
cTEMPORARY WORKAROUND check for derivative flag intdd:cando_txs
c
      cando_txs_deriv = .true.
      if (rtdb_get(rtdb,'intdd:cando_txs',MT_LOG,1,status))
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
     &        ' intdd_init: cando_sp set to always be ',def_cando_sp
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
        user_cando_nw = .true.
        def_cando_nw  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intdd_init: cando_nw set to always be ',def_cando_nw
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
        user_cando_txs = .true.
        def_cando_txs  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intdd_init: cando_txs set to always be ',def_cando_txs
          call util_flush(6)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,status)) then
        user_cando_hnd = .true.
        def_cando_hnd  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intdd_init: cando_hnd set to always be ',def_cando_hnd
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
     &          ' intdd_init: internal texas instability ',
     &          'possible cando',
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
     &      ('intdd_init: logic error with user cando settings',911, 
     &        INT_ERR)
      endif
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'intdd_init')
00100 continue
      if (.not.status) then
        write(luout,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(luout,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('intdd_init: basis handles hosed ',nbas, INT_ERR)
      endif
*      write(luout,*)' intdd_init: basis set handles valid '
c
c check for both sp and gc shells
c
      call int_bothsp_gc_check(bases,nbas,'intdd_init')
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('intdd_init: angm error',angm, INT_ERR)
        ang2use = max(ang2use,angm)
00300 continue
*
* test for higher than g functions  0123456
      if (ang2use.ge.7) call errquit
     &    ('only basis sets with s through g functions are allowed',
     &    911, INT_ERR)
*
c.. for second derivatives add 2
      call defNxyz(ang2use+2)
c
c initialize spcart stuff 
c
      if (.not.(spcart_init((ang2use+1),.true.,.false.))) then
        call errquit('intdd_init: spcart_init failed',911, INT_ERR)
      endif
c.. read in approximate memory from rtdb
c.  
c.. parameter is default value used in hf2d
c.. input memory in words should be scaled by 1/12 for same range
      intd_memthresh = intd_memp
      if (rtdb_get(rtdb,'intdd:approxmem',MT_INT,1,intd_memtmp)) then
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
     &        'der2_int')
      endif
      call hnd_initdd(bases,nbas)
c
c    confirm basis set maximum buffer size for second derivatives
c
      nbf2use = 0
      do ibas = 1,nbas
        if (.not.bas_nbf_cn_max(bases(ibas),nbf2use_test))
     &      call errquit('intdd_init: bas_nbf_cn_max failed',911,
     &             INT_ERR)
        nbf2use = max(nbf2use,nbf2use_test)
      enddo
      maxgdd  = 78*nbf2use*nbf2use*nbf2use*nbf2use  ! maxg for 2e2nd derivs
      nbf2use = maxgdd + maxgdd/10  ! add 10%
      isz_2e4c = max(isz_2e4c,nbf2use)
      mem_2e4c = max(mem_2e4c,isz_2e4c)
      isz_2e3c = maxgdd/nbf2use
      isz_2e2c = maxgdd/nbf2use/nbf2use
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
        if (.not.ecp_check_handle(ecp_bsh,'intdd_init')) call errquit
     &        ('intdd_init: ecp handle is invalid fatal error',911,
     &         INT_ERR)
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
        if (.not.so_check_handle(so_bsh,'intdd_init')) call errquit
     &        ('intdd_init: so handle is invalid fatal error',911,
     &         INT_ERR)
      endif
      if (any_so.or.any_ecp) then
        if (.not.int_ecp_init(ecp_bsh,so_bsh,2)) call errquit
     &        ('intdd_init: int_ecp_init failed ',911, INT_ERR)
      endif
      init_int   = 1
      init_intd  = 1
      init_intdd = 1
      end

C> @}
