c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Initialize the integral derivative drivers
C>
C> This is the main initialization routine for integral derivatives.
C> Default memory requirements, accuracy thresholds, and other 
C> initializations for all base integral codes are set here. 
C> This routine will read (from the rtdb) any integral 
C> settings changed by the user.
C>
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
#include "rel_consts.fh"
c
c This routine is a wrapper around the original intd_init which is now
c called intd_init_org. The reason for this extension is related to the
c basis set use in Douglas-Kroll. Douglas-Kroll uses a fitting basis
c set, which could be larger (# of exponents but also higher angular
c momenta), which also needs to be included in the initialization
c procedure. The wrapper sets the Douglas-Kroll stuff and adds the new
c fitting basis set to the list of active basis sets and calls
c intd_init_org.
c
c::passed
c:tex-\begin{verbatim}
      integer rtdb        !< [Input] run time data base handle
      integer nbas        !< [Input] number of basis sets to be used
      integer bases(nbas) !< [Input] basis set handles
c:tex-\end{verbatim}
c::local
      integer mybm_a, mybIm_a, ibas, nrbas     ! Temporary variables
      character*255 dk_bas_name, dk_trans_name ! strings used for Douglas-Kroll
c
c   Dyall-modified Dirac equation
c
      if (.not. rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir))
     &    dyall_mod_dir = .false.
      if (.not. rtdb_get(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
     &    atomic_nesc = .false.
      if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
     &    nesc_1e_approx = .false.
      if (.not. rtdb_get(rtdb,'ss_one_cent',mt_log,1,ss_one_cent))
     &    ss_one_cent = .false.
      if (.not. rtdb_get(rtdb,'do_ssss',mt_log,1,do_ssss)) 
     &    do_ssss = .false.
c
c     Setup array with basis sets from bases
c
      if (.not. ma_push_get(MT_INT,nbas+1,'mybases',mybm_a,mybIm_a))
     $     call errquit('int_init: ma_push_get failed mybases',0,
     &           INT_ERR)
      do 00500 ibas = 1, nbas
         int_mb(mybIm_a+ibas-1) = bases(ibas)
00500 continue
      nrbas = nbas
c
c Initialize the Douglas-Kroll variables
c
c form_type: 1. Free-particle projection operators
c            2. External potential projection operators without pxVp integrals
c            3. Complete spin-free Douglas-Kroll (i.e. including pxVp)
c
      if (.not. rtdb_get(rtdb,'doug_kroll',MT_LOG,1,doug_kroll))
     $    doug_kroll = .false.
      if (doug_kroll) then
c
c        Which form of Douglas-Kroll has to be done, default = 2
c
         if (.not. rtdb_get(rtdb,'doug_kroll:type',MT_INT,1,
     $      dk_form_type)) dk_form_type = 2
c
c        If one want one can specify an even tempered basis as fitting basis
c
         if (.not. rtdb_get(rtdb,'doug_kroll:eventemp',MT_LOG,1,
     $      dk_even_bas)) dk_even_bas = .false.
         if (.not. rtdb_get(rtdb,'doug_kroll:etratio',MT_DBL,1,
     $      dk_ratio)) dk_ratio = 2.0d0
c
c        Generate the fitting basis set and add it to the list of
c        used basis sets
c
         do 00400 ibas = 1, nbas
            if (.not. bas_name(bases(ibas),dk_bas_name,dk_trans_name))
     $         call errquit('int_init: failed in bas_name',0, INT_ERR)
            if (dk_bas_name(1:8) .eq. 'ao basis') goto 00004
00400    continue
         call errquit('int_init: could not find ao basis',0, BASIS_ERR)
00004    call bas_create_uncont(rtdb,bases(ibas),'auto D-K basis',
     $      dk_fitbas)
c
c        Add fitting basis set to list of loaded basis sets for int_init
c
         int_mb(mybIm_a+nbas) = dk_fitbas
         nrbas = nrbas + 1
      endif
      if (doug_kroll .or. dyall_mod_dir) call rel_init_con (rtdb)
      call intd_init_org(rtdb, nrbas, int_mb(mybIm_a))
      if (.not. ma_pop_stack(mybm_a))
     $     call errquit('int_init: ma_pop_stack mybm_a',0, MA_ERR)
c
      end
c
      subroutine intd_init_org(rtdb, nbas, bases)
c
c Initializes integral derivative code 
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
      logical  int_rel_init
      external int_rel_init
      logical  texas_check_basis_ok
      external texas_check_basis_ok
c::passed
      integer rtdb        ! [input] run time data base handle
      integer nbas        ! [input] number of basis sets to be used
      integer bases(nbas) ! [input] basis set handles
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
      user_cando_sim  = .false.
      def_cando_sp    = .false.
      def_cando_nw    = .false.
      def_cando_txs   = .false.
      def_cando_hnd   = .false.
      def_cando_sim   = .false.
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
c
      if (rtdb_get(rtdb,'int:cando_sim',MT_LOG,1,status)) then
        user_cando_sim = .true.
        def_cando_sim  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' intd_init: cando_sim set to always be ',def_cando_sim
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
          write(luout,*)' sim ', def_cando_sim
          call util_flush(6)
        endif
        call errquit
     &      ('intd_init: logic error with user cando settings',911, 
     &         INT_ERR)
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
        call errquit('intd_init: basis handles hosed ',nbas, BASIS_ERR)
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
     &         call errquit('intd_init: angm error',angm, INT_ERR)
        ang2use = max(ang2use,angm)
        if (dyall_mod_dir) then
          if(.not.bas_rel_high_ang(bases(ibas),angm))
     &        call errquit('int_init: angm error',angm, INT_ERR)
          ang2use = max(ang2use,angm+1)
        end if
00300 continue
*
* test for higher than h functions  0123456
      if (ang2use.ge.8) call errquit
     &    ('only basis sets with s through h functions are allowed',
     &    911 , INT_ERR)
*
c.. for first derivatives add 1
      call defNxyz(ang2use+1)
c
c initialize spcart stuff 
c
      if (.not.(spcart_init((ang2use+1),.true.,.false.))) then
        call errquit('intd_init: spcart_init failed',911, INT_ERR)
      endif
c
c Initialize DmD relativistic basis set stuff
c
      if (dyall_mod_dir) then
        if (.not. int_rel_init(rtdb,nbas,bases)) call errquit
     &      ('int_init: int_rel_init failed ',911, INT_ERR)
        if (.not. nesc_1e_approx) then
          user_cando_sp   = .true.
          def_cando_sp   = .false.
          user_cando_txs = .true.
          def_cando_txs  = .false.
          user_cando_sim = .true.
          def_cando_sim  = .false.
          if (ga_nodeid().eq.0 .and. oprint) then
            write(luout,*) ' intd_init: cando_sp,cando_txs,cando_sim',
     &          ' turned off for NESC2E'
            call util_flush(luout)
          end if 
        end if
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
      elseif (.not.(user_cando_nw.and.
     S       (.not.def_cando_nw))) then
c     empty if clause to avoid the call to nwcsim_init
      elseif (.not.(user_cando_sim.and.(.not.def_cando_sim))) then
         call nwcsim_init(rtdb,nbas,bases,1)
      endif
      call hnd_initd(bases,nbas)
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
     &        ('intd_init: ecp handle is invalid fatal error',911,
     &          INT_ERR)
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
     &        ('intd_init: so handle is invalid fatal error',911,
     &          INT_ERR)
      endif
      if (any_so.or.any_ecp) then
        if (.not.int_ecp_init(ecp_bsh,so_bsh,1)) call errquit
     &        ('intd_init: int_ecp_init failed ',911, INT_ERR)
      endif
      init_int = 1
      init_intd = 1
      end
      subroutine sp_der_init(nbas,bases)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
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
     &        call errquit('sp_der_init: bas_numcont failed',911,
     &              INT_ERR)
          do 00200 icont = 1,ncont
            if(.not.bas_continfo
     &          (bases(ibas),icont,mytype,nprim,ngen,sphorcart))
     &          call errquit('sp_der_init: bas_continfo failed',911,
     &              INT_ERR)
            
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
      call errquit('sp_der_init',911, INT_ERR)
      end

C> @}
