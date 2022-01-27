c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Main initialization routine for integrals
C>
C> This routine is the main initialization routine for the 
C> integral packages. It sets memory limits, accuracy thresholds,
C> and other initializations for all base integral codes. The
C> routine reads any integral settings changed by the user from
C> the runtime database.
C>
C> This routine is a wrapper around the original `int_init` which is now
C> called `int_init_org`. The reason for this extension is related to the
C> basis set use in Douglas-Kroll. Douglas-Kroll uses a fitting basis
C> set, which could be larger (# of exponents but also higher angular 
C> momenta), which also needs to be included in the initialization
C> procedure. The wrapper sets the Douglas-Kroll stuff and adds the new
C> fitting basis set to the list of active basis sets and calls 
C> `int_init_org`.
C>
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
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "apiP.fh"
#include "candoP.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
#include "util.fh"
#include "rel_nwc.fh"
#include "stdio.fh"
c
c This routine is a wrapper around the original int_init which is now
c called int_init_org. The reason for this extension is related to the
c basis set use in Douglas-Kroll. Douglas-Kroll uses a fitting basis
c set, which could be larger (# of exponents but also higher angular 
c momenta), which also needs to be included in the initialization
c procedure. The wrapper sets the Douglas-Kroll stuff and adds the new
c fitting basis set to the list of active basis sets and calls 
c int_init_org.
c
c ... jochen 01/22: the exact decoupling functionality needs 
c     the same setup as Douglas-Kroll, so we added that here
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
      logical ouseracc
c
      integer irtdb
      common /int_rtdb/irtdb   ! rtdb member variable handle
c
c     Set member irtdb
c
      irtdb = rtdb  
c
c     Put integral accuracy under user control
c
      ouseracc = .false.
      if (rtdb_get(rtdb,'int:acc_std',mt_dbl,1,val_int_acc_std)) then
         ouseracc = .true.
      else
c
cniri Reinstate the tighter tolerance. We encountered a subtle problem with a large system where 
cniri a rather diffuse augmented basis set was used. There were a number of dependencies.
cniri We managed to address the issue by ramping up this tolerance.
cniri val_int_acc_std = 1d-12
      val_int_acc_std = 1d-15


      end if
      if (rtdb_get(rtdb,'int:acc_high',mt_dbl,1,val_int_acc_high)) then 
         ouseracc = .true.
      else
         val_int_acc_high = val_int_acc_std**2
      end if
      if (ouseracc .and. ga_nodeid().eq.0 .and. 
     $     util_print('information',print_default)) then
         write(LuOut,606) val_int_acc_std, val_int_acc_high
 606     format(/'int_init: setting std/high accuracies to',1p,2d9.1)
         call util_flush(LuOut)
      end if
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
     &       INT_ERR)
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
         call errquit('int_init: could not find ao basis',0, INT_ERR)
00004    call bas_create_uncont(rtdb,bases(ibas),'auto D-K basis',
     $      dk_fitbas)
c
c        Add fitting basis set to list of loaded basis sets for int_init
c
         int_mb(mybIm_a+nbas) = dk_fitbas
         nrbas = nrbas + 1
      endif

c ... jochen 01/22: same idea as for douglas-kroll: we initialize
c     an additional uncontracted basis. 
c     in principle we could have used the corresponding D-K code
c     and variables, but I want to avoid potentially conflicting
c     variable definitions

      if (.not. rtdb_get(rtdb,'x2c',MT_LOG,1,exact_2c)) 
     &    exact_2c = .false.
c      write (luout,*) 'jochen: in int_init x2c = ',exact_2c
      if (exact_2c) then

        if (doug_kroll) call errquit(
     &      'int_init: x2c and DK not allowed',1,INPUT_ERR)
        
c       Generate the fitting basis set and add it to the list of 
c       used basis sets
        
        do 300 ibas = 1, nbas
          if (.not. bas_name(bases(ibas),dk_bas_name,dk_trans_name))
     &        call errquit('int_init: x2c failed in bas_name',0,
     &        INT_ERR)
          if (dk_bas_name(1:8) .eq. 'ao basis') goto 350
 300    continue
        call errquit('int_init: x2c could not find ao basis',0,
     &      INT_ERR)
 350    call bas_create_uncont(rtdb,bases(ibas),'auto X2C basis',
     &      x2c_bas)
       
c       Add fitting basis set to list of loaded basis sets
       
        int_mb(mybIm_a+nbas) = x2c_bas
        nrbas = nrbas + 1
      endif ! exact_2c
c ... jochen end

c     rel_init_con initializes some relativistic constants
      if (doug_kroll .or. dyall_mod_dir .or. exact_2c)
     &    call rel_init_con (rtdb)

      call int_init_org(rtdb, nrbas, int_mb(mybIm_a))
      if (.not. ma_pop_stack(mybm_a))
     $     call errquit('int_init: ma_pop_stack mybm_a',0, INT_ERR)
c
      end
c
C>
C> \brief Initialization routine for integrals
C>
C> This routine is the initialization routine for the 
C> integral packages. It sets memory limits, accuracy thresholds,
C> and other initializations for all base integral codes. The
C> routine reads any integral settings changed by the user from
C> the runtime database. Only data associated with relativistic
C> basis sets are handled in `int_init`.
C>
      subroutine int_init_org(rtdb, nbas, bases)
c
c initializes integral code to data structers for a integral computation
c
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "apiP.fh"
#include "candoP.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
#include "util.fh"
#include "stdio.fh"
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
      integer rtdb        !< [Input] run time data base handle
      integer nbas        !< [Input] number of basis sets to be used
      integer bases(nbas) !< [Input] basis set handles
c::local
      integer txs_mem_min ! memory from texas
      integer ibas, ang2use, angm, type
      logical status, oprint
      integer nqmax_texas  ! maximum number of quartets in texas blocking interface
      parameter (nqmax_texas = 10000)

c
c     print info/warnings unless print set to none. errors always print.
c
      oprint = util_print('information',print_low)
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
      if(init_int.eq.1 .and. oprint) then
        write(luout,*)' warning nested int_inits'
        write(luout,*)' int_init already called '
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
     &        ' int_init: cando_sp set to always be ',def_cando_sp
          call util_flush(luout)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
        user_cando_nw = .true.
        def_cando_nw  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' int_init: cando_nw set to always be ',def_cando_nw
          call util_flush(luout)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
        user_cando_txs = .true.
        def_cando_txs  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' int_init: cando_txs set to always be ',def_cando_txs
          call util_flush(luout)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,status)) then
        user_cando_hnd = .true.
        def_cando_hnd  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' int_init: cando_hnd set to always be ',def_cando_hnd
          call util_flush(luout)
        endif
      endif
c
      if (rtdb_get(rtdb,'int:cando_sim',MT_LOG,1,status)) then
        user_cando_sim = .true.
        def_cando_sim  = status
        if (ga_nodeid().eq.0 .and. oprint) then
          write(luout,*)
     &        ' int_init: cando_sim set to always be ',
     S          def_cando_sim
          call util_flush(luout)
        endif
      endif
c
      if (.not.user_cando_txs) then
        if (.not.texas_check_basis_ok(nbas,bases)) then
          user_cando_txs = .true.
          def_cando_txs = .false.
          if (ga_nodeid().eq.0 .and. oprint) then
            write(luout,*)
     &          ' int_init: internal texas instability possible cando_',
     &          'txs set to always be ',def_cando_txs
            call util_flush(luout)
          endif
        endif
      endif
* sanity checking: e.g., you only want to turn off a particular integral
* code never always turn it on.
*
      if (def_cando_sp.or.def_cando_nw.or.def_cando_txs.or.
     &    def_cando_hnd.or.def_cando_sim) then
        if (ga_nodeid().eq.0) then
          write(luout,*)' you are trying to turn an integral code on? '
          write(luout,*)' sp  ', def_cando_sp
          write(luout,*)' nw  ', def_cando_nw
          write(luout,*)' txs ', def_cando_txs
          write(luout,*)' hnd ', def_cando_hnd
          write(luout,*)' sim ', def_cando_sim
          call util_flush(luout)
        endif
        call errquit
     &      ('int_init: logic error with user cando settings',911,
     &       INT_ERR)
      endif
c
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'int_init')
00100 continue

      if (.not.status) then
        write(luout,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(luout,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('int_init: basis handles hosed ',nbas, INT_ERR)
      endif
*      write(luout,*)' int_init: basis set handles valid '
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
     &      call errquit('int_init: angm error',angm, INT_ERR)
        ang2use = max(ang2use,angm)
        if (dyall_mod_dir) then
          if(.not.bas_rel_high_ang(bases(ibas),angm))
     &        call errquit('int_init: angm error',angm, INT_ERR)
          ang2use = max(ang2use,angm+1)
        end if
00300 continue
*
* test for higher than L functions  0123456789
      if (ang2use.ge.9) call errquit
     &    ('only basis sets with s through l functions are allowed',
     &    911, BASIS_ERR)
*
      call defNxyz(ang2use)
c
c initialize spcart stuff 
c
      if (.not.(spcart_init(ang2use,.true.,.false.))) then
        call errquit('int_init: spcart_init failed',911, INT_ERR)
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
          if (ga_nodeid().eq.0 .and. oprint) then
            write(luout,*) ' int_init: cando_sp and cando_txs',
     &          ' turned off for NESC2E'
            call util_flush(luout)
          end if 
        end if
      end if
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exact_mem(rtdb,bases,nbas)
      call sp_init(nbas,bases)
      call init70               ! To generate tables etc.
      call hnd_init(bases,nbas)
      call int_acc_std()
* def u=f d=f -> f.and.!f -> f -> e = t
* no txs u=t d=f -> t.and.!f -> t -> e = f
      if (.not.(user_cando_txs.and.(.not.def_cando_txs))) then
        call texas_init(rtdb,nbas,bases,nqmax_texas,txs_mem_min,
     *                  'scfd_int')
      elseif (.not.(user_cando_nw.and.
     S       (.not.def_cando_nw))) then
         call igamma_init
      else
c     simint as last option (e.g. set int:cando_txs f set int:cando_nw f)
         call nwcsim_init(rtdb,nbas,bases,0)
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
     &        ('int_init: ecp handle is invalid fatal error',911,
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
        if (.not.so_check_handle(so_bsh,'int_init')) call errquit
     &        ('int_init: so handle is invalid fatal error',911, 
     &          INT_ERR)
      endif
      if (any_so.or.any_ecp) then
        if (.not.int_ecp_init(ecp_bsh,so_bsh,0)) call errquit
     &      ('int_init: int_ecp_init failed ',911, INT_ERR)
        if (dyall_mod_dir .or. doug_kroll) call errquit('int_init: '//
     &      'do not mix ECPs and all-electron relativistic methods',911,
     &       INT_ERR)
      endif
      init_int = 1
      end
      subroutine sp_init(nbas,bases)
      implicit none
      integer nbas
      integer bases(nbas)
      end
C>
C> \brief Check if the integrals for all basis sets can be evaluated by
C> the Texas integral package
C>
      logical function texas_check_basis_ok(nbas,bases)
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "util.fh"
c::passed
      integer nbas             !< [Input] number of basis sets
      integer bases(nbas)      !< [Input] basis set handles
*return:: texas_check_basis_ok !< [Output] true texas is okay to use
c::local
      integer ibas       ! specific basis index
      integer ishell     ! specific shell index
      integer my_highang ! high angular momentum specified in basis set
      integer my_nshell  ! number of contracted shells in basis set
      integer my_type    ! type of a given contraction
      integer my_nprim   ! number of primitives in a given contraction
      integer my_ngen    ! number of contractions in a given contraction
      integer my_spc            ! spherical/cartesian flag for a given contraction
      double precision exponent ! largest exponent for basis set
c
      texas_check_basis_ok = .true.  
c
      do ibas = 1,nbas
        if (.not.bas_high_angular(bases(ibas),my_highang))
     &      call errquit
     &      ('texas_check_basis_ok: bas_high_angular failed',911,
     &       BASIS_ERR)
        if (my_highang.ge.2) then
          if (.not.bas_numcont(bases(ibas),my_nshell))
     &        call errquit
     &        ('texas_check_basis_ok: bas_numcont failed',911,
     &       BASIS_ERR)
          do ishell = 1, my_nshell
            if (.not.bas_continfo
     &          (bases(ibas),ishell,my_type,my_nprim,my_ngen,my_spc))
     &          call errquit
     &          ('texas_check_basis_ok: bas_continfo failed',911,
     &       BASIS_ERR)
            if (my_type.ge.2) then
              if (.not.
     &            bas_cont_large_exponent(bases(ibas),ishell,exponent))
     &            call errquit
     &            ('texas_check_basis_ok: bas_cont_lrg_expnnt failed',
     &            911, BASIS_ERR)
              if (exponent.gt.10000.0d00) then
                texas_check_basis_ok = .false.
                return
              endif
            endif
          enddo
        endif
      enddo
      end
c
c     == access function for rtdb member variable handle ==
      integer function int_get_rtdb()
c
      implicit none
c
      integer irtdb
      common /int_rtdb/irtdb

      int_get_rtdb = irtdb
c   
      return
      end
C> @}
