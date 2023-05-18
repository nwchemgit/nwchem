      logical function nadct_trace(rtdb, nadct, trans)
c     
c Compute the terms that go into the nonadiabatic derivative coupling
c terms 
c                     f_I,J(R)^x = <I|d/dR_x J>
c
c Each term is computed by a trace of data from COLUMBUS with 
c gradient integrals from NWChem.
c 
c This is also computing gradients apparently
 
      implicit none

#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "global.fh"
#include "util.fh"
#include "schwarz.fh"
#include "dftps.fh"
#include "zora.fh"
#include "rel_consts.fh"
#include "case.fh"
#include "dft_scaleig.fh"
#ifdef SCALAPACK
      integer ga_cholesky, ga_llt_i
      external ga_cholesky, ga_llt_i
#endif
c
      integer rtdb
      logical nadct
      integer me, nproc, cinfo
      logical IOLGC, oprint, oprint_time, oprint_ioinfo
      logical ldirect
      integer iVxc_opt, lmiss, ivcoul_opt
      integer g_xcinv, g_oi, g_tmpb
      integer lsvals, isvals, n_dep
      double precision toll_s
      integer nat
      integer g_force
      integer lforce
      integer l_force, l_frc_nuc, l_frc_kin, l_frc_fds, l_frc_2el, 
     &        l_frc_sow, l_frc_sox, l_frc_soy, l_frc_soz, l_frc_atr
      integer k_force, k_frc_nuc, k_frc_kin, k_frc_fds, k_frc_2el, 
     &        k_frc_sow, k_frc_sox, k_frc_soy, k_frc_soz, k_frc_atr
      integer aodens, aodensa
      parameter(aodens=24, aodensa=25)
      integer ninfo, ninfmx, ntitmx, nsymx
      parameter(ninfo=5,ninfmx=10,ntitmx=30,nsymx=8)
      integer infocd(ninfmx)
      character*80 siftits(ntitmx)
      integer nsym, nmot, nmbpsy(nsymx), nsift
      integer infocda(ninfmx)
      character*80 siftitsa(ntitmx)
      integer nsifta
      LOGICAL :: trans ! is this a transition property run?

cgk debug
*      write(*,*)'gk: in nadct_trace'
cgk end
      nadct_trace=.true.
c  
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('dft cinfo', print_default)
      oprint_time = util_print('dft timings', print_high)
      oprint_ioinfo = util_print('io info', print_debug)
      if(ga_nodeid().eq.0.and.
     .     util_print('geometry', print_high)) then
        if(.not.geom_print(geom)) call errquit(
     .       ' nadct_trace: geom print failed',0, GEOM_ERR)
      endif
      IOLGC = .TRUE.
      if (noio.eq.1) IOLGC = .FALSE.
!#ifdef SCALAPACK
!#if defined(PARALLEL_DIAG)
!      call slexit
!      call slexit2
!      if (.not. rtdb_cget(rtdb, 'dft:scaladiag',1,yev))
!     &     yev='d'
!      if (.not. rtdb_get(rtdb, 'dft:scalanb',mt_int,1,nb))
!     &     nb=0
!#endif
!#endif
c     
c     Determine the characteristics of the AO Gaussian basis set.
c     This uses the old NWdft data structures, rather than tthe rtdb
c     
      if (.not. bas_numbf(AO_bas_han, nbf_ao)) then
         call errquit('Exiting from nadct_trace',1, BASIS_ERR)
      endif 
cgk debug
*     write(*,*)'in nadct_trac, AO_bas_han = ', AO_bas_han
cgk end
*     iter = 0
      if(nproc.gt.1) then
        write(6,'(a,a)')'WARNING: Writing aoints, calculating 1e ints ',
     &    'with 1 CPU'
        write(6,'(a,a)')'WARNING: we are currently supporting 1 ',
     &    'process for aoints'
      endif
c
c     Check zora options
      do_zora_get_correction=.false.
      if (rtdb_get(rtdb,'zora:calc',mt_log,1,
     &    do_zora_get_correction)) 
     &    do_zora_get_correction=.true.
c
      do_zora = .false.
      if (rtdb_get(rtdb,'zora:on',mt_log,1,do_zora)) 
     &     do_zora = .true.
c
c     == Get the grid cutoff or set default ==
      if(.not.rtdb_get(rtdb,'zora:cutoff',mt_dbl,1,zoracutoff))
     &  zoracutoff = 1.0d-08
c
c     Get the light speed and other constants, if specified
      call rel_init_con(rtdb)
      clight_au = cau
c
*c     Get and check the Coulomb Attenuation Method (CAM) parameters
*      if (.not.rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1,
*     &   cam_exch)) cam_exch=.false.
*      if (.not.rtdb_get(rtdb, 'dft:cam_omega', mt_dbl, 1,
*     &   cam_omega)) cam_omega=0.d0
*      if (.not.rtdb_get(rtdb, 'dft:cam_alpha', mt_dbl, 1,
*     &   cam_alpha)) cam_alpha=0.d0
*      if (.not.rtdb_get(rtdb, 'dft:cam_beta', mt_dbl, 1,
*     &   cam_beta)) cam_beta=0.d0
*c
*      ldirect = .false.
*      if (cam_exch .and. (.not.CDFIT) ) then
*       if (rtdb_get(rtdb, 'dft:direct', mt_log, 1, ldirect)) then
*        if (.not.ldirect) then
*          if (ga_nodeid().eq.0) then
*           write(luout,*) "CAM type calculations without CDFIT are only"
*           write(luout,*) "available with the direct option"
*          end if
*          call errquit('Please use the direct option',4, INPUT_ERR)
*        end if
*       end if
*      end if
c     
c     Select method to compute the XC potential.
      if (XC_bas_han.eq.-99999)then
         iVxc_opt = 0
         XCFIT    = .FALSE.
      else
         call errquit('nadct_trace: XC basis fitting not implemented',
     &     0,CAPMIS_ERR)
      endif

*      if (XCFIT)then
*c     
*c        Determine the characteristics of the XC Gaussian basis set.
*c     
*         if (.not. bas_numbf(XC_bas_han, nbf_xc)) then
*            call errquit('Exiting from scfdriver',4, BASIS_ERR)
*         endif 
*c     
*c           Compute the matrix inverse of the 2-ctr OIs.
*c     
**ga:1:0
*               if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'XCinv', 
*     &                             nbf_xc, 0, g_xcinv)) 
*     &            call errquit('nadct_trace: error creating g_xcinv',
*     &              0,GA_ERR)
*               call ga_zero(g_xcinv)
**ga:1:0
*            if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'g_oi', 
*     &                          nbf_xc, 0, g_oi)) 
*     &         call errquit('nadct_trace: error creating g_oi',0,
*     &       GA_ERR)
*            call ga_zero(g_oi)
*            call int_1e_ga(xc_bas_han, xc_bas_han,
*     .           g_oi, 'overlap', oskel)
*            if (oskel)call sym_symmetrize(
*     (           geom, xc_bas_han, .false., g_oi)
*#if defined(PARALLEL_DIAG)
*            call ga_chol(g_oi, g_xcinv, cinfo)
*#else 
*            call ga_chol_seq(g_oi, g_xcinv, cinfo)
*#endif
*            if (cinfo.ne.0)then
*               if (me.eq.0)then
*                  write(LuOut,*)' Problem in performing a Choleski '
*                  write(LuOut,*)' decomposition of the 2-ctr overlap '
*                  write(LuOut,*)' matrix using XC fitting basis. '
*                  write(LuOut,*)' Attempting a diag/inverse. '
*               endif
*            endif
*            if (cinfo.eq.0) then
*#if defined(PARALLEL_DIAG)
*               call ga_inverse(g_xcinv, g_oi)
*               if (odftps) call pstat_on(ps_bdgemm)
*               call ga_dgemm('T', 'N', nbf_xc, nbf_xc, nbf_xc, 1.d0, 
*     &              g_oi, g_oi, 0.d0, g_xcinv)
*               if (odftps) call pstat_off(ps_bdgemm)
*#else 
*               call ga_copy(g_xcinv, g_oi)
*               call ga_inv_seq(g_oi, g_xcinv)
*#endif
*            else
*               call dft_invdiag(g_oi, g_xcinv, nbf_xc)
*            endif
*            if (.not. ga_destroy(g_oi)) call errquit
*     &         ('nadct_trace: could not destroy g_oi', 0, GA_ERR)
*            if (me.eq.0.and.IOLGC) then
*               lmiss = 0
*               call dft_invio('XC', g_xcinv, nbf_xc, 'WRITE', lmiss)
*               if (lmiss.ne.0)call errquit
*     &         (' nadct_trace; dft_invio - abnormal write of XC ', 0,
*     &       DISK_ERR)
*            endif
*         if (IOLGC)then
*            if (.not. ga_destroy(g_xcinv)) call errquit
*     &         ('nadct_trace: could not destroy g_xcinv', 0, GA_ERR)
*         endif
*      endif
*c
c     Done with XC setup stuff....
c     Select method to compute Coulomb potential.
c     
      if (CD_bas_han.eq.-99999)then
         iVcoul_opt = 0
         CDFIT      = .FALSE.
      else
         call errquit('nadct_trace: charge density fitting not implement
     &ed', 0,CAPMIS_ERR)
      endif
c
c
*      if (CDFIT)then
*c     
*c        Determine the characteristics of the CD Gaussian basis set.
*c     
*         if (odftps) call pstat_on(ps_2cer)
*         if (.not. bas_numbf(CD_bas_han, nbf_cd)) then
*            call errquit('Exiting from scfdriver',5, BASIS_ERR)
*         endif 
*c     
*c        Compute the matrix inverse of the CD 2-ctr ERIs.
*c     
*         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERI', 
*     &                       0, nbf_cd, g_2ceri)) 
*     &      call errquit('nadct_trace: error creating g_2ceri',0,
*     &       GA_ERR)
*         call ga_zero(g_2ceri)
*         call dft_get2eri(CD_bas_han, g_2ceri,oskel)
*         if (odftps) call pstat_off(ps_2cer)
*         if (oskel)call 
*     .        sym_symmetrize(geom,cd_bas_han,.false.,g_2ceri)
*         call ga_sync()
*         if (.not. ga_duplicate(g_2ceri, g_cdinv, 'CD 2cERInv'))
*     &    call errquit('nadct_trace: error creating g_cdinv',0, 
*     &      GA_ERR)
*         if (odftps) call pstat_on(ps_inv)
*         call ga_zero(g_cdinv)
*         lmiss = 1
*c         if (lmiss.eq.1) then
*            call ga_zero(g_cdinv)
*            cinfo = 0
*#if defined(PARALLEL_DIAG)
*#ifdef SCALAPACK
*            call ga_copy(g_2ceri, g_cdinv)
*            call ga_sync()
*            cinfo= ga_cholesky('U',g_cdinv)
*#else
*            call ga_chol(g_2ceri, g_cdinv, cinfo)
*#endif
*#else 
*            call ga_chol_seq(g_2ceri, g_cdinv, cinfo)
*#endif
*            if (cinfo.ne.0)then
*               if (me.eq.0)then
*                  write(LuOut,*)' Problem in performing a Choleski '
*                  write(LuOut,*)' decomposition of the 2-ctr ERI '
*                  write(LuOut,*)' matrix using CD fitting basis. '
*                  write(LuOut,*)' Attempting a diag/inverse. '
*               endif
*            endif
*            if (cinfo.eq.0) then
*               g_tmpb = g_2ceri
*#if defined(PARALLEL_DIAG)
*#ifdef SCALAPACK
*               cinfo = ga_llt_i('U',g_cdinv,-1)
*            if (cinfo.ne.0)then
*               if (me.eq.0)then
*                  write(LuOut,*)' Problem in performing a Invers. '
*                  write(LuOut,*)' of the 2-ctr ERI '
*               endif
*               call ga_sync
*               call errquit('nadct_trace: inverse failed ',0,0)
*            endif
*               
*#else
*               call ga_inverse(g_cdinv, g_tmpb)
*#endif
*               if (odftps) call pstat_on(ps_bdgemm)
*#ifndef SCALAPACK
*               call ga_dgemm('T', 'N', nbf_cd, nbf_cd, nbf_cd, 1.d0, 
*     &              g_tmpb, g_tmpb, 0.d0, g_cdinv)
*#endif
*               if (odftps) call pstat_off(ps_bdgemm)
*#else 
*               call ga_copy(g_cdinv, g_tmpb)
*               call ga_inv_seq(g_tmpb, g_cdinv)
*#endif
*            else
*               call dft_invdiag(g_2ceri, g_cdinv, 
*     &                          nbf_cd)
*            endif
*            if (odftps) call pstat_off(ps_inv)
*#ifndef SCALAPACK            
*c
*c     second build of g_2ceri needed becuase previous calls destroyed it
*c     
*            if (odftps) call pstat_on(ps_2cer)
*            call ga_zero(g_2ceri)
*            call dft_get2eri(CD_bas_han, g_2ceri,oskel)
*            if (oskel)call 
*     .           sym_symmetrize(geom,cd_bas_han,.false.,g_2ceri)
*            if (odftps) call pstat_off(ps_2cer)
*#endif
*            if (IOLGC.and.(me.eq.0)) then
*               lmiss = 0
*               call dft_invio('CDI', g_cdinv, nbf_cd, 'WRITE', lmiss)
*               if (lmiss.ne.0)call errquit
*     &         (' nadct_trace; dft_invio - abnormal write of CDI ', 0
*     &         ,DISK_ERR)
*               lmiss = 0
*               call dft_invio('CD', g_2ceri, nbf_cd, 'WRITE', lmiss)
*               if (lmiss.ne.0)call errquit
*     &         (' nadct_trace; dft_invio - abnormal write of CD ', 0,
*     &       DISK_ERR)
*            endif
*            if (oprint_ioinfo)then
*               call ga_print(g_cdinv)
*               call ga_print(g_2ceri)
*            endif
*c         endif
*         if (IOLGC) then
*            if (.not. ga_destroy(g_cdinv)) call errquit
*     &         ('nadct_trace: could not destroy g_xcinv', 0, GA_ERR)
*            if (.not. ga_destroy(g_2ceri)) call errquit
*     &         ('nadct_trace: could not destroy g_xcinv', 0, GA_ERR)
*         endif
*      endif
*c
c     Done with CD setup stuff
c
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'ovl eig vals', lsvals, 
     &   isvals))
     &   call errquit('nadct_trace: cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      call dfill(nbf_ao, 0.0d0, dbl_mb(isvals), 1)
c
c     Set orbital overlap tolerance
c     
      call util_tolls(rtdb,.false.,toll_s,n_dep,dbl_mb(isvals),nbf_ao)
c
c     open aodens
c
      if(me.eq.0) then
        if(ga_nodeid().eq.0) then
          open(unit=aodens,status='old',form='unformatted',
     &      file='aodens')
          rewind aodens
        endif

cgk debug
*     write(*,*)'gk: calling rdhcid'
cgk end
        call rdhcid(aodens, infocd, nsift, siftits, nsym, nmot, nmbpsy)
cgk debug
!      write(*,*)'gk: back from rdhcid'
!      write(*,*)'gk: nsym = ', nsym
!      WRITE(*,*)"LB, info=",infocd
cgk end
      endif

      if (nmot.ne.nbf_ao) then
         call errquit('nadct_trace: aodens number of orbitals does not m
     &atch NWChem orbital number', nmot, BASIS_ERR)
      endif 

      if (nsym.ne.1) then
         call errquit('nadct_trace: only c1 symmetry implented...invalid
     & aodens, nsym=', nsym, BASIS_ERR)
      endif 

c     open aodensa
c
      IF (trans) THEN
        if(me.eq.0) then
          if(ga_nodeid().eq.0) then
            open(unit=aodensa,status='old',form='unformatted',
     &        file='aodens.asym')
            rewind aodens
          endif

          call rdhcid(aodensa, infocda,nsifta,siftitsa,nsym,nmot,nmbpsy)
        endif

        if (nmot.ne.nbf_ao) then
         call errquit('nadct_trace: aodens number of orbitals does not m
     &atch NWChem orbital number', nmot, BASIS_ERR)
        endif 

        if (nsym.ne.1) then
         call errquit('nadct_trace: only c1 symmetry implented...invalid
     & aodens, nsym=', nsym, BASIS_ERR)
        endif 
      ENDIF
c     set up force arrays
      if (.not. geom_ncent(geom,nat))
     $     call errquit('nadct_trace: could not get natoms',0, GEOM_ERR)

c     global force array
      if (.not. ga_create(mt_dbl, 3, nat, 'forces', 3, 0, g_force))
     $     call errquit('nadct_trace: failed to create force GA',0,
     &       GA_ERR)

      call ga_zero(g_force)
c     local replication for separate pieces
      lforce = nat * 3
      if (.not. ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     $     call errquit('nadct_trace:could not allocate l_force',lforce,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_nuc, 
     $     k_frc_nuc)) call errquit
     &    ('nadct_trace:could not allocate l_frc_nuc',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_kin,
     $     k_frc_kin)) call errquit
     &    ('nadct_trace:could not allocate l_frc_kin',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_fds,
     $     k_frc_fds)) call errquit
     &    ('nadct_trace:could not allocate l_frc_fds',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_2el,
     $     k_frc_2el)) call errquit
     &    ('nadct_trace:could not allocate l_frc_2el',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_sow,
     $     k_frc_sow)) call errquit
     &    ('nadct_trace:could not allocate l_frc_sow',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_sox,
     $     k_frc_sox)) call errquit
     &    ('nadct_trace:could not allocate l_frc_sox',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_soy,
     $     k_frc_soy)) call errquit
     &    ('nadct_trace:could not allocate l_frc_soy',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_soz,
     $     k_frc_soz)) call errquit
     &    ('nadct_trace:could not allocate l_frc_soz',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'forces', l_frc_atr,
     $     k_frc_atr)) call errquit
     &    ('nadct_trace:could not allocate l_frc_atr',lforce, MA_ERR)

      call dfill(lforce, 0.0D0, dbl_mb(k_force), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_nuc), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_kin), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_fds), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_sow), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_sox), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_soy), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_soz), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_atr), 1)
c
c     Calculate non spin-orbit energy gradient terms
c       trace( del(h)*D + del(g)*d - del(S)*F )
c
cgk debug
*     write(*,*)'gk: calling egrad_trace'
*     write(*,*)'gk:  nmot = ', nmot
*     write(*,*)'gk:  nbf_ao = ', nbf_ao
cgk end
      call egrad_trace(AO_bas_han, aodens, ninfo, infocd, nbf_ao,
     &  g_force, dbl_mb(k_frc_nuc), dbl_mb(k_frc_kin), 
     &  dbl_mb(k_frc_fds), dbl_mb(k_frc_2el), dbl_mb(k_frc_sow),
     &  dbl_mb(k_frc_sox), dbl_mb(k_frc_soy), dbl_mb(k_frc_soz),
     &  dbl_mb(k_force), nat, rtdb, nadct, aodensa, infocda, trans)
cgk debug
*     write(*,*)'gk: back from egrad_trace'
cgk end
      close(unit=aodens)
      IF(trans) CLOSE(UNIT=aodensa)

 1000 continue
cgk debug
*      write(*,*)'gk: leaving nadct_traces'
cgk end
      return
      end
