c
c     /* Main DFT driver */
c
      logical function wrt_dft_aoints(rtdb)
c     
c     $Id: 
c     
      implicit none
#include "errquit.fh"
c
      integer rtdb              ! rtdb handle [input]
      Logical IOLGC, oprint, oprint_time, oprint_ioinfo
c     
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "global.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "dftps.fh"
#include "bgj.fh"
#include "dft_scaleig.fh"
!#ifdef SCALAPACK
!      integer ga_cholesky,ga_llt_i
!      external ga_cholesky,ga_llt_i
!#endif
c     
#include "zora.fh"
#include "rel_consts.fh"
c
#include "case.fh"
#include "geomP.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_xcinv, g_oi
*     integer g_2ceri, gcdinv
      integer g_tmpb, g_tmp
      integer nvxc
      integer me, nproc, cinfo
      integer iter, nconstr
      integer iVxc_opt, lmiss, ivcoul_opt
      integer lsvals, isvals,  idum,n_dep
      double precision toll_s
      logical dft_scf,atmblk, lcdft
      external dft_scf
      logical dft_scf_so 
      external dft_scf_so 
      logical status
      logical ldirect
c     sifs
      integer sizeof_int, nmap, l_map, k_map, szintarr, szlabs
      integer k_lab
      integer l_lab
      integer  nmax
      parameter(nmap=2)
      integer*4 aoints
      parameter(aoints=23)
      integer*4 ifmt1, ibvtyp, n1max, ierr, ifmt, n2max
      INTEGER l1rec, l2rec
      integer lrecal
      integer*4 ibitv
      integer*4 mxengy
      parameter (mxengy=1)
      real*8 energy(mxengy)
      integer*4 ninfo, iunit2, aoint2
      parameter(ninfo=5)
      integer*4 info(ninfo)  ! not to be confused with NWChem info
      real*8 tol2e
      parameter(tol2e=1d-10) ! integrate with user parameters...

cgk debug
      integer mem1, max1e, mem2, max2e
cgk end
      ierr=0 !initialize
       

      wrt_dft_aoints = .false.
c     
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('dft cinfo', print_default)
      oprint_time = util_print('dft timings', print_high)
      oprint_ioinfo = util_print('io info', print_debug)
      if(ga_nodeid().eq.0.and.
     .     util_print('geometry', print_high)) then
        if(.not.geom_print(geom)) call errquit(
     .       ' wrt_dft_aoints: geom print failed',0, GEOM_ERR)
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
         call errquit('Exiting from wrt_dft_aoints',1, BASIS_ERR)
      endif 
*     iter = 0
      if(nproc.gt.1) then
        write(6,'(a,a)')'WARNING: Writing aoints, calculating 1e ints ',
     &    'with 1 CPU'
        write(6,'(a,a)')'WARNING: we are currently supporting 1 ',
     &    'process for aoints'
      endif
c
c     open the sifs file
c
      if(ga_nodeid().eq.0) then
        open(unit=aoints,status='unknown',form='unformatted',
     &    file='aoints')
        rewind aoints
*       open(unit=23,status='unknown',form='formatted',
*    &    file='testfile')
*       write(23,*) 'testing testfile'
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
         call errquit('wrt_dft_aoints: XC basis fitting not implemented'
     &,     0,CAPMIS_ERR)
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
*     &            call errquit('wrt_dft_aoints: error creating g_xcinv',
*     &              0,GA_ERR)
*               call ga_zero(g_xcinv)
**ga:1:0
*            if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'g_oi', 
*     &                          nbf_xc, 0, g_oi)) 
*     &         call errquit('wrt_dft_aoints: error creating g_oi',0,
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
*     &         ('wrt_dft_aoints: could not destroy g_oi', 0, GA_ERR)
*            if (me.eq.0.and.IOLGC) then
*               lmiss = 0
*               call dft_invio('XC', g_xcinv, nbf_xc, 'WRITE', lmiss)
*               if (lmiss.ne.0)call errquit
*     &         (' wrt_dft_aoints; dft_invio - abnormal write of XC ', 0,
*     &       DISK_ERR)
*            endif
*         if (IOLGC)then
*            if (.not. ga_destroy(g_xcinv)) call errquit
*     &         ('wrt_dft_aoints: could not destroy g_xcinv', 0, GA_ERR)
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
         call errquit('wrt_dft_aoints: charge density fitting not implem
     &ented', 0,CAPMIS_ERR)
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
*     &      call errquit('wrt_dft_aoints: error creating g_2ceri',0,
*     &       GA_ERR)
*         call ga_zero(g_2ceri)
*         call dft_get2eri(CD_bas_han, g_2ceri,oskel)
*         if (odftps) call pstat_off(ps_2cer)
*         if (oskel)call 
*     .        sym_symmetrize(geom,cd_bas_han,.false.,g_2ceri)
*         call ga_sync()
*         if (.not. ga_duplicate(g_2ceri, g_cdinv, 'CD 2cERInv'))
*     &    call errquit('wrt_dft_aoints: error creating g_cdinv',0, 
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
*               call errquit('wrt_dft_aoints: inverse failed ',0,0)
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
*     &         (' wrt_dft_aoints; dft_invio - abnormal write of CDI ', 0
*     &         ,DISK_ERR)
*               lmiss = 0
*               call dft_invio('CD', g_2ceri, nbf_cd, 'WRITE', lmiss)
*               if (lmiss.ne.0)call errquit
*     &         (' wrt_dft_aoints; dft_invio - abnormal write of CD ', 0,
*     &       DISK_ERR)
*            endif
*            if (oprint_ioinfo)then
*               call ga_print(g_cdinv)
*               call ga_print(g_2ceri)
*            endif
*c         endif
*         if (IOLGC) then
*            if (.not. ga_destroy(g_cdinv)) call errquit
*     &         ('wrt_dft_aoints: could not destroy g_xcinv', 0, GA_ERR)
*            if (.not. ga_destroy(g_2ceri)) call errquit
*     &         ('wrt_dft_aoints: could not destroy g_xcinv', 0, GA_ERR)
*         endif
*      endif
*c
c     Done with CD setup stuff
c
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'ovl eig vals', lsvals, 
     &   isvals))
     &   call errquit('wrt_dft_aoints: cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      call dfill(nbf_ao, 0.0d0, dbl_mb(isvals), 1)
c
c     Set orbital overlap tolerance
c     
      call util_tolls(rtdb,.false.,toll_s,n_dep,dbl_mb(isvals),nbf_ao)
c
c     set up SIFS buffer and other data for 1e ints.
c
      ibvtyp=0
      lrecal = 4096
C LB
c      write(*,*)'LB: before sifcfg(1), ierr=',ierr
c      WRITE(*,*)'n1max=',n1max,'n2max=',n2max,'ifmt1=',ifmt1,'ifmt',ifmt
c      WRITE(*,*)"LB, sizeof l1rec=",sizeof(l1rec)
c      WRITE(*,*)"LB, sizeof lrecal=",sizeof(lrecal)
C LB
      call sifcfg(1, lrecal, nbf_ao, ibvtyp, ifmt1, l1rec, n1max, ierr)
      lrecal = 4096
C LB
c      write(*,*)'LB: after sifcfg(1), ierr=',ierr
c      WRITE(*,*)'l1rec=',l1rec
C LB
      call sifcfg(2, lrecal, nbf, ibvtyp, ifmt, l2rec, n2max, ierr)
C LB
c      write(*,*)'LB: after sifcfg(2), ierr=',ierr
c      WRITE(*,*)'l2rec=',l2rec
C LB

      info(1) = 1 ! =FSPLIT; do not split aoints file
      info(2) = l1rec
      info(3) = n1max
      info(4) = l2rec ! supposed be something else but not ready yet
      info(5) = n2max ! supposed be something else but not ready yet
      nmax=max(n1max,n2max)

cgk debug
c      write(*,*)'gk: n1max = ', n1max
c      write(*,*)'gk: l1rec = ', l1rec
c      write(*,*)'gk: n2max = ', n1max
c      write(*,*)'gk: l2rec = ', l2rec
c      write(*,*)'gk: ifmt1 = ', ifmt1
c      write(*,*)'gk: ifmt = ', ifmt
cgk end
      if(ierr.ne.0) 
     &   call errquit('wrt_dft_aoints: error in sifcfg', ierr,
     &      MA_ERR)
*     if(n1max.le.max1e) 
*    &   call errquit('wrt_dft_aoints: sifcfg does not match max1e', 
*    &      max1e, MA_ERR)
c
c     Allocate memory for i4 arrays used in sifs routines in wrt_dft_aoints
c
c     get memory for integer map(npbf,2)
      sizeof_int=MA_sizeof(MT_INT, 1, MT_BYTE)
cgk debug
*     write (*,*) 'gk: MA_sizeof(MT_INT, 1, MT_BYTE) = ',sizeof_int
cgk end
      if(sizeof_int.eq.4) then
        szintarr=(nbf_ao+1)/2
        szlabs=(nmax+1)/2
      elseif(sizeof_int.eq.8) then
        szintarr=nbf_ao
        szlabs=nmax
      else
        call errquit('wrt_dft_aoints: unsupporte integer type size=',
     &    sizeof_int, MA_ERR)
      endif
      if (.not. MA_push_get(MT_Int,(szintarr+1)*nmap,
     &  'wrt_dft_aoints:map', l_map, k_map))
     $     call errquit('wrt_dft_aoints: ma failed for map', 
     &        szintarr*nmap, MA_ERR)
c     get memory for labels
      if (.not. MA_push_get(MT_Int,4*szlabs,'wrt_dft_aoints:ilab',l_lab,
     &        k_lab))
     $     call errquit('wrt_dft_aoints: ma failed for lab', 2*szlabs, 
     &        MA_ERR)

      call schwarz_tidy()
      call int_terminate()
*     call int_init(rtdb, 1, AO_bas_han)
*     call schwarz_init(geom, AO_bas_han)
*     call int_terminate()
*     call intd_init(rtdb, 1, AO_bas_han)
cgk debug
*     write(*,*)'gk: calling int_init'
cgk end
      call int_init(rtdb, 1, AO_bas_han)
cgk debug
*     write(*,*)'gk: back from int_init'
cgk end
      call int_mem_1e(max1e, mem1)
*     call int_mem(max1e, max2e, mem1, mem2)
cgk debug
*      write(*,*)'gk:  before int_1e_sifs, mem1 = ', mem1
*      write(*,*)'gk:  mem2 = ', mem2
*      write(*,*)'gk:  max1e = ', max1e
*      write(*,*)'gk:  max2e = ', max2e
*      call int_mem_print()
cgk end

cgk debug
*     write(*,*)'gk: calling int_1e_sifs'
*     write(*,*)'gk:    mxengy=', mxengy
*      write(*,*)'gk:    info=', info
cgk write

      call int_1e_sifs(AO_bas_han, aoints, energy, mxengy, 
     &   nbf_ao, nmap, int_mb(k_map),int_mb(k_map+szintarr*nmap+1),
     &   ifmt1, ibvtyp, ibitv, l1rec, n1max, int_mb(k_lab),
     &   ninfo, info)

cgk debug
*     write(*,*)'gk: back from int_1e_sifs'
*     write(*,*)'gk:  nbf_ao=',nbf_ao
cgk write

cgk debug
c      write(*,*)'gk: done with one electron integrals'
c      write(*,*)'gk: now let us do 2e ints.'
c      write(*,*)'LB: aoints=',aoints,'iunit2=',iunit2,'aoint2=',aoint2
c     &,'ierr=',ierr
c      WRITE(*,*)'LB: info =',info(1),info(2),info(3),info(4)
cgk end
c     open 2 electron itegrals file
      iunit2=aoints
      call sifo2f(aoints,iunit2,'aoints',info,aoint2,ierr)

      if(ierr.ne.0) 
     &   call errquit('wrt_dft_aoints: error in sifo2f', ierr,
     &      MA_ERR)

cgk debug
*     write(*,*)'gk: writing 2e ints to unit, aoint2 = ', aoint2
*     write(*,*)'gk: iunit2 = ', iunit2
cgk end

      call int_2e_sifs(geom, AO_bas_han, tol2e, oskel, aoint2, ninfo, 
     &     info, ifmt, ibvtyp, ibitv, int_mb(k_lab))

c
c     done with memory for sifs files
c
      if (.not. MA_pop_stack(l_lab)) call errquit
     $     ('wrt_dft_aoints: pop failed at l_lab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_map)) call errquit
     $     ('wrt_dft_aoints: pop failed at l_map', 0, GA_ERR)

      if(ga_nodeid().eq.0) then
        close(unit=aoints)
      endif

      wrt_dft_aoints=.true.

!      if (me.eq.0) then
!         if (CDFIT.and.IOLGC) then
!            call dft_invio('CD', idum, 1, 'KILL', lmiss)
!            call dft_invio('CDI', idum, 1, 'KILL', lmiss)
!         endif
!         if (XCFIT.and.IOLGC) then
!            call dft_invio('XC', idum, 1, 'KILL', lmiss)
!         endif
!      endif
!      if (CDFIT.and.(.not.IOLGC)) then
!         if (.not. ga_destroy(g_cdinv)) call errquit
!     &      ('wrt_dft_aoints: could not destroy g_cdinv', 0, GA_ERR)
!         if (.not. ga_destroy(g_2ceri)) call errquit
!     &      ('wrt_dft_aoints: could not destroy g_2ceri', 0, GA_ERR)
!      endif
!      if (XCFIT.and.(.not.IOLGC))then
!         if (.not. ga_destroy(g_xcinv)) call errquit
!     &      ('wrt_dft_aoints: could not destroy g_xcinv', 0, GA_ERR)
!         if (.not. ga_destroy(g_oi)) call errquit
!     &      ('wrt_dft_aoints: could not destroy g_oi', 0, GA_ERR)
!      endif
      return
      end
