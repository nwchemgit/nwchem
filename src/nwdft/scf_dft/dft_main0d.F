      logical function dft_main0d(rtdb)
c     
c     $Id: dft_main0d.F,v 1.43 2000-04-22 01:05:45 edo Exp $
c     
      implicit none
c
      integer rtdb              ! rtdb handle [input]
      Logical LResult, IOLGC, oprint, oprint_time, oprint_ioinfo
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
c     
#if defined(DFT_TIMER)
#include "timer.fh"
#endif
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_dens(2), g_evecs(2), g_vxc(4)
      integer g_svecs, g_xcinv, g_oi
      integer g_s, g_fock, g_tmpb, g_tmp
      integer nvxc
      integer me, nproc, info
      double precision Enuc, Etold
      integer iter
      integer iVxc_opt, lmiss, ivcoul_opt
      integer isvals, lsvals,  idum, i
      double precision toll_s
      logical dft_scf,l1ecache
      external dft_scf
      logical dft_scf_so 
      external dft_scf_so 
      character*80 theory
c     
#if defined(DFT_TIMER)
      call t_stat(11,-2)
      call dft_second(tp1)
#endif
c     
      dft_main0d = .false.
c     
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('dft info', print_default)
      oprint_time = util_print('dft timings', print_high)
      oprint_ioinfo = util_print('io info', print_debug)
      call schwarz_init(geom, AO_bas_han)
      l1ecache=.false.
      if ((nbf_ao*nbf_ao*3*8)/ga_nnodes() .lt. 2000000) then! Cache if < 2MB/node
         l1ecache=.true.
         call int_1e_cache_ga(AO_bas_han, oskel)
      endif
      IOLGC = .TRUE.
      if (noio.eq.1) IOLGC = .FALSE.
c     
c     Evaluate the nuclear replusion energy.
c     
      LResult = geom_nuc_rep_energy(geom, Enuc)
c     
c     Determine the characteristics of the AO Gaussian basis set.
c     This uses the old NWdft data structures, rather than tthe rtdb
c     
      if (.not. bas_numbf(AO_bas_han, nbf_ao)) then
         call errquit('Exiting from dft_main0d',1)
      endif 
c
c     Set orbital overlap tolerance
c     
      if (.not. rtdb_get(rtdb, 'dft:toll_s', MT_DBL, 1, toll_s))then
         toll_s = 1.d-6
         if (.not. rtdb_put(rtdb, 'dft:toll_s', MT_DBL, 1, toll_s))
     &      call errquit('dft_main0d: failed to put toll_s in rtdb', 0)
      endif
c
      iter = 0
      Etold = 0.D0
c     
c     Select method to compute the XC potential.
c     
      if (XC_bas_han.eq.-99999)then
         iVxc_opt = 0
         XCFIT    = .FALSE.
      else
         iVxc_opt = 1
         XCFIT    = .TRUE.
      endif

      if (XCFIT)then
c     
c        Determine the characteristics of the XC Gaussian basis set.
c     
         if (.not. bas_numbf(XC_bas_han, nbf_xc)) then
            call errquit('Exiting from scfdriver',4)
         endif 
         lmiss = 1
         if (IOLGC)then
*ga:1:0
            if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'XCinv', 
     &                          nbf_xc, 0, g_xcinv)) 
     &         call errquit('dft_main0d: error creating g_xcinv',0)
            call ga_zero(g_xcinv)
            if (me.eq.0) then
c
c              try to read old XC matrix on disk; use if appropriate
c
               call dft_invio('XC', g_xcinv, nbf_xc, 'READ', lmiss)
               if (lmiss.ne.0)then
                  if (oprint_ioinfo)write(LuOut,*)
     &            ' old XC READ miss; constructing new matrix. '
               endif
            endif
cedo            call ga_sync
            call ga_brdcst(Msg_XC_Ex, lmiss, mitob(1), 0)
         endif
         if (lmiss.eq.1) then
c     
c           Compute the matrix inverse of the 2-ctr OIs.
c     
            if (.not.IOLGC)then
*ga:1:0
               if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'XCinv', 
     &                             nbf_xc, 0, g_xcinv)) 
     &            call errquit('dft_main0d: error creating g_xcinv',0)
               call ga_zero(g_xcinv)
            endif
*ga:1:0
            if (.not. ga_create(mt_dbl, nbf_xc, nbf_xc, 'g_oi', 
     &                          nbf_xc, 0, g_oi)) 
     &         call errquit('dft_main0d: error creating g_oi',0)
            call ga_zero(g_oi)
            call int_1e_ga(xc_bas_han, XC_bas_han,
     .           g_oi, 'overlap', oskel)
#if defined(PARALLEL_DIAG)
            call ga_chol(g_oi, g_xcinv, info)
#else 
            call ga_chol_seq(g_oi, g_xcinv, info)
#endif
            if (info.ne.0)then
               if (me.eq.0)then
                  write(LuOut,*)' Problem in performing a Choleski '
                  write(LuOut,*)' decomposition of the 2-ctr overlap '
                  write(LuOut,*)' matrix using XC fitting basis. '
                  write(LuOut,*)' Attempting a diag/inverse. '
               endif
            endif
            if (info.eq.0) then
#if defined(PARALLEL_DIAG)
               call ga_inverse(g_xcinv, g_oi)
               call ga_dgemm('T', 'N', nbf_xc, nbf_xc, nbf_xc, 1.d0, 
     &              g_oi, g_oi, 0.d0, g_xcinv)
#else 
               call ga_copy(g_xcinv, g_oi)
               call ga_inv_seq(g_oi, g_xcinv)
#endif
            else
               call dft_invdiag(g_oi, g_xcinv, nbf_xc)
            endif
            if (.not. ga_destroy(g_oi)) call errquit
     &         ('dft_main0d: could not destroy g_oi', 0)
            if (me.eq.0.and.IOLGC) then
               lmiss = 0
               call dft_invio('XC', g_xcinv, nbf_xc, 'WRITE', lmiss)
               if (lmiss.ne.0)call errquit
     &         (' dft_main0d; dft_invio - abnormal write of XC ', 0)
            endif
         endif
         if (IOLGC)then
            if (.not. ga_destroy(g_xcinv)) call errquit
     &         ('dft_main0d: could not destroy g_xcinv', 0)
         endif
      endif
c
c     Done with XC setup stuff
c
c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(1) = '  XCinv   '
      call tcllct(t_pgm(1,1),tp2,tp1)
      call dft_second(tp1)
#endif
c     
c     Select method to compute Coulomb potential.
c     
      if (CD_bas_han.eq.-99999)then
         iVcoul_opt = 0
         CDFIT      = .FALSE.
      else
         iVcoul_opt = 1
         CDFIT      = .TRUE.
      endif

c      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2a'
      if (CDFIT)then
c     
c        Determine the characteristics of the CD Gaussian basis set.
c     
         if (.not. bas_numbf(CD_bas_han, nbf_cd)) then
            call errquit('Exiting from scfdriver',5)
         endif 
c     
c        Compute the matrix inverse of the CD 2-ctr ERIs.
c     
*ga:1:0
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERI', 
     &                       0, nbf_cd, g_2ceri)) 
     &      call errquit('dft_main0d: error creating g_2ceri',0)
!      g_2ceri = ga_create_atom_blocked(geom, CD_bas_han, 'CD2cer')
         call ga_zero(g_2ceri)
         call dft_get2eri(CD_bas_han, g_2ceri)
*ga:1:0
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERInv', 
     &                        0, nbf_cd,g_cdinv)) 
     &      call errquit('dft_main0d: error creating g_cdinv',0)
!      g_cdinv = ga_create_atom_blocked(geom, CD_bas_han, 'CD2ceriv')
         call ga_zero(g_cdinv)
         lmiss = 1
c         if (lmiss.eq.1) then
            call ga_zero(g_cdinv)
            info = 0
c            do i = 0, nproc
c               if (me.eq.i) call dft_dump_info(me)
c            enddo
#if defined(PARALLEL_DIAG)
            call ga_chol(g_2ceri, g_cdinv, info)
#else 
            call ga_chol_seq(g_2ceri, g_cdinv, info)
#endif
            if (info.ne.0)then
               if (me.eq.0)then
                  write(LuOut,*)' Problem in performing a Choleski '
                  write(LuOut,*)' decomposition of the 2-ctr ERI '
                  write(LuOut,*)' matrix using CD fitting basis. '
                  write(LuOut,*)' Attempting a diag/inverse. '
c      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2g'
c                  call errquit(' temp stop this should not happen',0)
               endif
            endif
            if (info.eq.0) then
               g_tmpb = g_2ceri
#if defined(PARALLEL_DIAG)
               call ga_inverse(g_cdinv, g_tmpb)
               call ga_dgemm('T', 'N', nbf_cd, nbf_cd, nbf_cd, 1.d0, 
     &              g_tmpb, g_tmpb, 0.d0, g_cdinv)
#else 
               call ga_copy(g_cdinv, g_tmpb)
               call ga_inv_seq(g_tmpb, g_cdinv)
#endif
            else
               call dft_invdiag(g_2ceri, g_cdinv, 
     &                          nbf_cd)
            endif
            call ga_zero(g_2ceri)
            call dft_get2eri(CD_bas_han, g_2ceri, geom)
            if (IOLGC.and.(me.eq.0)) then
               lmiss = 0
               call dft_invio('CDI', g_cdinv, nbf_cd, 'WRITE', lmiss)
               if (lmiss.ne.0)call errquit
     &         (' dft_main0d; dft_invio - abnormal write of CDI ', 0)
               if (oprint_ioinfo)then
                 write(LuOut,*)' g_cdinv = ', g_cdinv
                 call ga_print(g_cdinv)
               endif
               lmiss = 0
               call dft_invio('CD', g_2ceri, nbf_cd, 'WRITE', lmiss)
               if (lmiss.ne.0)call errquit
     &         (' dft_main0d; dft_invio - abnormal write of CD ', 0)
               if (oprint_ioinfo)then
                 write(LuOut,*)' g_2ceri = ',  g_2ceri
                 call ga_print(g_2ceri)
               endif
            endif
cedo            call ga_sync
c         endif
         if (IOLGC) then
c            call gadest(g_cdinv)
            if (.not. ga_destroy(g_cdinv)) call errquit
     &         ('dft_main0d: could not destroy g_xcinv', 0)
c            call gadest(g_2ceri)
            if (.not. ga_destroy(g_2ceri)) call errquit
     &         ('dft_main0d: could not destroy g_xcinv', 0)
         endif
cedo         call ga_sync
      endif
c
c     Done with CD setup stuff
c
cedo      call ga_sync
c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(2) = '  CDinv   '
      call tcllct(t_pgm(1,2),tp2,tp1)
      call dft_second(tp1)
#endif
c     
c     Compute the AO overlap matrix.
c     
      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call ga_zero(g_s)
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'ovl eig vals', lsvals, 
     &   isvals))
     &   call errquit('dft_main0d: cannot allocate ovl eig vals', 0)
      call dfill(nbf_ao, 0.0d0, dbl_mb(isvals), 1)
      call int_1e_ga(AO_bas_han, AO_bas_han, g_s, 'overlap', oskel)
      if (oskel)call sym_symmetrize(geom, AO_bas_han, .false., g_s)
c     
c     
c     orthog BS
c     
c     Diag Overlap Matrix
c     
      g_svecs = ga_create_atom_blocked(geom, AO_bas_han, 
     &                                 'AO ovl eig vecs')
      call ga_zero(g_svecs)
#if defined(PARALLEL_DIAG)
      call ga_diag_std(g_s, g_svecs, Dbl_MB(isvals))
#else 
      call ga_diag_std_seq(g_s, g_svecs, Dbl_MB(isvals))
#endif
      do i = me, nbf_ao-1, nproc
         if (Dbl_MB(isvals+i).lt.toll_s)  
     &      write(LuOut,*) ' singular eigenvalue in S matrix ',
     &      i, toll_s, Dbl_MB(isvals+i)
      enddo
c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(3) = '  S Diag   '
      call tcllct(t_pgm(1,3),tp2,tp1)
      call dft_second(tp1)
            ltime = 7
      if(me.eq.0) call ffflush(ltime)
#endif
c     
cedo      call ga_sync
      g_evecs(1) = ga_create_atom_blocked(geom, AO_bas_han, 
     &                                    'alpha evecs')
      call ga_zero(g_evecs(1))
      g_dens(1) = ga_create_atom_blocked(geom, AO_bas_han, 
     &                                   'alpha density')
      call ga_zero(g_dens(1))
      if (ipol.ne.1)then
         g_dens(2) = ga_create_atom_blocked(geom, AO_bas_han, 
     &                                      'beta density')
         call ga_zero(g_dens(2))
         g_evecs(2) = ga_create_atom_blocked(geom, AO_bas_han, 
     &                                       'beta evecs')
         call ga_zero(g_evecs(2))
         g_Vxc(2) = ga_create_atom_blocked(geom, AO_bas_han, 'Vxcs')
         call ga_zero(g_Vxc(2))
      endif 
      nVxc = ipol
      g_fock = ga_create_atom_blocked(geom, AO_bas_han, 'fock matrix')
      call ga_zero(g_fock)
      g_vxc(1) = g_fock
c      lkeeps=.false.
c
c     keeps hardwired right now ... better something smarter
c
      lkeeps=.false.
      if(lkeeps) then
c
c     get s-12 and s+12
c
cnso         call diis_gets12pm(geom,ao_bas_han,toll_s,dbl_mb(isvals), 
cnos     ,        g_svecs,g_s12m, g_s12p)
cnos         call ga_copy(g_s12m,g_evecs(1))
      else
c     
c     Begin iteration sequence with S**(-1/2) in g_evecs arrays
c     
      g_tmp = ga_create_atom_blocked(geom, AO_bas_han, 'temp')
         call diis_bld12(toll_s, dbl_mb(isvals), g_svecs, g_evecs(1), 
     &        g_tmp,  2)
         if (.not. ga_destroy(g_tmp)) call errquit
     &      ('dft_main0d: could not destroy g_tmp', 0)
      endif
      if (ipol.gt.1) call ga_copy(g_evecs(1), g_evecs(2))
c     
      if (oprint_time)
     &     call dft_tstamp(' Before call to DFT_SCF. ')
c     
c     if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2'
c      dft_main0d = dft_scf(rtdb, Etold, Enuc, iVcoul_opt, iVxc_opt, 
c     &     iter, g_dens, g_evecs, g_vxc, g_fock, 
c     &     g_svecs, isvals, g_xcinv, g_s)
c 
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0)
      if (theory .eq. 'dft') then
         dft_main0d = dft_scf(rtdb, Etold, Enuc, iVcoul_opt, iVxc_opt,
     &        iter, g_dens, g_evecs, g_vxc, g_fock,
     &        g_svecs, isvals, g_xcinv, g_s)
      else if (theory .eq. 'sodft') then
         dft_main0d = dft_scf_so(rtdb, Etold, Enuc, iVcoul_opt, 
     &        iVxc_opt, iter, g_dens, g_evecs, g_vxc, g_fock, 
     &        g_svecs, isvals, g_xcinv)
      endif 
cso 
c     
c     Clean up memory usage.
c     
      if (.not.ma_pop_stack(lsvals))
     &     call errquit('dft_main0d: cannot pop stack', 12)
cnos      if (.not.CDFIT)then
         if (.not. ga_destroy(g_s)) call errquit
     &      ('dft_main0d: could not destroy g_s', 0)
cnos      endif
      if (me.eq.0) then
         if (CDFIT.and.IOLGC) then
            call dft_invio('CD', idum, 1, 'KILL', lmiss)
            call dft_invio('CDI', idum, 1, 'KILL', lmiss)
         endif
         if (XCFIT.and.IOLGC) then
            call dft_invio('XC', idum, 1, 'KILL', lmiss)
         endif
      endif
      if (CDFIT.and.(.not.IOLGC)) then
         if (.not. ga_destroy(g_cdinv)) call errquit
     &      ('dft_main0d: could not destroy g_cdinv', 0)
         if (.not. ga_destroy(g_2ceri)) call errquit
     &      ('dft_main0d: could not destroy g_2ceri', 0)
      endif
      if (XCFIT.and.(.not.IOLGC))then
         if (.not. ga_destroy(g_xcinv)) call errquit
     &      ('dft_main0d: could not destroy g_xcinv', 0)
      endif
      if (.not. ga_destroy(g_fock)) call errquit
     &   ('dft_main0d: could not destroy g_fock', 0)
      if (ipol.eq.2)then
         if (.not. ga_destroy(g_vxc(2))) call errquit
     &      ('dft_main0d: could not destroy g_vxc(2)', 0)
      endif
      if (ipol.ne.1)then
         if (.not. ga_destroy(g_evecs(2))) call errquit
     &      ('dft_main0d: could not destroy g_evecs(2)', 0)
         if (.not. ga_destroy(g_dens(2))) call errquit
     &      ('dft_main0d: could not destroy g_dens(2)', 0)
      endif
      if (.not. ga_destroy(g_evecs(1))) call errquit
     &   ('dft_main0d: could not destroy g_evecs(1)', 0)
      if (.not. ga_destroy(g_dens(1))) call errquit
     &   ('dft_main0d: could not destroy g_dens(1)', 0)
      if (.not. ga_destroy(g_svecs)) call errquit
     &   ('dft_main0d: could not destroy g_svecs', 0)
      if(l1ecache) call int_1e_uncache_ga()
      call schwarz_tidy
cedo      call ga_sync
c     
#if defined(DFT_TIMER)
Cedo      if (oprint_time) then
         if (me.eq.0)then
C            ltime = LuOut
            ltime = 7
            write(ltime,*)
         endif  
         pgmname = ' SCFDRIVER'
         call t_stat(11,2)
cedo      endif
#endif
c     
      return
      end
