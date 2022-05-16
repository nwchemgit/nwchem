c
c     == Main DFT quasi-Newton-Raphson driver ==
      logical function dft_cg_solve(rtdb)
c     
c     $Id$
c     
      implicit none
c
#include "errquit.fh"
c
      integer rtdb              ! rtdb handle [input]
      Logical LResult, IOLGC, oprint, oprint_time, oprint_ioinfo
c     
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cscf.fh"
#include "dft_conv.fh"
#include "global.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "dftps.fh"
#include "bgj.fh"
#include "cosmo.fh"
c
#ifdef SCALAPACK
#include "dft_scaleig.fh"
      integer ga_cholesky,ga_llt_i
      external ga_cholesky,ga_llt_i
#endif
c     
#include "zora.fh"
#include "rel_consts.fh"
#include "case.fh"
#include "../nwdft/grid/cgridfile.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_dens(2), g_evecs(2), g_vxc(4)
c
      integer g_dens_at(2)
c
      integer g_svecs, g_xcinv, g_oi
      integer g_s, g_fock, g_tmpb, g_tmp
      integer nvxc
      integer me, nproc, info
      double precision Enuc, Etold
      integer nconstr
      integer iVxc_opt, lmiss, ivcoul_opt
      integer isvals, lsvals,  idum,n_dep
      integer ilo, ihi, ispin, ipol, noc(2)
      integer ao_bas_han
      double precision toll_s
c
      logical int_normalize
      external int_normalize
      external grid_reopen
      logical grid_reopen
c
      logical lsuccess, ltried_scf, ltried_cgmin, lcgmin
      logical dft_scf,atmblk, lcdft, dft_cscf_equals_scftype
      external dft_scf, dft_cscf_equals_scftype
      logical dft_scf_so 
      external dft_scf_so 
      logical  dft_uks_nr_solve, dft_roks_nr_solve
      external dft_uks_nr_solve, dft_roks_nr_solve
c
      character*80 theory
      logical status
      logical ldirect
      character*32 pname
      character*255 blob
c
      double precision eval_pr_tol_lo, eval_pr_tol_hi
      parameter (eval_pr_tol_lo = -1.5d0, eval_pr_tol_hi=0.5)
c
      double precision eone, etwo, enrep, energy
      double precision dft_time
      double precision edisp, dum
      logical disp
      double precision ecosmo,gcds
      logical do_cosmo_smd
      double precision esol, egas
      double precision gstote, gspol
c
      double precision tanalyze
c
      logical fon
      logical  xc_chkdispauto
      external xc_chkdispauto
c
      nexc = 1
c
c     == Preliminaries ==
      dft_cg_solve = .false.
      pname = "dft_cg_solve:"
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('dft info', print_default)
      oprint_time = util_print('dft timings', print_high)
      oprint_ioinfo = util_print('io info', print_debug)
      lsuccess = .false.
      ltried_scf = .false.
      ltried_cgmin = .false.
c
c    cg not ready for fon
c
      if (.not.rtdb_get(rtdb, 'dft:fon', mt_log, 1, fon)) fon=.false.
      if (fon) call errquit(pname//' cg not ready for FON ',0,
     E     CAPMIS_ERR)
c
c     cg not ready for smear
c
      if (rtdb_get(rtdb, 'dft:smear_sigma', mt_dbl, 1,
     &     ssmear))then
         if(ssmear.gt.1d-9)
     c        call errquit(pname//' cg not ready for SMEAR ',0,
     E        CAPMIS_ERR)
      endif

c
c     read MO analysis threshold
      if (.not.rtdb_get(rtdb,'movecs:tanalyze',MT_DBL,1,tanalyze))
     &        tanalyze = 0.025d0
c
c     == get occupancies ==
      if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,noc))
     &   call errquit(pname//'Error getting occupancies',
     &                0,RTDB_ERR)
c
c     == Evaluate the nuclear repulsion energy ==
      LResult = geom_nuc_rep_energy(geom, Enuc)
      if (.not.LResult) 
     & call errquit(pname//'Error in nuclear repulsion',
     &                1,BASIS_ERR)
c
c     == Dispersion corrections ==
      Edisp = 0.0d0
      if (.not.rtdb_get(rtdb, 'dft:disp', mt_log, 1, disp))
     &   disp=.false.
c
      if(disp.or.xc_chkdispauto())
     &      call xc_vdw(rtdb,geom,Edisp,dum,'energy')
c
c     == needed for dftmp2 ==
      if (.not. rtdb_cput(rtdb,'scf:scftype', 1, scftype))
     &   call errquit(pname//'rtdb_cput scftype failed',
     &                0, RTDB_ERR)
      if (.not.dft_cscf_equals_scftype("UHF")) then
        if (.not.rtdb_put(rtdb,'scf:nclosed',mt_int,1,noc))
     &     call errquit(pname//'Error putting occupancies',
     &                  0,RTDB_ERR)
      endif
c
c     == Do quasi-Newton-Raphson DFT ==
      call dft_cg_init(rtdb)
      call dft_cg_guess(rtdb)
      dft_time = -util_cpusec()
      grid_written=grid_reopen(geom)
c     grid is not there yet
      grid_written=.false.
      if (dft_cscf_equals_scftype("UHF")) then
        call uhf_init(rtdb)
        lsuccess = dft_uks_nr_solve(rtdb, energy, eone, etwo, exc,
     +                              enrep, edisp, nexc)
        call uhf_tidy(rtdb)
      else
        call rohf_init(rtdb)
        lsuccess = dft_roks_nr_solve(rtdb, energy, eone, etwo,
     +                               exc, enrep, edisp, nexc)
        call rohf_tidy
      endif
      call grid_cleanup(.true.)
c
      dft_time = dft_time + util_cpusec()
c
c     get the cosmo energy
      if(cosmo_on.and.cosmo_phase.eq.2) then
         if (.not.rtdb_get(
     $     rtdb,'cosmo:energy',mt_dbl,1,ecosmo))
     $     call errquit(
     $     'dft_cg_solve: rtdb get failed for ecosmo',
     $     0,rtdb_err)
         if (.not. rtdb_get(
     $     rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $     call errquit(
     $     'dft_cg_solve: cannot get do_cosmo_smd from rtdb',
     $     0,rtdb_err)
         if (do_cosmo_smd) then
           if (.not. rtdb_get
     $       (rtdb,'cosmo:gcds',mt_dbl,1,gcds))
     $       call errquit(
     $       'dft_cg_solve: cannot get gcds from rtdb',
     $       0,rtdb_err)
           if(.not.rtdb_get
     $       (rtdb,'cosmo:gspol',mt_dbl,1,gspol))
     $       call errquit(
     $       'dft_cg_solve: cannot get gspol from rtdb',
     $       0,rtdb_err)
         endif ! do_cosmo_smd
      endif
c
      if (ga_nodeid().eq.0) then
        if (.not.lsuccess) then
          write(luout,*)
          call util_print_centered(luout,
     +         'Calculation failed to converge', 20, .true.)
          write(luout,*)
        endif
        if (nexc.le.1)then
          write(LuOut,200)energy,eone,etwo,exc(1),enuc
        else
          write(LuOut,210)energy,eone,etwo,exc(1),exc(2),enuc
        end if
c
        if (abs(edisp).gt.0.0d0) write(luout,215) edisp
c
        if (cosmo_on.and.cosmo_phase.eq.2) then
          if (do_cosmo_smd) then
           write(LuOut,216) ecosmo+gcds
          else
           write(LuOut,216) ecosmo
          end if
         endif
c
        write(luout,220) rho_n
        write(luout,230) dft_time
      endif
c
c     print cosmo, smd output
      if (.not. rtdb_get(rtdb,'cosmo:egas',mt_dbl,1,egas)) egas = 0.d0
      if (.not. rtdb_get(rtdb,'cosmo:esol',mt_dbl,1,esol)) esol = 0.d0
      gstote = esol - gcds
      if (me.eq.0) then
        if(cosmo_on.and.cosmo_phase.eq.2) then
           call cosmo_smd_output(egas, esol, do_cosmo_smd, gcds,
     $            gstote, gspol)
        endif
      endif ! me.eq.0
c
      if (dft_cscf_equals_scftype("UHF")) then
        ipol = 2
        do ispin = 1, ipol
          do ilo = 1,max(1,min(nalpha,nbeta)-10)
            if (dbl_mb(k_eval+(ispin-1)*nbf+ilo-1) .ge. eval_pr_tol_lo)
     &          goto 900
          enddo
 900      do ihi = min(max(nalpha,nbeta)+10,nmo),nmo
            if (dbl_mb(k_eval+(ispin-1)*nbf+ihi-1) .ge. eval_pr_tol_hi)
     &          goto 910
          enddo
          ihi = max(ihi-1,1)
 910      continue
          call movecs_fix_phase(g_movecs(ispin))
          if (util_print('final vectors analysis',print_default)) then
            ilo = 1
            ihi = nmo
          if (ispin.eq.1) then
            blob = 'DFT Final Alpha Molecular Orbital Analysis'
          else
            blob = 'DFT Final Beta Molecular Orbital Analysis'
          endif
c
          call movecs_print_anal(basis, ilo, ihi, tanalyze,
     &         g_movecs(ispin),blob,
     &         .true., dbl_mb(k_eval+(ispin-1)*nbf),
     &         oadapt, int_mb(k_irs+(ispin-1)*nbf),
     &         .true., dbl_mb(k_occ+(ispin-1)*nbf))
          endif
        enddo
c
      else  ! closed shell
c
        do ilo = 1,max(1,min(nalpha,nbeta)-10)
          if (dbl_mb(k_eval+ilo-1) .ge. eval_pr_tol_lo) goto 920
        enddo
 920    do ihi = min(max(nalpha,nbeta)+10,nmo),nmo
          if (dbl_mb(k_eval+ihi-1) .ge. eval_pr_tol_hi) goto 930
        enddo
        ihi = max(ihi-1,1)
 930    continue
        call movecs_fix_phase(g_movecs)
        if (util_print('final vectors analysis',print_default)) then
          ilo = 1
          ihi = nmo
        blob = 'DFT Final Molecular Orbital Analysis'
        call movecs_print_anal(basis, ilo, ihi, 0.15d0,
     &         g_movecs,blob,
     &         .true., dbl_mb(k_eval),
     &         oadapt, int_mb(k_irs),
     &         .true., dbl_mb(k_occ))
        endif
       endif
c
c      == wrap up ==
       call dft_cg_answer(rtdb,energy,lsuccess)
       call dft_cg_tidy(rtdb)
       dft_cg_solve = lsuccess
c
      return
c
 200  format(//
     &     '         Total DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '    Exchange-Corr. energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12/)
c
 210  format(//
     &     '         Total DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '          Exchange energy =', f22.12/
     &     '       Correlation energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12/)
c
 215  format('    Dispersion correction =', f22.12/)
 216  format('             COSMO energy =', f22.12/)
c
 220  format(' Numeric. integr. density =', f22.12/)
 230  format('     Total iterative time =', f9.1,'s'//)
c
      end
