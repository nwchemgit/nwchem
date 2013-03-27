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
     +                              enrep, nexc)
        call uhf_tidy(rtdb)
      else
        call rohf_init(rtdb)
        lsuccess = dft_roks_nr_solve(rtdb, energy, eone, etwo,
     +                               exc, enrep, nexc)
        call rohf_tidy
      endif
      call grid_cleanup(.true.)
c
      dft_time = dft_time + util_cpusec()
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
        write(luout,220) rho_n
        write(luout,230) dft_time
      endif
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
          if (util_print('final vectors analysis',print_high)) then
            ilo = 1
            ihi = nmo
          endif
          if (ispin.eq.1) then
            blob = 'DFT Final Alpha Molecular Orbital Analysis'
          else
            blob = 'DFT Final Beta Molecular Orbital Analysis'
          endif
          call movecs_fix_phase(g_movecs(ispin))
          call movecs_print_anal(basis, ilo, ihi, 0.15d0,
     &         g_movecs(ispin),blob,
     &         .true., dbl_mb(k_eval+(ispin-1)*nbf),
     &         oadapt, int_mb(k_irs+(ispin-1)*nbf),
     &         .true., dbl_mb(k_occ+(ispin-1)*nbf))
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
        if (util_print('final vectors analysis',print_high)) then
          ilo = 1
          ihi = nmo
        endif
        call movecs_fix_phase(g_movecs)
        blob = 'DFT Final Molecular Orbital Analysis'
        call movecs_print_anal(basis, ilo, ihi, 0.15d0,
     &         g_movecs,blob,
     &         .true., dbl_mb(k_eval),
     &         oadapt, int_mb(k_irs),
     &         .true., dbl_mb(k_occ))
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
 220  format(' Numeric. integr. density =', f22.12/)
 230  format('     Total iterative time =', f9.1,'s'//)
      end
