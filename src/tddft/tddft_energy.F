      logical function tddft_energy(rtdb)
c
c $Id: tddft_energy.F,v 1.4 2002-04-10 19:17:18 sohirata Exp $
c 
c Main routine for TDDFT excitation energy calculations.
c
c Capability includes:
c (1) CIS, TDHF, TDDFT, Tamm-Dancoff TDDFT,
c (2) Spin-restricted singlets and triplets,
c (3) Spin-unrestricted doublets etc.,
c (4) Transition moments and oscillator strengths,
c (5) Characterizes symmetry of excited states,
c (6) Uses Davidson's trial vector algorithm.
c
c See, e.g., Bauernschmitt & Ahlrichs, CPL 256, 454 (1996);
c Hirata & Head-Gordon, CPL 302, 375 (1999), for formalisms.
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
c NOTE ON PARALLELIZATION: parallelized portion of this code
c is almost entirely hidden in shell_fock_build.F and GA routines.
c While we forbid nodes other than zero to print any routine
c outputs, the only other parts that explicitly refer to node id are
c tddft_epsilon and tddft_residuals, where we add or devide by
c (epsilon_a-epsilon_i) by node zero only.
c
c CAUTION: numfunc=30 hardwired from "cdft.fh".  It is nicer if it
c is read from runtime database.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer rtdb             ! Runtime database handle
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      logical oskel,oadapt     ! Symmetry control parameters
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer numfunc          ! Number of functionals implemented
      parameter (numfunc=30)   ! Must be consistent with "cdft.fh"
      character*255 movecs_out ! MO vector file name
      character*40 xname(numfunc)       ! Names of X functionals
      character*40 cname(numfunc)       ! Names of C functionals
      character*40 xcname(numfunc)      ! Names of XC functionals
      double precision xfac(numfunc)    ! Prefactors for X functionals
      double precision cfac(numfunc)    ! Prefactors for C functionals
      logical lxfac(numfunc)            ! True if local X functional
      logical nlxfac(numfunc)           ! True if gc X functional
      logical lcfac(numfunc)            ! True if local C functional
      logical nlcfac(numfunc)           ! True if gc C functional
      logical xccomb(numfunc)           ! True if composite XC functional
      logical lhf                       ! True if HF (no DFT)
      logical ldft                      ! True if pure DFT (no HF exch)
      logical lhybrid                   ! True if hybrid DFT and HF
      logical lhasdft                   ! True if it has DFT
      logical lhashf                    ! True if it has HF
      logical tda                       ! True if Tamm-Dancoff approximation
      logical singlet                   ! True if singlet excited states
      logical triplet                   ! True if triplet excited states
      integer nroots           ! Number of roots sought
      integer maxvecs          ! Maximum number of trial vectors
      integer maxiter          ! Maximum number of Davidson iterations
      integer target           ! Target excited state
      integer targetsym        ! Symmetry of target excited state
      double precision tol2e   ! 2-e integral tolerance
      integer nov(2)           ! Number of occupied virtual pairs
      integer g_trials(2)      ! GA handle for trial vectors (XY)
      integer g_prod_p(2)      ! GA handle for A+B product vectors (XY)
      integer g_prod_m(2)      ! GA handle for A-B product vectors (XY)
      integer g_movecs(2)      ! GA handle for MO vectors
      integer l_evl,k_evl      ! MA handles for eigenvalues (AB)
      integer l_occ,k_occ      ! MA handles for occupancies (AB)
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      double precision thresh  ! Convergence threshold for Davidson iteration
      double precision energy  ! Target excited state energy
      logical xc_active_save   ! Previous value for xc_active
      integer calc_type_save   ! Previous value for calc_type
      logical algorithm        ! =1 multiple fock build; =2: conventional
      integer mult             ! Ground state spin multiplicity
c
      logical dft_energy
      external dft_energy
c
      integer i
      integer memleft
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c -------------------------
c Ground-state HF/DFT first
c -------------------------
c
      call tddft_dftinit(rtdb,oskel,oadapt)
      if (.not.dft_energy(rtdb)) then
        tddft_energy=.false.
        return
      endif
      tddft_energy=.false.
c
c -------------------
c Print utility start
c -------------------
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'tddft')
c
c ---------------------
c Initialize parameters
c ---------------------
c
      call tddft_init(rtdb,geom,oskel,oadapt,ao_bas_han,
     1  nbf_ao,movecs_out,ipol,nocc,nmo,nfc,nfv,numfunc,
     2  xname,cname,xcname,xfac,cfac,lxfac,nlxfac,lcfac,nlcfac,xccomb,
     3  lhf,ldft,lhybrid,lhasdft,lhashf,
     4  tda,singlet,triplet,nroots,maxvecs,tol2e,thresh,
     5  maxiter,target,targetsym,algorithm,g_movecs,
     6  l_evl,k_evl,l_occ,k_occ,mult)
c
c ----------------------------------------------------
c Determine memory size and alter maxvecs if necessary
c ----------------------------------------------------
c
      call tddft_memory(rtdb,nroots,maxvecs,maxiter,
     1  ipol,nbf_ao,nocc,nmo,nfc,nfv,nov,tda,algorithm)
c
c ----------------------------------------------------
c Allocate global arrays for trial and product vectors
c ----------------------------------------------------
c
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nov(i),maxvecs,
     1    'trial vectors',-1,-1,g_trials(i)))
     2    call errquit('tddft_energy: failed to create g_trials',0)
        if (.not.ga_create(mt_dbl,nov(i),maxvecs,
     1    'product vectors',-1,-1,g_prod_p(i)))
     2    call errquit('tddft_energy: failed to create g_prod_p',0)
        if (.not.tda) then
          if (.not.ga_create(mt_dbl,nov(i),maxvecs,
     1      'product vectors',-1,-1,g_prod_m(i)))
     2      call errquit('tddft_energy: failed to create g_prod_m',0)
        endif
      enddo
      memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
      if (nodezero.and.util_print('trial vectors',print_high))
     1  write(LuOut,9010) memleft
 9010 format(2x,'Memleft=',i10,' doubles')
c
c -----------------------------------------------------
c Iterative solution of TDDFT matrix eigenvalue problem
c -----------------------------------------------------
c
c Here we first initialize 2-e integral evaluators
c
      call int_init(rtdb,1,ao_bas_han)
      call schwarz_init(geom,ao_bas_han)
      call scf_get_fock_param(rtdb,tol2e)
c
c Set parameters for fock_xc if lhasdft
c
      if (lhasdft) then
        if (.not.rtdb_get(rtdb,'bgj:xc_active',mt_log,1,
     1    xc_active_save)) xc_active_save=.false.
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type',mt_int,1,
     1    calc_type_save)) calc_type_save=0
        if (.not.rtdb_put(rtdb,'bgj:xc_active',mt_log,1,.true.))
     1    call errquit('tddft_energy: failed to set xc_active true',0)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,2))
     1    call errquit('tddft_energy: failed to set calc_type 2',0)
      endif
c
c Now we call tddft_davidson
c
      if (ipol.eq.1) then
c
c -------------------
c Restricted singlets
c -------------------
c
        if (singlet) then
c
c Generate initial guess trial vectors
c
          call tddft_iguess(rtdb,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1      nroots,maxvecs,g_trials,g_movecs,dbl_mb(k_evl),
     2      ntrials,nproducts)
          if (nodezero.and.util_print('trial vectors',print_high)) then
            do i=1,ntrials
              call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1          g_prod_p,g_prod_m,i,ntrials,nproducts,tda)
            enddo
          endif
c
c Davidson iteration
c
          if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1      call errquit('tddft_energy: failed to set triplet',0)
          call tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1      nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2      g_movecs,dbl_mb(k_evl),ntrials,nproducts,
     3      tda,oskel,.true.,.false.,tol2e,xfac(1),lhashf,thresh,
     4      maxiter,target,targetsym,algorithm,mult)
        endif
c
c -------------------
c Restricted triplets
c -------------------
c
        if (triplet) then
c
c Generate initial guess trial vectors
c
          call tddft_iguess(rtdb,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1      nroots,maxvecs,g_trials,g_movecs,dbl_mb(k_evl),
     2      ntrials,nproducts)
          if (nodezero.and.util_print('trial vectors',print_high)) then
            do i=1,ntrials
              call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1          g_prod_p,g_prod_m,i,ntrials,nproducts,tda)
            enddo
          endif
c
c Davidson iteration
c
          if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.true.))
     1      call errquit('tddft_energy: failed to set triplet',0)
          call tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1      nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2      g_movecs,dbl_mb(k_evl),ntrials,nproducts,
     3      tda,oskel,.false.,.true.,tol2e,xfac(1),lhashf,thresh,
     4      maxiter,target,targetsym,algorithm,mult)
        endif
c
c ------------
c Unrestricted
c ------------
c
      else
c
c Generate initial guess trial vectors
c
        call tddft_iguess(rtdb,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1    nroots,maxvecs,g_trials,g_movecs,dbl_mb(k_evl),
     2    ntrials,nproducts)
        if (nodezero.and.util_print('trial vectors',print_high)) then
          do i=1,ntrials
            call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1        g_prod_p,g_prod_m,i,ntrials,nproducts,tda)
          enddo
        endif
c
c Davidson iteration
c
        call tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1    nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2    g_movecs,dbl_mb(k_evl),ntrials,nproducts,
     3    tda,oskel,.false.,.false.,tol2e,xfac(1),lhashf,thresh,
     4    maxiter,target,targetsym,algorithm,mult)
      endif
c
c And we close the 2-e integral evaluators
c
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()
c
c -----------------
c Print utility end
c -----------------
c
      call util_print_pop
c
c ------
c Return
c ------
c
      if (lhasdft) then
        if (.not.rtdb_put(rtdb,'bgj:xc_active',mt_log,1,
     1    xc_active_save))
     2    call errquit('tddft_energy: failed to restore xc_active',0)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,
     1    calc_type_save))
     2    call errquit('tddft_energy: failed to restore calc_type',0)
      endif
      if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1  call errquit('tddft_energy: failed to unset triplet',0)
      if (rtdb_get(rtdb,'tddft:energy',mt_dbl,1,energy))
     1  tddft_energy=.true.
      call tddft_finalize(geom,ao_bas_han,ipol,g_movecs,l_evl,k_evl,
     1  l_occ,k_occ)
      do i=1,ipol
        if (.not.ga_destroy(g_trials(i))) call errquit
     1    ('tddft_energy: failed to destroy g_trials',0)
        if (.not.ga_destroy(g_prod_p(i))) call errquit
     1    ('tddft_energy: failed to destroy g_prod_p',0)
        if (.not.tda) then
          if (.not.ga_destroy(g_prod_m(i))) call errquit
     1      ('tddft_energy: failed to destroy g_prod_q',0)
        endif
      enddo
      return
      end
