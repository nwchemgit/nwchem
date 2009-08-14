      logical function tddft_energy(rtdb)
c
c $Id: tddft_energy.F,v 1.13 2008-12-06 23:10:14 niri Exp $
c 
c Main routine for TDDFT excitation energy calculations.
c
c Capability includes:
c (1) CIS, TDHF, TDDFT, Tamm-Dancoff TDDFT,
c (2) Spin-restricted singlets and triplets,
c (3) Spin-unrestricted doublets etc.,
c (4) Transition moments and oscillator strengths,
c (5) Characterizes symmetry of excited states,
c (6) Uses Davidson's trial vector algorithm,
c (7) Spin (S**2) expectation values,
c (8) Incore/disk-based multiple/single tensor contraction,
c (9) Numerical gradients, optimization, and frequencies.
c
c See, e.g., Bauernschmitt & Ahlrichs, CPL 256, 454 (1996);
c Hirata & Head-Gordon, CPL 302, 375 (1999), for formalisms.
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
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
#include "dra.fh"
#include "xc.fh"
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
      character*255 movecs_out ! MO vector file name
      character*255 civecs_tmp ! CI vector file name
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
      character*4 targetsym    ! Symmetry of target excited state
      logical symmetry         ! Symmetry restriction
      double precision tol2e   ! 2-e integral tolerance
      integer nov(2)           ! Number of occupied virtual pairs
      integer g_trials(2)      ! GA handle for trial vectors (XY)
      integer g_prod_p(2)      ! GA handle for A+B product vectors (XY)
      integer g_prod_m(2)      ! GA handle for A-B product vectors (XY)
      integer d_trials(2)      ! DRA handle for trial vectors (XY)
      integer d_prod_p(2)      ! DRA handle for A+B product vectors (XY)
      integer d_prod_m(2)      ! DRA handle for A-B product vectors (XY)
      integer g_buffer1(2)     ! GA handle for a buffer vector for disk algorithms
      integer g_buffer2(2)     ! GA handle for a buffer vector for disk algorithms
      integer g_movecs(2)      ! GA handle for MO vectors
      integer l_evl,k_evl      ! MA handles for eigenvalues (AB)
      integer l_occ,k_occ      ! MA handles for occupancies (AB)
      integer l_irs(2),k_irs(2)! MA handles for irreps (AB)
      integer ntrials          ! Current number of trial vectors
      integer nproducts        ! Current number of product vectors
      double precision thresh  ! Convergence threshold for Davidson iteration
      double precision energy  ! Target excited state energy
      logical xc_active_save   ! Previous value for xc_active
      integer calc_type_save   ! Previous value for calc_type
      integer algorithm        ! =1 incore-multiple fock build; =2: incore-single fock build
                               ! =3 disk-multiple fock build;   =4: disk-single fock build
      integer mult             ! Ground state spin multiplicity
c
      logical dft_energy
      external dft_energy
c
      integer i
      integer memleft
      integer dummy
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
     5  maxiter,target,targetsym,symmetry,algorithm,g_movecs,
     6  l_evl,k_evl,l_occ,k_occ,l_irs,k_irs,mult)
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
      if ((algorithm.eq.1).or.(algorithm.eq.2)) then
        do i=1,ipol
          if (.not.ga_create(mt_dbl,nov(i),maxvecs,
     1      'trial vectors',-1,maxvecs,g_trials(i)))
     2      call errquit('tddft_energy: failed to create g_trials',0,
     &       GA_ERR)
          if (.not.ga_duplicate(g_trials(i),g_prod_p(i),
     1      'product vectors'))
     2      call errquit('tddft_energy: failed to create g_prod_p',0,
     &       GA_ERR)
          if (.not.tda) then
            if (.not.ga_duplicate(g_trials(i),g_prod_m(i),
     1        'product vectors'))
     2        call errquit('tddft_energy: failed to create g_prod_m',0,
     &       GA_ERR)
          endif
        enddo
c
c -----------------------------------------------------------
c Allocate disk resident arrays for trial and product vectors
c -----------------------------------------------------------
c
      else if ((algorithm.eq.3).or.(algorithm.eq.4)) then
        do i=1,ipol
          if (.not.ga_create(mt_dbl,nov(i),1,
     1      'buffer vector 1',-1,1,g_buffer1(i)))
     2      call errquit('tddft_energy: failed to create g_buffer1',0,
     &       GA_ERR)
          if (.not.ga_duplicate(g_buffer1(i),g_buffer2(i),
     1      'buffer vector 2'))
     2      call errquit('tddft_energy: failed to create g_buffer2',0,
     &       GA_ERR)
        enddo
c
c Six DRA's for 10 GB of total disk space and 10 MB of memory buffer
c see also tddft_memory.F when the following line is to be modified
c
        dummy=dra_init(6,1.0d10,1.0d10,1.0d7)
        call util_file_name('tvector1',.true.,.false.,civecs_tmp)
        dummy=dra_create(mt_dbl,nov(1),maxvecs,
     1    'trial vectors',civecs_tmp,DRA_RW,-1,maxvecs,d_trials(1))
        call util_file_name('pvector1',.true.,.false.,civecs_tmp)
        dummy=dra_create(mt_dbl,nov(1),maxvecs,
     1    'product vectors',civecs_tmp,DRA_RW,-1,maxvecs,d_prod_p(1))
        if (.not.tda) then
          call util_file_name('mvector1',.true.,.false.,civecs_tmp)
          dummy=dra_create(mt_dbl,nov(1),maxvecs,'product vectors',
     1      civecs_tmp,DRA_RW,-1,maxvecs,d_prod_m(1))
        endif
        if (ipol.eq.2) then
          call util_file_name('tvector2',.true.,.false.,civecs_tmp)
          dummy=dra_create(mt_dbl,nov(2),maxvecs,
     1      'trial vectors',civecs_tmp,DRA_RW,-1,maxvecs,d_trials(2))
          call util_file_name('pvector2',.true.,.false.,civecs_tmp)
          dummy=dra_create(mt_dbl,nov(2),maxvecs,
     1      'product vectors',civecs_tmp,DRA_RW,-1,maxvecs,d_prod_p(2))
          if (.not.tda) then
            call util_file_name('mvector2',.true.,.false.,civecs_tmp)
            dummy=dra_create(mt_dbl,nov(2),maxvecs,'product vectors',
     1        civecs_tmp,DRA_RW,-1,maxvecs,d_prod_m(2))
          endif
        endif
      else
        call errquit('tddft_energy: unknown algorithm',0, INPUT_ERR)
      endif
      memleft=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
      if (nodezero.and.util_print('trial vectors',print_high))
     1  write(LuOut,9010) memleft
 9010 format(2x,'Memleft=',i10,' doubles')
c
c -----------------------------------------------------
c Iterative solution of TDDFT matrix eigenvalue problem
c -----------------------------------------------------
c
c Set parameters for fock_xc if lhasdft
c
      if (lhasdft) then
        if (.not.rtdb_get(rtdb,'bgj:xc_active',mt_log,1,
     1    xc_active_save)) xc_active_save=.false.
        if (.not.rtdb_get(rtdb,'fock_xc:calc_type',mt_int,1,
     1    calc_type_save)) calc_type_save=0
        if (.not.rtdb_put(rtdb,'bgj:xc_active',mt_log,1,.true.))
     1    call errquit('tddft_energy: failed to set xc_active true',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,2))
     1    call errquit('tddft_energy: failed to set calc_type 2',0,
     &       RTDB_ERR)
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
          call tddft_iguess(rtdb,geom,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1      nroots,maxvecs,g_trials,g_buffer1,g_buffer2,d_trials,
     2      g_movecs,dbl_mb(k_evl),ntrials,nproducts,algorithm,
     3      targetsym,symmetry,l_irs,k_irs,tda)
          if (util_print('trial vectors',print_debug)) then
            do i=1,ntrials
              call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1          g_prod_p,g_prod_m,g_buffer1,d_trials,d_prod_p,d_prod_m,
     2          i,ntrials,nproducts,tda,algorithm)
            enddo
          endif
c
c Davidson iteration
c
          if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1      call errquit('tddft_energy: failed to set triplet',0,
     &       RTDB_ERR)
          call tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1      nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2      g_buffer1,g_buffer2,d_trials,d_prod_p,d_prod_m,
     3      g_movecs,dbl_mb(k_evl),ntrials,nproducts,
     4      tda,oskel,.true.,.false.,tol2e,xfac(1),lhashf,thresh,
     5      maxiter,target,targetsym,algorithm,mult,l_irs,k_irs)
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
          call tddft_iguess(rtdb,geom,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1      nroots,maxvecs,g_trials,g_buffer1,g_buffer2,d_trials,
     2      g_movecs,dbl_mb(k_evl),ntrials,nproducts,algorithm,
     3      targetsym,symmetry,l_irs,k_irs,tda)
          if (util_print('trial vectors',print_debug)) then
            do i=1,ntrials
              call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1          g_prod_p,g_prod_m,g_buffer1,d_trials,d_prod_p,d_prod_m,
     2          i,ntrials,nproducts,tda,algorithm)
            enddo
          endif
c
c Davidson iteration
c
          if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.true.))
     1      call errquit('tddft_energy: failed to set triplet',0,
     &       RTDB_ERR)
          call tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1      nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2      g_buffer1,g_buffer2,d_trials,d_prod_p,d_prod_m,
     3      g_movecs,dbl_mb(k_evl),ntrials,nproducts,
     4      tda,oskel,.false.,.true.,tol2e,xfac(1),lhashf,thresh,
     5      maxiter,target,targetsym,algorithm,mult,l_irs,k_irs)
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
        call tddft_iguess(rtdb,geom,nbf_ao,ipol,nocc,nmo,nfc,nfv,
     1    nroots,maxvecs,g_trials,g_buffer1,g_buffer2,d_trials,
     2    g_movecs,dbl_mb(k_evl),ntrials,nproducts,algorithm,
     3    targetsym,symmetry,l_irs,k_irs,tda)
        if (util_print('trial vectors',print_debug)) then
          do i=1,ntrials
            call tddft_dump(ipol,nocc,nmo,nfc,nfv,g_trials,
     1        g_prod_p,g_prod_m,g_buffer1,d_trials,d_prod_p,d_prod_m,
     2        i,ntrials,nproducts,tda,algorithm)
          enddo
        endif
c
c Davidson iteration
c
        call tddft_davidson(rtdb,geom,ao_bas_han,nbf_ao,ipol,
     1    nocc,nmo,nfc,nfv,nroots,maxvecs,g_trials,g_prod_p,g_prod_m,
     2    g_buffer1,g_buffer2,d_trials,d_prod_p,d_prod_m,
     3    g_movecs,dbl_mb(k_evl),ntrials,nproducts,
     4    tda,oskel,.false.,.false.,tol2e,xfac(1),lhashf,thresh,
     5    maxiter,target,targetsym,algorithm,mult,l_irs,k_irs)
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
     2    call errquit('tddft_energy: failed to restore xc_active',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,
     1    calc_type_save))
     2    call errquit('tddft_energy: failed to restore calc_type',0,
     &       RTDB_ERR)
      endif
      if(.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1  call errquit('tddft_energy: failed to unset triplet',0,
     &       RTDB_ERR)
      if (rtdb_get(rtdb,'tddft:energy',mt_dbl,1,energy))
     1  tddft_energy=.true.
      call tddft_finalize(geom,ao_bas_han,ipol,g_movecs,l_evl,k_evl,
     1  l_occ,k_occ,l_irs,k_irs)
      if ((algorithm.eq.1).or.(algorithm.eq.2)) then
        do i=1,ipol
          if (.not.ga_destroy(g_trials(i))) call errquit
     1      ('tddft_energy: failed to destroy g_trials',0, GA_ERR)
          if (.not.ga_destroy(g_prod_p(i))) call errquit
     1      ('tddft_energy: failed to destroy g_prod_p',0, GA_ERR)
          if (.not.tda) then
            if (.not.ga_destroy(g_prod_m(i))) call errquit
     1        ('tddft_energy: failed to destroy g_prod_q',0, GA_ERR)
          endif
        enddo
      else
        do i=1,ipol
          if (.not.ga_destroy(g_buffer1(i))) call errquit
     1      ('tddft_energy: failed to destroy g_buffer1',0, GA_ERR)
          if (.not.ga_destroy(g_buffer2(i))) call errquit
     1      ('tddft_energy: failed to destroy g_buffer2',0, GA_ERR)
          dummy=dra_delete(d_trials(i))
          dummy=dra_delete(d_prod_p(i))
          if (.not.tda) then
            dummy=dra_delete(d_prod_m(i))
          endif
        enddo
        dummy=dra_terminate()
      endif
      return
      end
