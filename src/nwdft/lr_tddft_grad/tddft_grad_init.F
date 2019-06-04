c
c     Initialize TDDFT gradient calculation
c     HvD 4/2006, NG 11/2012
c
      subroutine tddft_grad_init(rtdb,ipol,tda,nat,
     +   nao,nocc,nmo,nfc,nfv,naoc,nav,nov,fn_movecs,numfunc,
     +   xname,cname,xcname,xfac,cfac,lxfac,nlxfac,lcfac,
     +   nlcfac,xccomb,lhf,ldft,lhybrid,lhasdft,lhashf,
     +   nroots, iptr_roots, ihdl_roots, fn_civecs,
     +   nsroots,iptr_sroots,ihdl_sroots,fn_civecs_singlet,
     +   ntroots,iptr_troots,ihdl_troots,fn_civecs_triplet,
     +   osinglet,otriplet,ntotroots,ihdl_bfao,ihdl_geom,tol2e,
     +   solve_thresh,solve_maxit,solve_maxsub,
     +   fn_solve,solve_when,kfac,title)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
#include "global.fh"
#include "case.fh"
c
      integer rtdb        ! Runtime database handle
c
      integer ipol        ! =1 (restricted), =2 (unrestricted)
      logical tda         ! True if Tamm-Dancoff approximation
      integer nat         ! the number of atoms
      integer nsroots     ! the number of singlet roots for which the 
c                         ! gradient should be calculated
      integer iptr_sroots ! a list of the actual roots
      integer ihdl_sroots ! the memory handle of the above
      integer ntroots     ! the number of triplet roots for which the 
c                         ! gradient should be calculated
      integer iptr_troots ! a list of the actual roots
      integer ihdl_troots ! the memory handle of the above
      integer nroots      ! the number of unrestricted roots for which the 
c                         ! gradient should be calculated
      integer iptr_roots  ! a list of the actual roots
      integer ihdl_roots  ! the memory handle of the above
      integer ihdl_bfao   ! the AO basis set handle
      integer ihdl_geom   ! the geometry handle
c
      integer nrootscalc  ! number of calculated roots in TDDFT
c
      integer len_fn_civecs           ! length of civecs file name
      character*(*) fn_civecs         ! file name for TDDFT vectors
      character*(*) fn_civecs_singlet ! file name for TDDFT singlet vectors
      character*(*) fn_civecs_triplet ! file name for TDDFT triplet vectors
c
      integer n,numfunc
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
      double precision xca,shift
c
      logical lhf                       ! True if HF (no DFT)
      logical ldft                      ! True if pure DFT (no HF exch)
      logical lhybrid                   ! True if hybrid DFT and HF
      logical lhasdft                   ! True if it has DFT
      logical lhashf                    ! True if it has HF
      logical use_nwxc                  ! True if using NWXC library
c
      logical cs00
      logical lb94
      logical ncap
      integer iAOacc
c
      logical osinglet
      logical otriplet
      integer ntotroots
c
      integer nao         ! the number of AOs
      integer nmo(2)      ! the number of MOs
      integer nfc(2)      ! the number of frozen core orbitals
      integer nfv(2)      ! the number of frozen virtual orbitals
      integer nocc(2)     ! the number of occupied orbitals
      integer naoc(2)     ! the number of active occupied orbitals
      integer nav(2)      ! the number of active virtual orbitals
      integer nov(2)      ! the dimension of the various arrays
c
      double precision tol2e ! the 2-electron integral tolerance
c
      double precision solve_thresh ! convergence threshold for solver
c
      integer solve_maxit  ! the maximum number of iterations in solver
      integer solve_maxsub ! the maximum number expansion vectors per root
c
      character*(*) fn_solve   ! file name to store linear sys solutions
      character*(*) solve_when ! when are the solutions to be stored
c
      double precision kfac ! the weight of the Hartree-Fock exchange
                            ! contributions
      character*(*) title      ! Title
c
c     Local variables
c
      double precision eps
      parameter (eps = 1.0e-8)
      character*(*) fn_movecs  ! MO vector file name
      character*255 basisname  ! Basis set name
      character*255 scftype    ! Type of SCF
      integer ipol_verify      ! check on ipol
      integer nao_verify       ! check on nao
      integer i                ! counters
      character*18 ctype
      integer itol2e           ! log10 of tol2e
c
      character*32 pname
c
c     Functions
c
      logical  int_normalize, movecs_read_header
      external int_normalize, movecs_read_header
c
      logical xc_ngot3rd
      external xc_ngot3rd
c
c     preliminaries: general parameters
      nao     = 0
      nmo(1)  = 0
      nmo(2)  = 0
      nfc(1)  = 0
      nfc(2)  = 0
      nfv(1)  = 0
      nfv(2)  = 0
      nocc(1) = 0
      nocc(2) = 0
      naoc(1) = 0
      naoc(2) = 0
      nav(1)  = 0
      nav(2)  = 0
      nov(1)  = 0
      nov(2)  = 0
      pname = "tddft_grad_init: "
c
c     get the title 
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     get the geometry
      if (.not.geom_create(ihdl_geom,'geometry'))
     +  call errquit(pname//'geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,ihdl_geom,'geometry'))
     +  call errquit(pname//'geom_rtdb_load failed',0, RTDB_ERR)
      if (.not.geom_ncent(ihdl_geom,nat)) 
     +  call errquit(pname//'geom_ncent failed',0,GEOM_ERR)
c
c     get the basis set
      if (.not.bas_create(ihdl_bfao,'ao basis'))
     +  call errquit(pname//'bas_create failed',0, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,ihdl_geom,ihdl_bfao,'ao basis'))
     +  call errquit(pname//'bas_rtdb_load failed',0, BASIS_ERR)
      if (.not.int_normalize(rtdb,ihdl_bfao))
     +  call errquit(pname//'int_normalize failed',0, BASIS_ERR)
      if (.not.bas_numbf(ihdl_bfao,nao))
     +  call errquit(pname//'bas_numbf failed',0, BASIS_ERR)
c
c     get XC functionals
      if (.not.rtdb_get(rtdb,'dft:xfac',mt_dbl,numfunc,xfac))
     1  call errquit(pname//'failed to read xfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:cfac',mt_dbl,numfunc,cfac))
     1  call errquit(pname//'failed to read cfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:lxfac',mt_log,numfunc,lxfac))
     1  call errquit(pname//'failed to read lxfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:nlxfac',mt_log,numfunc,nlxfac))
     1  call errquit(pname//'failed to read nlxfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:lcfac',mt_log,numfunc,lcfac))
     1  call errquit(pname//'failed to read lcfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:nlcfac',mt_log,numfunc,nlcfac))
     1  call errquit(pname//'failed to read nlcfac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:xccomb',mt_log,numfunc,xccomb))
     1  call errquit(pname//'failed to read xccomb',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:cname',numfunc,cname))
     1  call errquit(pname//'failed to read xname',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:xcname',numfunc,xcname))
     1  call errquit(pname//'failed to read cname',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:xname',numfunc,xname))
     1  call errquit(pname//'failed to read xcname',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:lb94',mt_log,1,lb94))
     1  call errquit(pname//'failed to read lb94',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:cs00',mt_log,1,cs00))
     1  call errquit(pname//'failed to read cs00',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:ncap',mt_log,1,ncap))
     1  call errquit(pname//'failed to read ncap',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:delta_ac',mt_dbl,1,shift))
     1  call errquit(pname//'failed to read delta_ac',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:cam_exch',mt_log,1,cam_exch))
     1   cam_exch=.false.
c
c     get calculation settings: ipol, nocc, nmo, tol2e
      if (.not.rtdb_get(rtdb,'dft:ipol',mt_int,1,ipol))
     +  call errquit(pname//'failed to read ipol',0, RTDB_ERR)
      if (.not.rtdb_get(rtdb,'dft:noc',mt_int,2,nocc))
     +  call errquit(pname//'failed to read nocc',0, RTDB_ERR)
      if (.not.rtdb_cget(rtdb,'dft:output vectors',1,fn_movecs))
     +  call errquit(pname//' no final MO found',0, RTDB_ERR)
      if (.not.movecs_read_header(fn_movecs,title,basisname,
     +  scftype,nao_verify,ipol_verify,nmo,2)) 
     +  call errquit(pname//'failed to read MO file header',0, DISK_ERR)
      if (nao.ne.nao_verify.or.ipol.ne.ipol_verify) 
     +  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
c
      if (.not.rtdb_get(rtdb,'dft:itol2e',mt_int,1,itol2e))
     +  call errquit(pname//'failed to read itol2e',0,RTDB_ERR)
      tol2e=10.0d0**(-itol2e)
c
c     Default to 14 as in nwdft/input_dft/tol_inp.F
      if (.not.rtdb_get(rtdb,'dft:iAOacc',mt_int,1,iAOacc))iAOacc = 14
c
c     get TDDFT settings: tda, nfc, nfv    
      if (.not.rtdb_get(rtdb,'tddft:tda',mt_log,1,tda)) 
     +  call errquit(pname//'failed to read tda',0,RTDB_ERR)
      if(.not.geom_freeze(rtdb, ihdl_geom, 'tddft', nfc(1))) then
        call errquit(pname//'failed to get number of cores',0,GEOM_ERR)
      endif
      if (.not.rtdb_get(rtdb,'tddft:frozen virtual',mt_int,1,nfv(1))) 
     +    nfv(1) = 0
      if (ipol.eq.2) then
        nfc(2) = nfc(1)
        nfv(2) = nfv(1)
      endif
c
c     get calculation type
      xca=0.0d0
      do n=2,numfunc
        xca=xca+dabs(xfac(n))
      enddo
      do n=1,numfunc
        xca=xca+dabs(cfac(n))
      enddo
      if ((xca.gt.1.0d-8).and.(dabs(xfac(1)).le.1.0d-8)) then
c       Pure DFT
        ldft=.true.
        lhf=.false.
        lhybrid=.false.
        lhasdft=.true.
        lhashf=.false.
      else if ((xca.gt.1.0d-8).and.(dabs(xfac(1)).gt.1.0d-8)) then
c       Hybrid DFT & HF
        ldft=.false.
        lhf=.false.
        lhybrid=.true.
        lhasdft=.true.
        lhashf=.true.
      else if ((xca.le.1.0d-8).and.(dabs(xfac(1)-1.0d0)).le.1.0d-8) then
c       Pure HF
        ldft=.false.
        lhf=.true.
        lhybrid=.false.
        lhasdft=.false.
        lhashf=.true.
      else
c       Screened HF
        ldft=.false.
        lhf=.false.
        lhybrid=.false.
        lhasdft=.false.
        lhashf=.true.
      endif
      if (tda) then
        if (ldft.or.lhybrid) then
          ctype='Tamm-Dancoff TDDFT'
        else if (lhf) then
          ctype='CIS               '
        else
          ctype='Screened CIS      '
        endif
      else
        if (ldft.or.lhybrid) then
          ctype='TDDFT             '
        else if (lhf) then
          ctype='TDHF              '
        else
          ctype='Screened TDHF     '
        endif
      endif
c
      use_nwxc = util_module_avail("nwxc")
      if (use_nwxc) then
        call nwxc_rtdb_load(rtdb,"dft",use_nwxc)
      endif
      if (use_nwxc) then
         call nwxc_getvals("nwxc_has_hfx",lhashf)
        if (lhashf) then
          call nwxc_getwght("nwxc_wght_hfx",xfac(1))
        endif
        call  nwxc_getvals("nwxc_has_xc",lhasdft)
        ldft    = lhasdft.and.(.not.lhashf)
        lhf     = (.not.lhasdft).and.lhashf
        lhybrid = lhasdft.and.lhashf
c       call nwxc_getvals("nwxc_has_mp2c",dftmp2)
c       if (dftmp2) then
c          call nwxc_getwght("nwxc_wght_mp2c",mp2fac)
c         if (.not.rtdb_put(rtdb,'dft:mp2fac', mt_dbl, 1, mp2fac))
c    &      call errquit('dft_inpana: rtdb_put failed', 2902, RTDB_ERR)
c       endif
        call nwxc_getvals("nwxc_has_cam",cam_exch)
        if (cam_exch) then
          call nwxc_get_cam(cam_alpha,cam_beta,cam_omega,cam_srhf)
        else
          cam_alpha = 0.0d0
          cam_beta  = 0.0d0
          cam_omega = 0.0d0
        endif
        call nwxc_print()
c       if (cam_exch.and.(ga_nodeid().eq.0)) then
c           write(LuOut,*)
c           write(LuOut,8202)
c    &         'Range-Separation Parameters        '
c           write(LuOut,8203)
c           write(LuOut,8201)cam_alpha,cam_beta,cam_omega,cam_srhf
c       end if ! cam_exch
      endif
c
      kfac   = xfac(1)
      lhashf = kfac.ne.0.0d0
c
c     TDDFT gradient settings
      nsroots = 0
      ihdl_sroots = 0
      iptr_sroots = 0
c
      ntroots = 0
      ihdl_troots = 0
      iptr_troots = 0
c
      nroots = 0
      ihdl_roots = 0
      iptr_roots = 0
c
c     are singlets and triplets calculated as part of the TDDFT energy ?
      if (.not.rtdb_get(rtdb,'tddft:singlet',mt_log,1,osinglet))
     &            osinglet=.true.
      if (.not.rtdb_get(rtdb,'tddft:triplet',mt_log,1,otriplet))
     &            otriplet=.true.
c
c     calculated roots in the main TDDFT calculation
      nrootscalc = 0
      if (.not.rtdb_get(rtdb,'tddft:nroots',mt_int,1,nrootscalc))
     +  call errquit(pname//'no excited states calculated',nrootscalc,
     +                  INPUT_ERR)
c
c     gradients can be calculated for only one state at a time
      if (ipol.eq.1) then  ! restricted
c
        if (osinglet) then
           nsroots = 1  ! only one state at a time
           if (.not.ma_push_get(mt_int,nsroots,'singlet roots',
     +                         ihdl_sroots,iptr_sroots))
     +     call errquit(pname//'failed to allocate iptr_sroots',nsroots,
     +                  MA_ERR)
           if (.not.rtdb_get(rtdb,'tddft_grad:isinglet_roots',mt_int,
     +                      nsroots,int_mb(iptr_sroots)))
     +     call errquit(pname//'failed to read isroots',0,RTDB_ERR)
        endif ! singlet
c
        if (otriplet) then
           ntroots = 1  ! only one state at a time
           if (.not.ma_push_get(mt_int,ntroots,'triplet roots',
     +                         ihdl_troots,iptr_troots))
     +     call errquit(pname//'failed to allocate troots',ntroots,
     +               MA_ERR)
           if (.not.rtdb_get(rtdb,'tddft_grad:itriplet_roots',mt_int,
     +                      ntroots,int_mb(iptr_troots)))
     +     call errquit(pname//'failed to read itroots',0,RTDB_ERR)
        endif ! triplet
c
      else   ! unrestricted
          nroots = 1  ! only one state at a time
          if (.not.ma_push_get(mt_int,nroots,'roots',
     +                         ihdl_roots,iptr_roots))
     +    call errquit(pname//'failed to allocate roots',nroots,MA_ERR)
          if (.not.rtdb_get(rtdb,'tddft_grad:iroots',mt_int,
     +                      nroots,int_mb(iptr_roots))) 
     +    call errquit(pname//'failed to read iroots',0,RTDB_ERR)
c
      endif  ! ipol
c
c     CI Vectors files for restricted singlets, triplets and unrestricted
      if (rtdb_cget(rtdb,'tddft:civecs',1,fn_civecs)) then
        len_fn_civecs = inp_strlen(fn_civecs)
        if (nsroots.gt.0) then
          fn_civecs_singlet = fn_civecs(1:len_fn_civecs)//"_singlet"
          if (ga_nodeid().eq.0) write(luout,*) fn_civecs_singlet
        else if (ntroots.gt.0) then
          fn_civecs_triplet = fn_civecs(1:len_fn_civecs)//"_triplet"
          if (ga_nodeid().eq.0) write(luout,*) fn_civecs_triplet
        endif
      else
        call errquit(pname//'no saved tddft results?',0,CALC_ERR)
      endif
c
c     Get the settings for the linear system solver
      if (.not.rtdb_get(rtdb,'tddft_grad:solve_thresh',mt_dbl,1,
     +                  solve_thresh)) solve_thresh=1.0d-4
c
c Daniel (3-29-13): For a large system, the original default value (10) 
c of solve_maxsub might not work efficiently.  I'd suggest using a 
c larger value for the default to aid in convergence.
c
      if (.not.rtdb_get(rtdb,'tddft_grad:solve_maxsub',mt_int,1,
     +                  solve_maxsub)) solve_maxsub=40
      if (.not.rtdb_get(rtdb,'tddft_grad:solve_maxit',mt_int,1,
     +                  solve_maxit)) solve_maxit=250
c
      if (.not.rtdb_cget(rtdb,'tddft:solve_file',1,fn_solve))
     +                  fn_solve = ''
      if (.not.rtdb_cget(rtdb,'tddft:solve_when',1,solve_when))
     +                  solve_when = 'never'

c     Modify the runtime database to suppress subsequent solutions
c     from being stored.
      if (solve_when.eq.'once') then
        if (.not.rtdb_cget(rtdb,'tddft:solve_when',1,'never')) 
     +     call errquit(pname//'cannot update solve_when',0, RTDB_ERR)
      endif
c
      do i = 1,ipol
        naoc(i)= (nocc(i)-nfc(i))
        nav(i) = (nmo(i)-nocc(i)-nfv(i))
        nov(i) = naoc(i)*nav(i)
      enddo
c
c     Checks
      ntotroots = nsroots+ntroots+nroots
      if (ntotroots.gt.nov(1)+nov(2)) ntotroots=nov(1)+nov(2)
c
      otriplet = .false.
      if (nroots.eq.0) then
        if (nsroots.eq.0.and.ntroots.gt.0) then
          otriplet = .true.
        elseif (nsroots.gt.0.and.ntroots.eq.0) then
          otriplet = .false.
        elseif (nsroots.eq.0.and.ntroots.eq.0) then
          call errquit(pname//'no excited states',0,INPUT_ERR)
        else
          call errquit(pname//'one spin state at a time',0,INPUT_ERR)
        endif
      else
       if (nsroots.gt.0.or.ntroots.gt.0) 
     +  call errquit(pname//'invalid combination of states',0,INPUT_ERR)
      endif
c
c     Check for 3rd derivatives
      if (xc_ngot3rd()) then
        call errquit(pname//'XC 3rd derivatives not implemented', 
     1    0, CAPMIS_ERR)
      endif
c
      end
c $Id$
