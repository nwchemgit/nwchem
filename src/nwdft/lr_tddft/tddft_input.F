      subroutine tddft_input(rtdb)
c
c $Id$
c
c Input parser for TDDFT called by input/input_parse.F.
c Also sets default values for input parameters.
c Modified from mp2_grad/mp2_input.F by So Hirata Feb, 2002.
c
c     TDDFT
c        [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c                 [virtual <integer nfzv default 0>]]
c        [(CIS||RPA) default RPA]
c        [NROOTS <integer nroots default 1>]
c        [ECUT   <not default: double ecut default -ihuge>]
c        [EWIN]
c        [ALPHA <integer astart> <integer aend>]
c        [BETA <integer bstart> <integer bend>]
c        [MAXVECS <integer maxvecs default 1000>]
c        [(SINGLET||NOSINGLET) default SINGLET]
c        [(TRIPLET||NOTRIPLET) default TRIPLET]
c        [THRESH <double thresh default 1e-4>]
c        [MAXITER <integer maxiter default 100>]
c        [TARGET <integer target default 1>]
c        [TARGETSYM <character targetsym default 'none'>]
c        [SYMMETRY]
c        [ALGORITHM <integer algorithm default 1>]
c        [CIVECS <default '*.civecs']
c        [TRANSDEN <default '*.tdens']
c        [PRINT (none||low||medium||high||debug)
c          <string list_of_names ...>]
c        [CDSPECTRUM <not default: compute rotatory strengths>]
c        [GIAO <not default: calc. magnetic trans. dipoles with GIAOs>]
c        [VELOCITY]  <not default: compute velocity gaude dipoles>]
c     END
c     
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*20 test
      integer nroots,maxvecs,maxiter,target,algorithm
      character*4 targetsym
      character*(nw_max_path_len) fn_civecs
      character*(nw_max_path_len) fn_transden
      integer ihuge
      parameter (ihuge=1000000000)
      logical singlet,triplet
      logical tda,symmetry, cdspectrum, velocity, lgiao, gmh
      double precision thresh
      double precision ecut
      integer vshift
      integer occtop(2)
      logical lecut
      double precision ewinl,ewinh
      logical lewin
      logical lowin
      integer owstart(2)
      integer owend(2)
      double precision diff_max
      logical lcivecs
      logical ltransden
c
c ------------------
c Set default values
c ------------------
c
c     FREEZE (Frozen cores/virtuals)
c
c     no action is taken
c
c     CIS (Tamm-Dancoff approximation)
c
      tda=.false.
      if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     NROOTS (the number of excited states sought) 
c
      nroots=1
      if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     ECUT (energy cutoff)
c
      ecut  = 0.d0  ! default
      lecut = .false. ! default
      if (.not.rtdb_put(rtdb,'tddft:ecut',mt_dbl,1,ecut))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tddft:lecut',mt_log,1,lecut))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     VSHIFT (reference virtual state)
c
      vshift = 0
      if (.not.rtdb_put(rtdb,'tddft:vshift',mt_int,1,vshift))
     &  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     OCCSTART 
c
      occtop(1) = 0
      occtop(2) = 0
      if (.not.rtdb_put(rtdb,'tddft:occtop',mt_int,2,occtop))
     &  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     OWIN
c
      lowin = .false.
      owend(1) = 0
      owend(2) = 0
      owstart(1) = 0
      owstart(2) = 0
      if (.not.rtdb_put(rtdb,'tddft:lecut',mt_log,1,lowin))
     1  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tddft:owstart',mt_int,2,owstart))
     &  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tddft:owend',mt_int,2,owend))
     &  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     EWIN
c
      lewin = .false.
      ewinh = 0.d0
      ewinl = 0.d0
      if (.not.rtdb_put(rtdb,'tddft:lecut',mt_log,1,lewin))
     1  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tddft:ewinh',mt_dbl,1,ewinh))
     &  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tddft:ewinl',mt_dbl,1,ewinl))
     &  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     DIFF_MAX (largest eigenvalue difference) : For solver
c
      diff_max = 100000.d0 ! default
      if (.not.rtdb_put(rtdb,'tddft:diff_max',mt_dbl,1,diff_max))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     MAXVECS (the maximum number of trial vectors)
c
      maxvecs=ihuge
      if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     SINGLET (singlet excited states sought for a closed shell)
c
      singlet=.true.
      if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     TRIPLET (triplet excited states sought for a closed shell)
c
      triplet=.true.
      if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     THRESH (convergence threshold for Davidson iteration)
c
      thresh=1.0d-4
      if (.not.rtdb_put(rtdb,'tddft:thresh',mt_dbl,1,thresh))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      maxiter=100
      if (.not.rtdb_put(rtdb,'tddft:maxiter',mt_int,1,maxiter))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     TARGET (the target excited state for, e.g., geometry optimization)
c
      target=1
      if (.not.rtdb_put(rtdb,'tddft:target',mt_int,1,target))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     TARGETSYM (the irrep of the target excited state)
c
      targetsym='none'
      if (.not.rtdb_cput(rtdb,'tddft:targetsym',1,targetsym))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
c
      symmetry=.false.
      if (.not.rtdb_put(rtdb,'tddft:symmetry',mt_log,1,symmetry))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     ALGORITHM (algorithm: =1 incore-multiple fock build; =2: incore-single fock build
c                           =3 disk-multiple fock build;   =4: disk-single fock build)
c
      algorithm=0
      if (.not.rtdb_put(rtdb,'tddft:algorithm',mt_int,1,algorithm))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     CIVECS (filename for solution vectors)
c
      call util_file_name('civecs',.false.,.false.,fn_civecs)
      if (.not.rtdb_cput(rtdb,'tddft:civecs',1,fn_civecs))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     TRANSDEN
c
      call util_file_name('tdens',.false.,.false.,fn_transden)
      if (.not.rtdb_cput(rtdb,'tddft:transden',1,fn_transden))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     CDSPECTRUM
c
      cdspectrum  = .false.
      if (.not.rtdb_put(rtdb,'tddft:cdspectrum',mt_log,1,cdspectrum))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
      velocity  = .false.
      if (.not.rtdb_put(rtdb,'tddft:velocity',mt_log,1,velocity))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
      lgiao = .false.
      if (.not.rtdb_put(rtdb,'tddft:giao',mt_log,1,lgiao))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)

c
c     GENERALIZED MULLIKEN-HUSH H_IJ
c
      gmh = .false.
      if (.not.rtdb_put(rtdb,'tddft:gmh',mt_log,1,gmh))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c ----------
c Read input
c ----------
c
 10   if (.not. inp_read()) 
     1  call errquit('tddft_input: failed reading input',0, INPUT_ERR)
      if (.not. inp_a(test)) 
     1  call errquit('tddft_input: failed reading keyword',0, INPUT_ERR)
c
c     FREEZE (Frozen cores/virtuals)
c
      if (inp_compare(.false.,test,'freeze')) then
        call freeze_input(rtdb,'tddft')
c
c     CIS (Tamm-Dancoff approximation)
c
      else if (inp_compare(.false.,test,'cis')) then
        tda=.true.
        if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     TDA (Tamm-Dancoff approximation)
c
      else if (inp_compare(.false.,test,'tda')) then
        tda=.true.
        if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     RPA (Full linear response)
c
      else if (inp_compare(.false.,test,'rpa')) then
        tda=.false.
        if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     NROOTS (the number of excited states sought)
c
      else if (inp_compare(.false.,test,'nroots')) then
        if (.not.inp_i(nroots)) then
          write(LuOut,*) 'tddft_input: nroots value not found; ',
     1      'default value of 1 will be used'
          nroots=1
        endif
        if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     ECUT (energy cutoff for occupied states)
c
      else if (inp_compare(.false.,test,'ecut')) then
        lecut = .true.
        lewin = .false.
        lowin = .false.
        if (lowin .or. lewin) then
          write(LuOut,*) 'tddft_input: cannot use energy cutoff and ',
     &      'orbital window or energy window simultaneously.'
          call errquit('tddft_input: orbital window error',0,INPUT_ERR)
        endif
        if (.not.inp_f(ecut)) then
          write(LuOut,*) 'tddft_input: ecut value not found; ',
     1      'default value of 0.d0 will be used'
          ecut  = 0.d0  ! default
        endif
        if (.not.rtdb_put(rtdb,'tddft:ecut',mt_dbl,1,ecut))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:lecut',mt_log,1,lecut))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     VSHIFT (reference virtual state)
c
      else if (inp_compare(.false.,test,'vshift')) then
        if (.not.inp_i(vshift)) vshift = 0
        if (.not.rtdb_put(rtdb,'tddft:vshift',mt_int,1,vshift))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     OCCTOP
c
      else if (inp_compare(.false.,test,'occtop')) then
        if (.not.inp_i(occtop(1))) occtop(1) = 0
        if (.not.inp_i(occtop(2))) occtop(2) = 0
        if (.not.rtdb_put(rtdb,'tddft:occtop',mt_int,2,occtop))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     EWIN (energy window for occupied states)
c
      else if (inp_compare(.false.,test,'ewin')) then
        lewin = .true.
        lecut = .false.
        lowin = .false.
        if (lecut .or. lowin) then
          write(LuOut,*) 'tddft_input: cannot use energy cutoff and ',
     &      'orbital or energy cutoff simultaneously.'
          call errquit('tddft_input: orbital window error',0,INPUT_ERR)
        endif
        if (.not.inp_f(ewinl)) lewin = .false.
        if (.not.inp_f(ewinh)) lewin = .false.
        if (.not.rtdb_put(rtdb,'tddft:ewinl',mt_dbl,1,ewinl))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:ewinh',mt_dbl,1,ewinh))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:lewin',mt_log,1,lewin))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c      DIFF_MAX (maximum eigenvalue difference): For Solver
c
      else if (inp_compare(.false.,test,'diff_max')) then
        if (.not.inp_f(diff_max)) then
          write(LuOut,*) 'tddft_input: diff_max value not found; ',
     1      'default value of 100000.0 will be used'
          diff_max = 100000.d0
        endif
        if (.not.rtdb_put(rtdb,'tddft:diff_max',mt_dbl,1,diff_max))
     &    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     ALPHA: alpha orbital window
c
      else if (inp_compare(.false.,test,'alpha')) then
        lowin = .true.
        lewin = .false.
        lecut = .false.
        if (lecut .or. lewin) then
          write(LuOut,*) 'tddft_input: cannot use energy cutoff and ',
     &      'orbital or energy window simultaneously.'
          call errquit('tddft_input: orbital window error',0,INPUT_ERR)
        endif
        if (.not.rtdb_get(rtdb,'tddft:owstart',mt_int,2,owstart))
     &    call errquit('tddft_input: failed reading rtdb',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'tddft:owend',mt_int,2,owend))
     &    call errquit('tddft_input: failed reading rtdb',0,RTDB_ERR)
c
c       Set defaults
        if (.not.inp_i(owstart(1))) owstart(1) = 1
        if (.not.inp_i(owend(1))) owend(1) = 1
c
        if (.not.rtdb_put(rtdb,'tddft:lowin',mt_log,1,lowin))
     &    call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:owstart',mt_int,2,owstart))
     &    call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:owend',mt_int,2,owend))
     &    call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     BETA: beta orbital window
c
      else if (inp_compare(.false.,test,'beta')) then
        lowin = .true.
        lewin = .false.
        lecut = .false.
        if (lecut .or. lewin) then
          write(LuOut,*) 'tddft_input: cannot use energy cutoff and ',
     &      'orbital or energy window simultaneously.'
          call errquit('tddft_input: orbital window error',0,INPUT_ERR)
        endif
        if (.not.rtdb_get(rtdb,'tddft:owstart',mt_int,2,owstart))
     &    call errquit('tddft_input: failed reading rtdb',0,RTDB_ERR)
        if (.not.rtdb_get(rtdb,'tddft:owend',mt_int,2,owend))
     &    call errquit('tddft_input: failed reading rtdb',0,RTDB_ERR)
c
c       Set defaults
        if (.not.inp_i(owstart(2))) owstart(2) = 1
        if (.not.inp_i(owend(2))) owend(2) = 1
c
        if (.not.rtdb_put(rtdb,'tddft:lowin',mt_log,1,lowin))
     &    call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:owstart',mt_int,2,owstart))
     &    call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tddft:owend',mt_int,2,owend))
     &    call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     MAXVECS (the maximum number of trial vectors)
c
      else if (inp_compare(.false.,test,'maxvecs')) then
        if (.not.inp_i(maxvecs)) then
          write(LuOut,*) 'tddft_input: maxvecs value not found; ',
     1      'will be determined based on available memory size'
          maxvecs=ihuge
        endif
        if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     NOSINGLET (no singlet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'nosinglet')) then
        singlet=.false.
        if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     SINGLET (singlet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'singlet')) then
        singlet=.true.
        if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     NOTRIPLET (no triplet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'notriplet')) then
        triplet=.false.
        if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     TRIPLET (triplet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'triplet')) then
        triplet=.true.
        if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     THRESH (convergence threshold for Davidson iteration)
c
      else if (inp_compare(.false.,test,'thresh')) then
        if (.not.inp_f(thresh)) then
          write(LuOut,*) 'tddft_input: thresh value not found; ',
     1      'default value of 1e-4 will be used'
          thresh=1.0d-4
        endif
        if (.not.rtdb_put(rtdb,'tddft:thresh',mt_dbl,1,thresh))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      else if (inp_compare(.false.,test,'maxiter')) then
        if (.not.inp_i(maxiter)) then
          write(LuOut,*) 'tddft_input: maxiter value not found; ',
     1      'default value of 100 will be used'
          maxiter=100
        endif
        if (.not.rtdb_put(rtdb,'tddft:maxiter',mt_int,1,maxiter))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     TARGET (the target excited state for, e.g., geometry optimization)
c
      else if (inp_compare(.false.,test,'target')) then
        if (.not.inp_i(target)) then
          write(LuOut,*) 'tddft_input: target value not found; ',
     1      'default value of 1 will be used'
          target=1
        endif
        if (target.gt.nroots) call errquit
     1    ('tddft_input: an illegal value for target',target, INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tddft:target',mt_int,1,target))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     TARGETSYM (the symmetry of the target excited state)
c
      else if (inp_compare(.false.,test,'targetsym')) then
        if (.not.inp_a(targetsym)) then
          write(LuOut,*) 'tddft_input: targetsym value not found; ',
     1      'no symmetry information will be used in specifying target'
          targetsym='none'
        endif
        if (.not.rtdb_cput(rtdb,'tddft:targetsym',1,targetsym))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
c
      else if (inp_compare(.false.,test,'symmetry')) then
        symmetry=.true.
        if (.not.rtdb_put(rtdb,'tddft:symmetry',mt_log,1,symmetry))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     ALGORITHM (algorithm: =1 incore-multiple fock build; =2: incore-single fock build
c                           =3 disk-multiple fock build;   =4: disk-single fock build)
c
      else if (inp_compare(.false.,test,'algorithm')) then
        if (.not.inp_i(algorithm)) then
          write(LuOut,*) 'tddft_input: algorithm value not found; ',
     1      'an optimal choice will be made later'
          algorithm=0
        endif
        if (.not.rtdb_put(rtdb,'tddft:algorithm',mt_int,1,algorithm))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     CDSPECTRUM
c
      else if (inp_compare(.false.,test,'cdspectrum')) then
        cdspectrum=.true.
        if (.not.rtdb_put(rtdb,'tddft:cdspectrum',mt_log,1,cdspectrum))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     VELOCITY
c
      else if (inp_compare(.false.,test,'velocity')) then
        velocity=.true.
        if (.not.rtdb_put(rtdb,'tddft:velocity',mt_log,1,velocity))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     GIAO
c
      else if (inp_compare(.false.,test,'giao')) then
        lgiao=.true.
        if (.not.rtdb_put(rtdb,'tddft:giao',mt_log,1,lgiao))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     GENERALIZED MULLIKEN-HUSH H_IJ COUPLING
c
      else if (inp_compare(.false.,test,'gmh')) then
       gmh=.true.
       if (.not.rtdb_put(rtdb,'tddft:gmh',mt_log,1,gmh))
     1    call errquit('tddft_input: failed writing to rtdb',0,
     &       RTDB_ERR)
c
c     CIVECS (filename for solution vectors)
c
      lcivecs = .false.  ! default   
      else if (inp_compare(.false.,test,'civecs')) then
       lcivecs = .true.
       if (.not.inp_a(fn_civecs)) then
         write(LuOut,*) 'tddft_input: filename not found; ',
     1      'default name will be used'
         call util_file_name('civecs',.false.,.false.,fn_civecs)
       endif
       if (.not.rtdb_cput(rtdb,'tddft:civecs',1,fn_civecs))
     1  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
       if (.not.rtdb_put(rtdb,'tddft:lcivecs',mt_log,1,lcivecs))
     1  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     TRANSDEN (filename for transition density)
c
      ltransden = .false.   
      else if (inp_compare(.false.,test,'transden')) then
       ltransden = .true.
       if (.not.inp_a(fn_transden)) then
         write(LuOut,*) 'tddft_input: filename not found; ',
     1      'default name will be used'
         call util_file_name('tdens',.false.,.false.,fn_transden)
       endif
       if (.not.rtdb_cput(rtdb,'tddft:transden',1,fn_transden))
     1  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
       if (.not.rtdb_put(rtdb,'tddft:ltransden',mt_log,1,ltransden))
     1  call errquit('tddft_input: failed writing to rtdb',0,RTDB_ERR)
c
c     PRINT
c
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'tddft')
c
c     END
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit('tddft_input: unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
