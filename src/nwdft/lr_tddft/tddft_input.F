      subroutine tddft_input(rtdb)
c
c $Id: tddft_input.F,v 1.9 2003-10-17 22:58:29 carlfahl Exp $
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
c        [MAXVECS <integer maxvecs default 1000>]
c        [(SINGLET||NOSINGLET) default SINGLET]
c        [(TRIPLET||NOTRIPLET) default TRIPLET]
c        [THRESH <double thresh default 1e-4>]
c        [MAXITER <integer maxiter default 100>]
c        [TARGET <integer target default 1>]
c        [TARGETSYM <character targetsym default 'none'>]
c        [SYMMETRY]
c        [ALGORITHM <integer algorithm default 1>]
c        [VECTOR <character vector default jobname'.tddft']
c        [PRINT (none||low||medium||high||debug)
c          <string list_of_names ...>]
c
c ... jochen:
c     added keyword for CD spectra
c     [CDSPECTRUM <not default: compute rotatory strengths>]
c     [VELOCITY]  <not default: compute velocity rep. of trans. dipoles>]
c     /jochen
c
c     END
c     
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*20 test
      integer nroots,maxvecs,maxiter,target,algorithm
      character*4 targetsym
C     character*80 vector
      integer ihuge
      parameter (ihuge=1000000000)
      logical singlet,triplet
      logical tda,symmetry, cdspectrum, velocity
      double precision thresh
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
c     VECTOR (filename for solution vectors)
c
C     call util_file_name('tddft',.false.,.false.,vector)
C     if (.not.rtdb_cput(rtdb,'tddft:vector',1,vector))
C    1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c
c     CIS (Tamm-Dancoff approximation)
c
c ... jochen:
      cdspectrum  = .false.
      if (.not.rtdb_put(rtdb,'tddft:cdspectrum',mt_log,1,cdspectrum))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
      velocity  = .false.
      if (.not.rtdb_put(rtdb,'tddft:velocity',mt_log,1,velocity))
     1  call errquit('tddft_input: failed writing to rtdb',0, RTDB_ERR)
c     /jochen
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
c ... jochen:
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
c     /jochen
c
c     VECTOR (filename for solution vectors)
c
C     else if (inp_compare(.false.,test,'vector')) then
C       if (.not.inp_a(vector)) then
C         write(LuOut,*) 'tddft_input: filename not found; ',
C    1      'default name will be used'
C         call util_file_name('tddft',.false.,.false.,vector)
C       endif
C       if (.not.rtdb_cput(rtdb,'tddft:vector',1,vector))
C    1    call errquit('tddft_input: failed writing to rtdb',0)
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
