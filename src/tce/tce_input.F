#ifndef __LINE__
#define __LINE__ 0
#endif
      subroutine tce_input(rtdb)
c
c $Id: tce_input.F,v 1.37 2009/03/05 23:48:58 d3p975 Exp $
c
c Input parser for TCE module for various many-electron theories.
c Also sets default values for input parameters.
c Modified from tddft/tddft_input.F by So Hirata Oct, 2002.
c
c     TCE
c        [(DFT||HF||SCF) default HF]
c        [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c                 [virtual <integer nfzv default 0>]]
c        [(LCCD||CCD||CCSD||LCCSD||CCSDT||CCSDTQ|| \ 
c          CCSD(T)||CCSD[T]||QCISD||CISD||CISDT||CISDTQ|| \
c          MBPT2||MBPT3||MBPT4||MP2||MP3||MP4|| \
c          CR-CCSD(T)||CR-CCSD[T]||LR-CCSD(T)||LR-CCSD(TQ)||CCSD(2)_T||CCSD(2)||
c          CCSDT(2)_Q) default CCSD]
c        [THRESH <double thresh default 1e-6>]
c        [MAXITER <integer maxiter default 100>]
c        [PRINT (none||low||medium||high||debug)]
c        [IO (fortran||eaf||ga||sf||replicated||dra||ga_eaf) default ga]
c        [DIIS <integer diis default 5>]
c        [EOMSOL <integer default 1 >]
c        [DIIS2 <integer diis default 5>]
c        [DIIS3 <integer diis default 5>]
c        [NROOTS <integer nroots default 0>]
c        [TARGET <integer target default 1>]
c        [TARGETSYM <character targetsym default 'none'>]
c        [SYMMETRY]
c        [DIPOLE]
c        [TILESIZE <no default (automatically adjusted)>]
c        [FRAGMENT <default -1 (off)>]
c        [(NO)FOCK <logical recompf default .true.>]
c        [ACTIVE_OA <default 0>]
c        [ACTIVE_OB <default 0>]
c        [ACTIVE_VA <default 0>]
c        [ACTIVE_VB <default 0>]
c        [T3A_LVL   <default 0>]
c     END
c
c     TASK TCE ENERGY
c
c     ... or ...
c
c     UCCSDT or UCC or UCCSD(T) etc.
c        [(DFT||HF||SCF) default HF]
c        [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c                 [virtual <integer nfzv default 0>]]
c        [THRESH <double thresh default 1e-6>]
c        [MAXITER <integer maxiter default 100>]
c        [PRINT (none||low||medium||high||debug)]
c        [IO (fortran||c||ga||sf||replicated) default ga]
c        [DIIS <integer diis default 5>]
c        [NROOTS <integer nroots default 0>]
c        [TARGET <integer target default 1>]
c        [TARGETSYM <character targetsym default 'none'>]
c        [SYMMETRY]
c        [DIPOLE]
c        [TILESIZE <no default (automatically adjusted)>]
c        [FRAGMENT <default -1 (off)>]
c        [(NO)FOCK <logical recompf default .true.>]
c        [ACTIVE_OA <default 0>]
c        [ACTIVE_OB <default 0>]
c        [ACTIVE_VA <default 0>]
c        [ACTIVE_VB <default 0>]
c     END
c
c     TASK UCCSDT ENERGY
c
c     ... etc.
c     
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer rtdb
      character*20 test
      integer maxiter
      character*10 model
      character*10 model2e
      character*10 module
      double precision thresh
      double precision maxdiff ! new
      character*10 ioalgchar
      integer ioalg
      integer reference
      integer diis,diis2,diis3
      integer eomsol
c --- level shift --
      double precision zlshift,zlshiftl,zlshift2(2),zlshift3(2)
c ------------------
      integer nroots
      integer target
      integer tilesize
      integer fragment
      character*4 targetsym
      logical symmetry
      logical left
c --- density matrix    
      logical idens
      character*256 file_densmat
c<-d3p975
      integer multipole
      logical recompf
      character*10 perturbative
      character*10 ccsd_var 
      integer oactive(2)
      integer vactive(2)
      integer numact
c --- 4 index transform. ---
      integer maxs,ichopx,i4im,idiskx
c --- EOM solver
      integer hbard
c
c -------------------------------------
c What input block are we dealing with?
c -------------------------------------
c
      if (.not.rtdb_cget(rtdb,'tce:module',1,module))
     1  call errquit('tce_input: failed reading from rtdb',0,
     2  RTDB_ERR)
c
c ------------------
c Set default values
c ------------------
c
c     DFT, HF, or SCF (reference wavefunction)
c
      reference=1
      if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     FREEZE (frozen cores/virtuals)
c
c     no action is taken
c
c     MODEL (the name of CC model requested)
c
c     no action is taken
c
c     THRESH (convergence threshold for Davidson iteration)
c
      thresh=1.0d-7
      if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     LEVEL SHIFT (for singles and doubles)
c
      zlshift=0.0d0
      if (.not.rtdb_put(rtdb,'tce:zlshift',mt_dbl,1,zlshift))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      zlshiftl=0.0d0
      if (.not.rtdb_put(rtdb,'tce:zlshiftl',mt_dbl,1,zlshiftl))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      zlshift2(1)=0.0d0
      zlshift2(2)=0.0d0
      if (.not.rtdb_put(rtdb,'tce:zlshift2',mt_dbl,2,zlshift2))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      zlshift3(1)=0.0d0
      zlshift3(2)=0.0d0
      if (.not.rtdb_put(rtdb,'tce:zlshift3',mt_dbl,2,zlshift3))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      maxiter=100
      if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     IO (I/O method, 0 = Fortran Direct Access, 
c                     1 = C Low-Level I/O,
c                     2 = GA Library,
c                     3 = SF library,
c                     4 = Replicated C Low-Level I/O)
c
      ioalg=2
      if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     DIIS (the vector space size in DIIS)
c
      diis=5
      if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      diis2=5
      if (.not.rtdb_put(rtdb,'tce:diis2',mt_int,1,diis2))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      diis3=5
      if (.not.rtdb_put(rtdb,'tce:diis3',mt_int,1,diis3))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     EOMCC SOLVER
c
      eomsol=1
      if (.not.rtdb_put(rtdb,'tce:eoms',mt_int,1,eomsol))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     DIMENSION OF THE EOM ITERATIVE SPACE
c
      hbard=500
      if (.not.rtdb_put(rtdb,'tce:hbard',mt_int,1,hbard))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     NROOTS (the number of excited state roots)
c
      nroots=0
      if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     2e STORAGE
c
        model2e='2espin'
      if (.not.rtdb_cput(rtdb,'tce:model2e',1,model2e))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c 4ind. transfromation
c
       maxs=30
      if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,maxs))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR) 
c
       ichopx=1
      if (.not.rtdb_put(rtdb,'tce:ichopx',mt_int,1,ichopx))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
       i4im=1
      if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,i4im))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
      idiskx=0
      if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,idiskx))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     TARGET (the target excited state for, e.g., geometry optimization)
c
      target=1
      if (.not.rtdb_put(rtdb,'tce:target',mt_int,1,target))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     TARGETSYM (the irrep of the target excited state)
c
      targetsym='none'
      if (.not.rtdb_cput(rtdb,'tce:targetsym',1,targetsym))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
c
      symmetry=.false.
      if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     IDENS (one particle reduced density matrix)
c
      idens=.false.
      if (.not.rtdb_put(rtdb,'tce:densmat',mt_log,1,idens))
     1   call errquit('tce_input: failed writing to rtdb',0,
     1   rtdb_err)
c
c     DIPOLE (dipole moments & dipole transition moments)
c
      left=.false.
      if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     MULTIPOLE LMAX (multipole moments highest angular momentum)
c
      multipole=0
      if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     FRAGMENT (fragment MO calculations)
c
      fragment=-1
      if (.not.rtdb_put(rtdb,'tce:fragment',mt_int,1,fragment))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     (NO)FOCK (recompute fock for, e.g., DFT, ROHF refs)
c
      recompf=.true.
      if (.not.rtdb_put(rtdb,'tce:recompf',mt_log,1,recompf))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     ACTIVE_OA,OB (Number of active occupied orbitals)
c
      oactive(1)=0
      if (.not.rtdb_put(rtdb,'tce:active_oa',mt_int,1,oactive(1)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      oactive(2)=0
      if (.not.rtdb_put(rtdb,'tce:active_ob',mt_int,1,oactive(2)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     ACTIVE_VA,VB (Number of active virtual orbitals)
c
      vactive(1)=0
      if (.not.rtdb_put(rtdb,'tce:active_va',mt_int,1,vactive(1)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      vactive(2)=0
      if (.not.rtdb_put(rtdb,'tce:active_vb',mt_int,1,vactive(2)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c    ACTIVE EXCITATION LEVEL (number of active orbitals in T3)
c
      numact=0
      if (.not.rtdb_put(rtdb,'tce:act_excit_lvl',mt_int,1,numact))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c ----------
c Read input
c ----------
c
 10   if (.not. inp_read()) 
     1  call errquit('tce_input: failed reading input',0,
     2  RTDB_ERR)
      if (.not. inp_a(test)) 
     1  call errquit('tce_input: failed reading keyword',0,
     2  RTDB_ERR)
c
c     DFT, HF, or SCF (reference wavefunction)
c
      if (inp_compare(.false.,test,'dft')) then
        reference=0
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'hf')) then
        reference=1
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'scf')) then
        reference=1
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     FREEZE (frozen cores/virtuals)
c
      else if (inp_compare(.false.,test,'freeze')) then
        call freeze_input(rtdb,'tce')
c
c     STORAGE OF 2-e INTEGRALS
c
      else if (inp_compare(.false.,test,'2eorb')) then
        if (module.eq.'tce') then
        model2e='2eorb'
        if (.not.rtdb_cput(rtdb,'tce:model2e',1,model2e))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
c
c     MODEL (the name of theory requested)
c
      else if (inp_compare(.false.,test,'multi')) then
        if (module.eq.'tce') then
        model='multi'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
          call errquit('tce_input: multiple theory inputs',0,
     2    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'eionly')) then
        if (module.eq.'tce') then
        model='eionly'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
          call errquit('tce_input: multiple theory inputs',0,
     2    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccd')) then
        if (module.eq.'tce') then
        model='ccd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     2  INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lccd')) then
        if (module.eq.'tce') then
        model='lccd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsd')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lccsd')) then
        if (module.eq.'tce') then
        model='lccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lccsd(t)')) then
        if (module.eq.'tce') then
        model='lccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cr-lccsd(t)')) then
        if (module.eq.'tce') then
        model='lccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='cr_(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsdta')) then
        if (module.eq.'tce') then
        model='ccsdta'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsdt')) then
        if (module.eq.'tce') then
        model='ccsdt'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsdtq')) then
        if (module.eq.'tce') then
        model='ccsdtq'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cc2')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        ccsd_var='cc2'
        if (.not.rtdb_cput(rtdb,'tce:ccsdvar',1,ccsd_var))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lr-ccsd')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        ccsd_var='lr-ccsd'
        if (.not.rtdb_cput(rtdb,'tce:ccsdvar',1,ccsd_var))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsd(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsd[t]')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='[t]'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'qcisd(t)')) then
        if (module.eq.'tce') then
        model='qcisd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cr-qcisd(t)')) then
        if (module.eq.'tce') then
        model='qcisd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='cr_(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
c
c
c
      else if (inp_compare(.false.,test,'lambda-ccsd(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='lambda(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
c      else if (inp_compare(.false.,test,'lambda-ccsd[t]')) then
c        if (module.eq.'tce') then
c        model='ccsd'
c        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
c     1    call errquit('tce_input: failed writing to rtdb',0,
c     2    RTDB_ERR)
c        perturbative='lambda[t]'
c        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
c     1    call errquit('tce_input: failed writing to rtdb',0,
c     2    RTDB_ERR)
c        left=.true.
c        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
c     1    call errquit('tce_input: failed writing to rtdb',0,
c     2    RTDB_ERR)
c        else
c        call errquit('tce_input: multiple theory inputs',0,
c     1    INPUT_ERR)
c        endif
c
c
c
      else if (inp_compare(.false.,test,'cr-ccsd(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='cr_(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lr-ccsd(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='lr_(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'creomsd(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='creom_(t)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'r-creom1(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='emb1'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'r-creom2(t)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='emb2'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lr-ccsd(tq)-1')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='lr_(tq1)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'lr-ccsd(tq)-1p')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='lr_(tq1p)'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cr-ccsd[t]')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='cr_[t]'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsd(2)_t')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='2_t'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsd(2)')) then
        if (module.eq.'tce') then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='2_tq'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'ccsdt(2)_q')) then
        if (module.eq.'tce') then
        model='ccsdt'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='2_q'
        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'qcisd')) then
        if (module.eq.'tce') then
        model='qcisd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cis')) then
        if (module.eq.'tce') then
           model='cis'
           if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1     call errquit('tce_input: failed writing to rtdb',0,
     1     rtdb_err)
        else
           call errquit('tce_input: multiple theory inputs',0,
     1     input_err)
        end if
      else if (inp_compare(.false.,test,'cisd')) then
        if (module.eq.'tce') then
        model='cisd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cisdt')) then
        if (module.eq.'tce') then
        model='cisdt'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'cisdtq')) then
        if (module.eq.'tce') then
        model='cisdtq'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mbpt2')) then
        if (module.eq.'tce') then
        model='mbpt2'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mbpt3')) then
        if (module.eq.'tce') then
        model='mbpt3'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mbpt4')) then
        if (module.eq.'tce') then
        model='mbpt4'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mbpt4(sdq)')) then
        if (module.eq.'tce') then
        model='mbpt4sdq'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mbpt4sdq(t)')) then
        if (module.eq.'tce') then
        model='mbpt4sdq_t'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mp2')) then
        if (module.eq.'tce') then
        model='mbpt2'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mp3')) then
        if (module.eq.'tce') then
        model='mbpt3'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mp4sdq')) then
        if (module.eq.'tce') then
        model='mbpt4sdq'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1       INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mp4sdq(t)')) then
        if (module.eq.'tce') then
        model='mbpt4sdq_t'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1       INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mp4')) then
        if (module.eq.'tce') then
        model='mbpt4'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
c
c     THRESH (convergence threshold for Davidson iteration)
c
      else if (inp_compare(.false.,test,'thresh')) then
        if (.not.inp_f(thresh)) then
          write(LuOut,*) 'tce_input: thresh value not found; ',
     1      'default value of 1e-6 will be used'
          thresh=1.0d-6
        endif
        if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     LEVEL SHIFT
c
      else if (inp_compare(.false.,test,'lshift')) then
        if (.not.inp_f(zlshift)) then
          write(LuOut,*) 'tce_input: lshift value not found; ',
     1      'default value of 0.0d0 will be used'
          zlshift=0.0d0
        endif
        if (.not.rtdb_put(rtdb,'tce:zlshift',mt_dbl,1,zlshift))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'lshiftl')) then
        if (.not.inp_f(zlshiftl)) then
          write(LuOut,*) 'tce_input: lshiftl value not found; ',
     1      'default value of 0.0d0 will be used'
          zlshiftl=0.0d0
        endif
        if (.not.rtdb_put(rtdb,'tce:zlshiftl',mt_dbl,1,zlshiftl))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'lshift2')) then
        if (.not.inp_f(zlshift2(1))) then
          write(LuOut,*) 'tce_input: lshift2(1) value not found; ',
     1      'default value of 0.0d0 will be used'
          zlshift2(1)=0.0d0
        endif
        if (.not.inp_f(zlshift2(2))) then
          write(LuOut,*) 'tce_input: lshift2(2) value not found; ',
     1      'default value of 0.0d0 will be used'
          zlshift2(2)=0.0d0
        endif
        if (.not.rtdb_put(rtdb,'tce:zlshift2',mt_dbl,2,zlshift2))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'lshift3')) then
        if (.not.inp_f(zlshift3(1))) then
          write(LuOut,*) 'tce_input: lshift3(1) value not found; ',
     1      'default value of 0.0d0 will be used'
          zlshift3(1)=0.0d0
        endif
        if (.not.inp_f(zlshift3(2))) then
          write(LuOut,*) 'tce_input: lshift3(2) value not found; ',
     1      'default value of 0.0d0 will be used'
          zlshift3(2)=0.0d0
        endif
        if (.not.rtdb_put(rtdb,'tce:zlshift3',mt_dbl,2,zlshift3))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      else if (inp_compare(.false.,test,'maxiter')) then
        if (.not.inp_i(maxiter)) then
          write(LuOut,*) 'tce_input: maxiter value not found; ',
     1      'default value of 100 will be used'
          maxiter=100
        endif
        if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     IOALGORITHM (I/O method)
c
      else if (inp_compare(.false.,test,'io')) then
        if (.not.inp_a(ioalgchar)) then
          write(LuOut,*) 'tce_input: ioalgorithm value not found; ',
     1      'default GA fully incore algorithm will be used'
          ioalg=2
        else
          if (ioalgchar.eq.'fortran') then
            ioalg=0
          else if (ioalgchar.eq.'eaf') then
            ioalg=1
          else if (ioalgchar.eq.'ga') then
            ioalg=2
          else if (ioalgchar.eq.'sf') then
            ioalg=3
          else if (ioalgchar.eq.'replicated') then
            ioalg=4
          else if (ioalgchar.eq.'dra') then
            ioalg=5
          else if (ioalgchar.eq.'ga_eaf') then
            ioalg=6
          endif
        endif
        if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c EOMCC SOLVER
c
      else if (inp_compare(.false.,test,'eomsol')) then
        if (.not.inp_i(eomsol)) then
          write(LuOut,*) 'tce_input: eomsol value not found; ',
     1      'default value of 1 will be used'
          eomsol=1
        endif
        if (.not.rtdb_put(rtdb,'tce:eoms',mt_int,1,eomsol))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     DIIS (the vector space size in DIIS)
c
      else if (inp_compare(.false.,test,'diis')) then
        if (.not.inp_i(diis)) then
          write(LuOut,*) 'tce_input: diis value not found; ',
     1      'default value of 5 will be used'
          diis=5
        endif
        if (.not.rtdb_put(rtdb,'tce:diis',mt_int,1,diis))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'diis2')) then
        if (.not.inp_i(diis2)) then
          write(LuOut,*) 'tce_input: diis2 value not found; ',
     1      'default value of 5 will be used'
          diis2=5
        endif
        if (.not.rtdb_put(rtdb,'tce:diis2',mt_int,1,diis2))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'diis3')) then
        if (.not.inp_i(diis3)) then
          write(LuOut,*) 'tce_input: diis3 value not found; ',
     1      'default value of 5 will be used'
          diis3=5
        endif
        if (.not.rtdb_put(rtdb,'tce:diis3',mt_int,1,diis3))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c DIMENSION OF EOMCC ITERATIVE SPACE
c
      else if (inp_compare(.false.,test,'hbard')) then
        if (.not.inp_i(hbard)) then
          write(LuOut,*) 'tce_input: hbard value not found; ',
     1      'default value of 500 will be used'
          hbard=500
        endif
        if (.not.rtdb_put(rtdb,'tce:hbard',mt_int,1,hbard))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     NROOTS (the number of excited state root)
c
      else if (inp_compare(.false.,test,'nroots')) then
        if (.not.inp_i(nroots)) then
          write(LuOut,*) 'tce_input: nroots value not found; ',
     1      'default value of 0 will be used'
          nroots=0
        endif
        if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     MAXDIFF (for EOM codes)
c
      else if (inp_compare(.false.,test,'maxdiff')) then
        if (.not.inp_f(maxdiff)) then
c          write(LuOut,*) 'tce_input: maxdiff value not found; ',
c     1      'default value of 1e-6 will be used'
          maxdiff=0.5d0
        endif
        if (.not.rtdb_put(rtdb,'tce:maxdiff',mt_dbl,1,maxdiff))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     2e STORAGE
c
      else if (inp_compare(.false.,test,'attilesize')) then
        if (.not.inp_i(maxs)) then
          write(LuOut,*) 'tce_input: attilesize value not found; ',
     1      'default value of 30 will be used'
          maxs=30
        endif
        if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,maxs))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
      else if (inp_compare(.false.,test,'split')) then
        if (.not.inp_i(ichopx)) then
          write(LuOut,*) 'tce_input: split value not found; ',
     1      'default value of 1 will be used'
          ichopx=1
        endif
        if (.not.rtdb_put(rtdb,'tce:ichopx',mt_int,1,ichopx))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
      else if (inp_compare(.false.,test,'2emet')) then
        if (.not.inp_i(i4im)) then
          write(LuOut,*) 'tce_input: 2emet value not found; ',
     1      'default value of 1 will be used'
          i4im=1
        endif
        if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,i4im))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
      else if (inp_compare(.false.,test,'idiskx')) then
        if (.not.inp_i(idiskx)) then
          write(LuOut,*) 'tce_input: idiskx value not found; ',
     1      'default value of 0 will be used'
          idiskx=0
        endif
        if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,idiskx))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     TARGET (the target excited state for, e.g., geometry optimization)
c
      else if (inp_compare(.false.,test,'target')) then
        if (.not.inp_i(target)) then
          write(LuOut,*) 'tce_input: target value not found; ',
     1      'default value of 1 will be used'
          target=1
        endif
        if (target.gt.nroots) call errquit
     1    ('tce_input: an illegal value for target',target,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:target',mt_int,1,target))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     TARGETSYM (the symmetry of the target excited state)
c
      else if (inp_compare(.false.,test,'targetsym')) then
        if (.not.inp_a(targetsym)) then
          write(LuOut,*) 'tce_input: targetsym value not found; ',
     1      'no symmetry information will be used in specifying target'
          targetsym='none'
        endif
        if (.not.rtdb_cput(rtdb,'tce:targetsym',1,targetsym))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
c
      else if (inp_compare(.false.,test,'symmetry')) then
        symmetry=.true.
        if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     IDENS (one particle reduced density matrix)
c
      else if (inp_compare(.false.,test,'densmat')) then
         idens=.true.
         left =.true.
         if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1      call errquit('tce_input: failed writing to rtdb',0,
     2      RTDB_ERR)
         if (.not.rtdb_put(rtdb,'tce:densmat',mt_log,1,idens))
     1      call errquit('tce_input: failed writing to rtdb',0,
     2      RTDB_ERR)
         if (.not.inp_a(file_densmat)) then
            call util_file_name('densmat', .false.,.false.,file_densmat)
         endif
         if (.not.rtdb_cput(rtdb,'tce:file_densmat',1,file_densmat))
     1      call errquit('tce_input: rtdb_cput failed - file_densmat',0,
     1      RTDB_ERR)
c
c     MULTIPOLE (multipole moments)
c
      else if (inp_compare(.false.,test,'multipole')) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.inp_i(multipole)) then
          write(LuOut,*) 'tce_input: multipole value not found; ',
     1      'all available multipoles (L=1,2,3) will be calculated'
          multipole=3
        endif
        if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     DIPOLE (dipole moments & dipole transition moments)
c     QUADRUPOLE (quadrupole moments & quadrupole transition moments)
c     OCTUPOLE (octupole moments & octupole transition moments)
c
      else if (inp_compare(.false.,test,'dipole')) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        multipole=max(multipole,1)
        if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'quadrupole')) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        multipole=max(multipole,2)
        if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'octupole')) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        multipole=max(multipole,3)
        if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     PROPERTY INPUT SUB-BLOCK
c
      else if (inp_compare(.false.,test,'tceprop')) then
        call tce_prop_input(rtdb)
c
c     TILESIZE (the maximum tile size)
c
      else if (inp_compare(.false.,test,'tilesize')) then
        if (.not.inp_i(tilesize))
     1    call errquit('tce_input: no tilesize given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:tilesize',mt_int,1,tilesize))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     FRAGMENT (if excited state calc, give an atom in an excited fragment)
c
      else if (inp_compare(.false.,test,'fragment')) then
        if (.not.inp_i(fragment)) then
          write(LuOut,*) 'tce_input: fragment value not found; ',
     1      'default value of 0 will be used'
          fragment=0
        endif
        if (.not.rtdb_put(rtdb,'tce:fragment',mt_int,1,fragment))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     (NO)FOCK (recompute fock for, e.g., DFT, ROHF refs)
c
      else if (inp_compare(.false.,test,'fock')) then
        recompf=.true.
        if (.not.rtdb_put(rtdb,'tce:recompf',mt_log,1,recompf))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'nofock')) then
        recompf=.false.
        if (.not.rtdb_put(rtdb,'tce:recompf',mt_log,1,recompf))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     ACTIVE_OA (Number of active occupied orbitals)
c
      else if (inp_compare(.false.,test,'active_oa')) then
        if (.not.inp_i(oactive(1)))
     1    call errquit('tce_input: no active_oa given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_oa',mt_int,1,oactive(1)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     ACTIVE_OB (Number of active occupied orbitals)
c
      else if (inp_compare(.false.,test,'active_ob')) then
        if (.not.inp_i(oactive(2)))
     1    call errquit('tce_input: no active_ob given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_ob',mt_int,1,oactive(2)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     ACTIVE_VA (Number of active virtual orbitals)
c
      else if (inp_compare(.false.,test,'active_va')) then
        if (.not.inp_i(vactive(1)))
     1    call errquit('tce_input: no active_va given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_va',mt_int,1,vactive(1)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     ACTIVE_VB (Number of active virtual orbitals)
c
      else if (inp_compare(.false.,test,'active_vb')) then
        if (.not.inp_i(vactive(2)))
     1    call errquit('tce_input: no active_vb given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_vb',mt_int,1,vactive(2)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     ACTIVE_EXCIT_LVL (T3 active excitation level)
c
      else if (inp_compare(.false.,test,'t3a_lvl')) then
        if (.not.inp_i(numact))
     1    call errquit('tce_input: no t3a_lvl given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:act_excit_lvl',mt_int,1,numact))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     PRINT
c
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'tce')
c
c     END
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit('tce_input: unknown directive',0,INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
c
c     This is the TCE property input block ("tceprop")
c
      subroutine tce_prop_input(rtdb)
c
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer rtdb
      integer n_a,k_a,l_a
      integer n_b,k_b,l_b
      integer n_c,k_c,l_c
      integer n_i,k_i,l_i
      integer i,icount
      character*10 module
      character*20 test
      character*20 beta_opt
      character*20 gamm_opt
      character*20 disp_opt
      character*20 beta_type
      character*20 gamm_type
      character*20 disp_type
      logical lineresp ! T(1) response equations - real frequency
      logical leftresp ! L(1) response equations 
      logical quadresp ! T(2) response equations
      logical imagresp ! T(1) response equatiosn - imaginary frequency
      logical status
c
      lineresp = .false.
      leftresp = .false.
      quadresp = .false.
      imagresp = .false.
c
c -------------------------------------
c What input block are we dealing with?
c -------------------------------------
c
      if (.not.rtdb_cget(rtdb,'tce:module',1,module)) then
        call errquit('tce_prop_input: line ',__LINE__,RTDB_ERR)
      endif
c
c ----------
c Read input
c ----------
c
 100  if (.not. inp_read()) then
        call errquit('tce_prop_input: line ',__LINE__,RTDB_ERR)
      endif
      if (.not. inp_a(test)) then
        call errquit('tce_prop_input: line ',__LINE__,RTDB_ERR)
      endif
c
c     POLARIZABILITY
c
      if (inp_compare(.false.,test,'polarizability').or.
     1    inp_compare(.false.,test,'polar').or.
     2    inp_compare(.false.,test,'alpha')) then
        lineresp = .true.
c
c     HYPERPOLARIZABILITY
c
      elseif (inp_compare(.false.,test,'hyperpolarizability').or.
     1    inp_compare(.false.,test,'hyperpolar').or.
     2    inp_compare(.false.,test,'beta')) then
        lineresp = .true.
        leftresp = .true.
        if (.not.inp_a(beta_opt)) then
          beta_type = 'static'
        else
          if (inp_compare(.false.,beta_opt,'shg')) then
            beta_type = 'SHG'
          elseif (inp_compare(.false.,beta_opt,'or')) then
            beta_type = 'OR'
          elseif (inp_compare(.false.,beta_opt,'eope')) then
            beta_type = 'EOPE'
          else
            call errquit('tce_prop_input: invalid option for beta',
     1                   __LINE__,RTDB_ERR)
          endif
        endif
c
c     SECOND HYPERPOLARIZABILITY
c
      elseif (inp_compare(.false.,test,'cubicpolarizability').or.
     1        inp_compare(.false.,test,'cubicpolar').or.
     2        inp_compare(.false.,test,'gamma')) then
        lineresp = .true.
        leftresp = .true.
        quadresp = .true.
        if (.not.inp_a(gamm_opt)) then
          gamm_type = 'static'
        else
          if (inp_compare(.false.,gamm_opt,'thg')) then
            gamm_type = 'THG'
          elseif (inp_compare(.false.,gamm_opt,'efish')) then
            gamm_type = 'EFISH'
          elseif (inp_compare(.false.,gamm_opt,'dfwm')) then
            gamm_type = 'DFWM'
          elseif (inp_compare(.false.,gamm_opt,'oke')) then
            gamm_type = 'OKE'
          elseif (inp_compare(.false.,gamm_opt,'cars')) then
            gamm_type = 'CARS'
          else
            call errquit('tce_prop_input: invalid option for gamma',
     1                   __LINE__,RTDB_ERR)
          endif
        endif
c
c     DISPERSION
c
      elseif (inp_compare(.false.,test,'dispersion').or.
     1        inp_compare(.false.,test,'disp').or.
     2        inp_compare(.false.,test,'c6')) then
        imagresp = .true.
        if (.not.inp_a(disp_opt)) then
          disp_type = 'derevianko'
        else
          if (inp_compare(.false.,disp_opt,'custom')) then
            disp_type = 'custom'
          elseif (inp_compare(.false.,disp_opt,'gauss-legendre')) then
            disp_type = 'gauss-leg'
          elseif (inp_compare(.false.,disp_opt,'gauss-chebyshev')) then
            disp_type = 'gauss-cheb'
          else
            call errquit('tce_prop_input: invalid option for disp',
     1                   __LINE__,RTDB_ERR)
          endif
        endif
c
c     AFREQ (omega for polarizability)
c
      elseif (inp_compare(.false.,test,'afreq')) then
        lineresp = .true.
        if (inp_i(n_a)) then
          if (.not.ma_push_get(mt_dbl,n_a,'afreq',l_a,
     1                         k_a)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          icount = 0
          do i = 0, n_a-1
            status = inp_f(dbl_mb(k_a+i))
            if (status) then
              icount = icount+1
            else
              write(6,'(a,a,i4,a)')
     1           'Response property input found ',
     2           'fewer frequencies than expected, only ',icount,
     3           'will be used'
              n_a = icount
              if (icount.eq.0) then
                if (.not.ma_pop_stack(l_a)) then
                  call errquit('tce_prop_input: ma_pop_stack',__LINE__,
     1                         MA_ERR)
                endif
              endif
              goto 300
            endif
          enddo
        else
          n_a = 1
          if (.not.ma_push_get(mt_dbl,n_a,'afreq',l_a,k_a)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          dbl_mb(k_a) = 0.0d0
        endif
 300    continue
c
c     BFREQ (omega for first hyperpolarizability)
c
      elseif (inp_compare(.false.,test,'bfreq')) then
        lineresp = .true.
        leftresp = .true.
        if (inp_i(n_b)) then
          if (.not.ma_push_get(mt_dbl,n_b,'bfreq',l_b,
     1                         k_b)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          icount = 0
          do i = 0, n_b-1
            status = inp_f(dbl_mb(k_b+i))
            if (status) then
              icount = icount+1
            else
              write(6,'(a,a,i4,a)')
     1           'Response property input found ',
     2           'fewer frequencies than expected, only ',icount,
     3           'will be used'
              n_b = icount
              if (icount.eq.0) then
                if (.not.ma_pop_stack(l_b)) then
                  call errquit('tce_prop_input: ma_pop_stack',__LINE__,
     1                         MA_ERR)
                endif
              endif
              goto 400
            endif
          enddo
        else
          n_b = 1
          if (.not.ma_push_get(mt_dbl,n_b,'bfreq',l_b,k_b)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          dbl_mb(k_b) = 0.0d0
        endif
 400    continue
c
c     CFREQ (omega for second hyperpolarizability)
c
      elseif (inp_compare(.false.,test,'cfreq')) then
        lineresp = .true.
        leftresp = .true.
        quadresp = .true.
        if (inp_i(n_c)) then
          if (.not.ma_push_get(mt_dbl,n_c,'cfreq',l_c,
     1                         k_c)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          icount = 0
          do i = 0, n_b-1
            status = inp_f(dbl_mb(k_c+i))
            if (status) then
              icount = icount+1
            else
              write(6,'(a,a,i4,a)')
     1           'Response property input found ',
     2           'fewer frequencies than expected, only ',icount,
     3           'will be used'
              n_c = icount
              if (icount.eq.0) then
                if (.not.ma_pop_stack(l_c)) then
                  call errquit('tce_prop_input: ma_pop_stack',__LINE__,
     1                         MA_ERR)
                endif
              endif
              goto 500
            endif
          enddo
        else
          n_c = 1
          if (.not.ma_push_get(mt_dbl,n_c,'cfreq',l_c,k_c)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          dbl_mb(k_c) = 0.0d0
        endif
 500    continue
c
c     IFREQ (omega for imaginary polarizability)
c
      elseif (inp_compare(.false.,test,'ifreq')) then
        imagresp = .true.
        if (inp_i(n_i)) then
          if (.not.ma_push_get(mt_dbl,n_i,'ifreq',l_i,
     1                         k_i)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          icount = 0
          do i = 0, n_i-1
            status = inp_f(dbl_mb(k_i+i))
            if (status) then
              icount = icount+1
            else
              write(6,'(a,a,i4,a)')
     1           'Response property input found ',
     2           'fewer frequencies than expected, only ',icount,
     3           'will be used'
              n_i = icount
              if (icount.eq.0) then
                if (.not.ma_pop_stack(l_c)) then
                  call errquit('tce_prop_input: ma_pop_stack',__LINE__,
     1                         MA_ERR)
                endif
              endif
              goto 600
            endif
          enddo
        else
          n_i = 1
          if (.not.ma_push_get(mt_dbl,n_i,'ifreq',l_i,k_i)) then
             call errquit('tce_prop_input: ma_push_get',__LINE__,MA_ERR)
          endif
          dbl_mb(k_i) = 0.0d0
        endif
 600    continue
c
c     END
c
      else if (inp_compare(.false.,test,'end')) then
        goto 200
      else
        call errquit('tce_prop_input: unknown directive',0,INPUT_ERR)
      endif
      goto 100
 200  return
c
c -------------------
c Push values to RTDB
c -------------------
c
      if (.not.rtdb_put(rtdb,'tce:lineresp',mt_log,1,lineresp)) then
        call errquit('tce_prop_input: rtdb_put',__LINE__,RTDB_ERR)
      endif
      if (.not.rtdb_put(rtdb,'tce:leftresp',mt_log,1,leftresp)) then
        call errquit('tce_prop_input: rtdb_put',__LINE__,RTDB_ERR)
      endif
c
c ------
c Return
c ------
c
      end

