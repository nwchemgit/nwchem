#ifndef __LINE__
#define __LINE__ 0
#endif
      subroutine tce_input(rtdb)
!
! $Id$
!
! Input parser for TCE module for various many-electron theories.
! Also sets default values for input parameters.
! Modified from tddft/tddft_input.F by So Hirata Oct, 2002.
!
!     TCE
!        [(DFT||HF||SCF) default HF]
!        [FREEZE [[core] (atomic || <integer nfzc default 0>)]
!                 [virtual <integer nfzv default 0>]]
!        [(LCCD||CCD||CCSD||LCCSD||CCSDT||CCSDTQ||
!          CCSD(T)||CCSD[T]||QCISD||CISD||CISDT||CISDTQ||
!          MBPT2||MBPT3||MBPT4||MP2||MP3||MP4|| \
!          CR-CCSD(T)||CR-CCSD[T]||LR-CCSD(T)||LR-CCSD(TQ)||CCSD(2)_T||CCSD(2)||
!          CCSDT(2)_Q) default CCSD]
!        [THRESH <double thresh default 1e-6>]
!        [MAXITER <integer maxiter default 100>]
!        [PRINT (none||low||medium||high||debug)]
!        [IO (fortran||eaf||ga||sf||replicated||dra||ga_eaf) default ga]
!        [DIIS <integer diis default 5>]
!        [EOMSOL <integer default 1 >]
!        [DIIS2 <integer diis default 5>]
!        [DIIS3 <integer diis default 5>]
!        [NROOTS <integer nroots default 0>]
!        [TARGET <integer target default 1>]
!        [TARGETSYM <character targetsym default 'none'>]
!        [SYMMETRY]
!        [DIPOLE]
!        [TILESIZE <no default (automatically adjusted)>]
!        [FRAGMENT <default -1 (off)>]
!        [(NO)FOCK <logical recompf default .true.>]
!        [ACTIVE_OA <default 0>]
!        [ACTIVE_OB <default 0>]
!        [ACTIVE_VA <default 0>]
!        [ACTIVE_VB <default 0>]
!        [T3A_LVL   <default 0>]
!     END
!
!     TASK TCE ENERGY
!
!     ... or ...
!
!     UCCSDT or UCC or UCCSD(T) etc.
!        [(DFT||HF||SCF) default HF]
!        [FREEZE [[core] (atomic || <integer nfzc default 0>)]
!                 [virtual <integer nfzv default 0>]]
!        [THRESH <double thresh default 1e-6>]
!        [MAXITER <integer maxiter default 100>]
!        [PRINT (none||low||medium||high||debug)]
!        [IO (fortran||c||ga||sf||replicated) default ga]
!        [DIIS <integer diis default 5>]
!        [NROOTS <integer nroots default 0>]
!        [TARGET <integer target default 1>]
!        [TARGETSYM <character targetsym default 'none'>]
!        [SYMMETRY]
!        [DIPOLE]
!        [TILESIZE <no default (automatically adjusted)>]
!        [FRAGMENT <default -1 (off)>]
!        [(NO)FOCK <logical recompf default .true.>]
!        [ACTIVE_OA <default 0>]
!        [ACTIVE_OB <default 0>]
!        [ACTIVE_VA <default 0>]
!        [ACTIVE_VB <default 0>]
!     END
!
!     TASK UCCSDT ENERGY
!
!     ... etc.
!     
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
! --- level shift --
      double precision zlshift,zlshiftl,zlshift2(2),zlshift3(2)
! ------------------
      integer nroots
      integer target
      integer tilesize
      integer fragment
      character*4 targetsym
      logical symmetry
      logical left
! --- density matrix    
      logical idens
      character*256 file_densmat
!<-d3p975
      integer multipole
      logical recompf
      character*10 perturbative
      character*10 ccsd_var 
      integer oactive(2)
      integer vactive(2)
      integer numact
! --- ccsd_act/eomccsd_act ---
      integer uact,oact
      double precision emin_act,emax_act
! --- 4 index transform. ---
      integer maxs,ichopx,i4im,idiskx
! --- EOM solver
      integer hbard

! --- TCE_CUDA
      integer icuda
!kbn --- EA/IPCCSD
#ifdef EACCSD
      logical eaccsd
#endif
#ifdef IPCCSD
      logical ipccsd
#endif

!
! -------------------------------------
! What input block are we dealing with?
! -------------------------------------
!
      if (.not.rtdb_cget(rtdb,'tce:module',1,module))
     1  call errquit('tce_input: failed reading from rtdb',0,
     2  RTDB_ERR)
!
! ------------------
! Set default values
! ------------------
!
!     DFT, HF, or SCF (reference wavefunction)
!
      reference=1
      if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     FREEZE (frozen cores/virtuals)
!
!     no action is taken
!
!     MODEL (the name of CC model requested)
!
!     no action is taken
!
!     THRESH (convergence threshold for Davidson iteration)
!
      thresh=1.0d-7
      if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     LEVEL SHIFT (for singles and doubles)
!
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
!
!     MAXITER (the maximum number of Davidson iterations)
!
      maxiter=100
      if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     IO (I/O method, 0 = Fortran Direct Access, 
!                     1 = C Low-Level I/O,
!                     2 = GA Library,
!                     3 = SF library,
!                     4 = Replicated C Low-Level I/O)
!
      ioalg=2
      if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     DIIS (the vector space size in DIIS)
!
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
!
!     EOMCC SOLVER
!
      eomsol=1
      if (.not.rtdb_put(rtdb,'tce:eoms',mt_int,1,eomsol))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     DIMENSION OF THE EOM ITERATIVE SPACE
!
      hbard=500
      if (.not.rtdb_put(rtdb,'tce:hbard',mt_int,1,hbard))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     NROOTS (the number of excited state roots)
!
      nroots=0
      if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     2e STORAGE
!
      model2e='default'
      if (.not.rtdb_cput(rtdb,'tce:model2e',1,model2e))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
! 4ind. transfromation
!
       maxs=30
      if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,maxs))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR) 
!
       ichopx=1
      if (.not.rtdb_put(rtdb,'tce:ichopx',mt_int,1,ichopx))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
       i4im=1
      if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,i4im))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
      idiskx=0
      if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,idiskx))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     TARGET (the target excited state for, e.g., geometry optimization)
!
      target=1
      if (.not.rtdb_put(rtdb,'tce:target',mt_int,1,target))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     TARGETSYM (the irrep of the target excited state)
!
      targetsym='none'
      if (.not.rtdb_cput(rtdb,'tce:targetsym',1,targetsym))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
!
      symmetry=.false.
      if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     IDENS (one particle reduced density matrix)
!
      idens=.false.
      if (.not.rtdb_put(rtdb,'tce:densmat',mt_log,1,idens))
     1   call errquit('tce_input: failed writing to rtdb',0,
     1   rtdb_err)
!
!     DIPOLE (dipole moments & dipole transition moments)
!
      left=.false.
      if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     MULTIPOLE LMAX (multipole moments highest angular momentum)
!
      multipole=0
      if (.not.rtdb_put(rtdb,'tce:multipole',mt_int,1,multipole))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     FRAGMENT (fragment MO calculations)
!
      fragment=-1
      if (.not.rtdb_put(rtdb,'tce:fragment',mt_int,1,fragment))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     (NO)FOCK (recompute fock for, e.g., DFT, ROHF refs)
!
      recompf=.true.
      if (.not.rtdb_put(rtdb,'tce:recompf',mt_log,1,recompf))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     ACTIVE_OA,OB (Number of active occupied orbitals)
!
      oactive(1)=0
      if (.not.rtdb_put(rtdb,'tce:active_oa',mt_int,1,oactive(1)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      oactive(2)=0
      if (.not.rtdb_put(rtdb,'tce:active_ob',mt_int,1,oactive(2)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!     ccsd_act/eomccsd_act
!
      oact=0
      if (.not.rtdb_put(rtdb,'tce:oact',mt_int,1,oact))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      uact=0
      if (.not.rtdb_put(rtdb,'tce:uact',mt_int,1,uact))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
      emin_act=0.0d0
      if (.not.rtdb_put(rtdb,'tce:eactmin',mt_dbl,1,emin_act))
     1  call errquit('tce_input: rtdb eactmin problem',0,
     2  RTDB_ERR)
      emax_act=0.0d0
      if (.not.rtdb_put(rtdb,'tce:eactmax',mt_dbl,1,emax_act))
     1  call errquit('tce_input: rtdb eactmax problem',0,
     2  RTDB_ERR)
!
!     ACTIVE_VA,VB (Number of active virtual orbitals)
!
      vactive(1)=0
      if (.not.rtdb_put(rtdb,'tce:active_va',mt_int,1,vactive(1)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
      vactive(2)=0
      if (.not.rtdb_put(rtdb,'tce:active_vb',mt_int,1,vactive(2)))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
!    ACTIVE EXCITATION LEVEL (number of active orbitals in T3)
!
      numact=0
      if (.not.rtdb_put(rtdb,'tce:act_excit_lvl',mt_int,1,numact))
     1  call errquit('tce_input: failed writing to rtdb',0,
     2  RTDB_ERR)
!
! ----------
! Read input
! ----------
!
   10 if (.not. inp_read()) 
     1  call errquit('tce_input: failed reading input',0,
     2  RTDB_ERR)
      if (.not. inp_a(test)) 
     1  call errquit('tce_input: failed reading keyword',0,
     2  RTDB_ERR)
!
!     DFT, HF, or SCF (reference wavefunction)
!
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
!
!     FREEZE (frozen cores/virtuals)
!
      else if (inp_compare(.false.,test,'freeze')) then
        call freeze_input(rtdb,'tce')
!
!     STORAGE OF 2-e INTEGRALS
!
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
      else if (inp_compare(.false.,test,'2espin')) then
        if (module.eq.'tce') then
        model2e='2espin'
        if (.not.rtdb_cput(rtdb,'tce:model2e',1,model2e))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
!
!     MODEL (the name of theory requested)
!
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
! ccsd_act/eomccsd_act
      else if (inp_compare(.false.,test,'ccsd_act')) then
        if (module.eq.'tce') then
        model='ccsd_act'
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
! ccsd_act/eomccsd_act
      else if (inp_compare(.false.,test,'crsd(t)ac')) then
        if (module.eq.'tce') then
        model='ccsd_act'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='cr_(t)a'
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
        if (.not.rtdb_put(rtdb,'tce:nts',mt_log,1,.false.))
     1    call errquit('tce_input: failed writing nts to rtdb',0,
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
!
!
!
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
!      else if (inp_compare(.false.,test,'lambda-ccsd[t]')) then
!        if (module.eq.'tce') then
!        model='ccsd'
!        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
!     1    call errquit('tce_input: failed writing to rtdb',0,
!     2    RTDB_ERR)
!        perturbative='lambda[t]'
!        if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbative))
!     1    call errquit('tce_input: failed writing to rtdb',0,
!     2    RTDB_ERR)
!        left=.true.
!        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
!     1    call errquit('tce_input: failed writing to rtdb',0,
!     2    RTDB_ERR)
!        else
!        call errquit('tce_input: multiple theory inputs',0,
!     1    INPUT_ERR)
!        endif
!
!
!
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
! ccsd_act/eomccsd-act
      else if (inp_compare(.false.,test,'creom(t)ac')) then
        if (module.eq.'tce') then
        model='ccsd_act'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        perturbative='creom(t)a'
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

#ifdef MRCC_METHODS
!kbn mrcc-r-1 -3
!     BWCCSD
      else if (inp_compare(.false.,test,'bwccsd')) then
        if (module.eq.'tce') then
        model='bwccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:mrcc',mt_int,1,1))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
      else if (inp_compare(.false.,test,'mkccsd')) then
!     MkCCSD
        if (module.eq.'tce') then
        model='mkccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:mrcc',mt_int,1,1))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        else
        call errquit('tce_input: multiple theory inputs',0,
     1    INPUT_ERR)
        endif
#endif

!
!     THRESH (convergence threshold for Davidson iteration)
!
      else if (inp_compare(.false.,test,'thresh')) then
        if (.not.inp_f(thresh)) then
          write(LuOut,*) 'tce_input: thresh value not found; ',
     1      'default value of 1e-6 will be used'
          thresh=1.0d-6
        endif
        if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     LEVEL SHIFT
!
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

!
!      TCE_CUDA Number of CUDA devices per node
!
      else if (inp_compare(.false.,test,'cuda')) then
#if defined(TCE_CUDA) || defined(TCE_HIP)
        if (.not.inp_i(icuda))
     1    call errquit('tce_input: no icuda',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:cuda',mt_int,1,icuda))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)

#else
        call errquit('cuda option needs TCE_CUDA compiled code',
     A       0, INPUT_ERR)
#endif
!
!     MAXITER (the maximum number of Davidson iterations)
!
      else if (inp_compare(.false.,test,'maxiter')) then
        if (.not.inp_i(maxiter)) then
          write(LuOut,*) 'tce_input: maxiter value not found; ',
     1      'default value of 100 will be used'
          maxiter=100
        endif
        if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     IOALGORITHM (I/O method)
!
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
!
! EOMCC SOLVER
!
      else if (inp_compare(.false.,test,'eomsol')) then
        if (.not.inp_i(eomsol)) then
          write(LuOut,*) 'tce_input: eomsol value not found; ',
     1      'default value of 1 will be used'
          eomsol=1
        endif
        if (.not.rtdb_put(rtdb,'tce:eoms',mt_int,1,eomsol))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     DIIS (the vector space size in DIIS)
!
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
!
! DIMENSION OF EOMCC ITERATIVE SPACE
!
      else if (inp_compare(.false.,test,'hbard')) then
        if (.not.inp_i(hbard)) then
          write(LuOut,*) 'tce_input: hbard value not found; ',
     1      'default value of 500 will be used'
          hbard=500
        endif
        if (.not.rtdb_put(rtdb,'tce:hbard',mt_int,1,hbard))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     NROOTS (the number of excited state root)
!
      else if (inp_compare(.false.,test,'nroots')) then
        if (.not.inp_i(nroots)) then
          write(LuOut,*) 'tce_input: nroots value not found; ',
     1      'default value of 0 will be used'
          nroots=0
        endif
        if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!kbn EACCSD
#ifdef EACCSD
      else if (inp_compare(.false.,test,'eaccsd')) then
        eaccsd=.true.
        if (.not.rtdb_put(rtdb,'tce:eaccsd',mt_log,1,eaccsd))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
#endif
!
!kbn IPCCSD
#ifdef IPCCSD
      else if (inp_compare(.false.,test,'ipccsd')) then
        ipccsd=.true.
        if (.not.rtdb_put(rtdb,'tce:ipccsd',mt_log,1,ipccsd))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
#endif
!
!     MAXDIFF (for EOM codes)
!
      else if (inp_compare(.false.,test,'maxdiff')) then
        if (.not.inp_f(maxdiff)) then
!          write(LuOut,*) 'tce_input: maxdiff value not found; ',
!     1      'default value of 1e-6 will be used'
          maxdiff=0.5d0
        endif
        if (.not.rtdb_put(rtdb,'tce:maxdiff',mt_dbl,1,maxdiff))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     2e STORAGE
!
      else if (inp_compare(.false.,test,'attilesize')) then
        if (.not.inp_i(maxs)) then
          write(LuOut,*) 'tce_input: attilesize value not found; ',
     1      'default value of 30 will be used'
          maxs=30
        endif
        if (.not.rtdb_put(rtdb,'tce:maxs',mt_int,1,maxs))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
      else if (inp_compare(.false.,test,'split')) then
        if (.not.inp_i(ichopx)) then
          write(LuOut,*) 'tce_input: split value not found; ',
     1      'default value of 1 will be used'
          ichopx=1
        endif
        if (.not.rtdb_put(rtdb,'tce:ichopx',mt_int,1,ichopx))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
      else if (inp_compare(.false.,test,'2emet')) then
        if (.not.inp_i(i4im)) then
          write(LuOut,*) 'tce_input: 2emet value not found; ',
     1      'default value of 1 will be used'
          i4im=1
        endif
        if (.not.rtdb_put(rtdb,'tce:i4im',mt_int,1,i4im))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
      else if (inp_compare(.false.,test,'idiskx')) then
        if (.not.inp_i(idiskx)) then
          write(LuOut,*) 'tce_input: idiskx value not found; ',
     1      'default value of 0 will be used'
          idiskx=0
        endif
        if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,idiskx))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     TARGET (the target excited state for, e.g., geometry optimization)
!
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
!
!     TARGETSYM (the symmetry of the target excited state)
!
      else if (inp_compare(.false.,test,'targetsym')) then
        if (.not.inp_a(targetsym)) then
          write(LuOut,*) 'tce_input: targetsym value not found; ',
     1      'no symmetry information will be used in specifying target'
          targetsym='none'
        endif
        if (.not.rtdb_cput(rtdb,'tce:targetsym',1,targetsym))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
!
      else if (inp_compare(.false.,test,'symmetry')) then
        symmetry=.true.
        if (.not.rtdb_put(rtdb,'tce:symmetry',mt_log,1,symmetry))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     IDENS (one particle reduced density matrix)
!
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
!
!     MULTIPOLE (multipole moments)
!
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
!
!     DIPOLE (dipole moments & dipole transition moments)
!     QUADRUPOLE (quadrupole moments & quadrupole transition moments)
!     OCTUPOLE (octupole moments & octupole transition moments)
!
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
!
!     PROPERTY INPUT SUB-BLOCK
!
      else if (inp_compare(.false.,test,'tceprop')) then
        call tce_prop_input(rtdb)
!
!     TILESIZE (the maximum tile size)
!
      else if (inp_compare(.false.,test,'tilesize')) then
        if (.not.inp_i(tilesize))
     1    call errquit('tce_input: no tilesize given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:tilesize',mt_int,1,tilesize))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     FRAGMENT (if excited state calc, give an atom in an excited fragment)
!
      else if (inp_compare(.false.,test,'fragment')) then
        if (.not.inp_i(fragment)) then
          write(LuOut,*) 'tce_input: fragment value not found; ',
     1      'default value of 0 will be used'
          fragment=0
        endif
        if (.not.rtdb_put(rtdb,'tce:fragment',mt_int,1,fragment))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     (NO)FOCK (recompute fock for, e.g., DFT, ROHF refs)
!
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
!
!    ccsd_act/eomccsd_act
!
      else if (inp_compare(.false.,test,'oact')) then
        if (.not.inp_i(oact))
     1    call errquit('tce_input: no oact given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:oact',mt_int,1,oact))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'uact')) then
        if (.not.inp_i(uact))
     1    call errquit('tce_input: no uact given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:uact',mt_int,1,uact))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
      else if (inp_compare(.false.,test,'emin_act')) then
        if (.not.inp_f(emin_act))
     1    call errquit('tce_input: no emin_act given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:eactmin',mt_dbl,1,emin_act))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'emax_act')) then
        if (.not.inp_f(emax_act))
     1    call errquit('tce_input: no emax_act given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:eactmax',mt_dbl,1,emax_act))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     ACTIVE_OA (Number of active occupied orbitals)
!
      else if (inp_compare(.false.,test,'active_oa')) then
        if (.not.inp_i(oactive(1)))
     1    call errquit('tce_input: no active_oa given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_oa',mt_int,1,oactive(1)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     ACTIVE_OB (Number of active occupied orbitals)
!
      else if (inp_compare(.false.,test,'active_ob')) then
        if (.not.inp_i(oactive(2)))
     1    call errquit('tce_input: no active_ob given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_ob',mt_int,1,oactive(2)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     ACTIVE_VA (Number of active virtual orbitals)
!
      else if (inp_compare(.false.,test,'active_va')) then
        if (.not.inp_i(vactive(1)))
     1    call errquit('tce_input: no active_va given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_va',mt_int,1,vactive(1)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     ACTIVE_VB (Number of active virtual orbitals)
!
      else if (inp_compare(.false.,test,'active_vb')) then
        if (.not.inp_i(vactive(2)))
     1    call errquit('tce_input: no active_vb given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:active_vb',mt_int,1,vactive(2)))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
!
!     ACTIVE_EXCIT_LVL (T3 active excitation level)
!
      else if (inp_compare(.false.,test,'t3a_lvl')) then
        if (.not.inp_i(numact))
     1    call errquit('tce_input: no t3a_lvl given',0,INPUT_ERR)
        if (.not.rtdb_put(rtdb,'tce:act_excit_lvl',mt_int,1,numact))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)

      else if (inp_compare(.false.,test,'tcc_spaces')) then
        if (.not.rtdb_put(rtdb,'tce:ltcc',mt_log,1,.true.))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)

!
!     PRINT
!
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'tce')
!
!     END
!
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit('tce_input: unknown directive',0,INPUT_ERR)
      endif
      goto 10
!
! ------
! Return
! ------
!
   20 return
      end
!
!     This is the TCE property input block ("tceprop")
!
      subroutine tce_prop_input(rtdb)
!
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
      logical status
!
      lineresp = .false.
      leftresp = .false.
      quadresp = .false.
!
! -------------------------------------
! What input block are we dealing with?
! -------------------------------------
!
      if (.not.rtdb_cget(rtdb,'tce:module',1,module)) then
        call errquit('tce_prop_input: line ',__LINE__,RTDB_ERR)
      endif
!
! ----------
! Read input
! ----------
!
  100 if (.not. inp_read()) then
        call errquit('tce_prop_input: line ',__LINE__,RTDB_ERR)
      endif
      if (.not. inp_a(test)) then
        call errquit('tce_prop_input: line ',__LINE__,RTDB_ERR)
      endif
!
!     POLARIZABILITY
!
      if (inp_compare(.false.,test,'polarizability').or.
     1    inp_compare(.false.,test,'polar').or.
     2    inp_compare(.false.,test,'alpha')) then
        lineresp = .true.
!
!     HYPERPOLARIZABILITY
!
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
!
!     SECOND HYPERPOLARIZABILITY
!
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
!
!     AFREQ (omega for polarizability)
!
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
  300   continue
!
!     BFREQ (omega for first hyperpolarizability)
!
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
  400   continue
!
!     CFREQ (omega for second hyperpolarizability)
!
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
  500   continue
!
!     END
!
      else if (inp_compare(.false.,test,'end')) then
        goto 200
      else
        call errquit('tce_prop_input: unknown directive',0,INPUT_ERR)
      endif
      goto 100
  200 return
!
! -------------------
! Push values to RTDB
! -------------------
!
      if (.not.rtdb_put(rtdb,'tce:lineresp',mt_log,1,lineresp)) then
        call errquit('tce_prop_input: rtdb_put',__LINE__,RTDB_ERR)
      endif
      if (.not.rtdb_put(rtdb,'tce:leftresp',mt_log,1,leftresp)) then
        call errquit('tce_prop_input: rtdb_put',__LINE__,RTDB_ERR)
      endif
!
! ------
! Return
! ------
!
      end

