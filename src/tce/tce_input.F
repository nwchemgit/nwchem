      subroutine tce_input(rtdb)
c
c $Id: tce_input.F,v 1.16 2003-10-19 00:52:07 sohirata Exp $
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
c          MBPT2||MBPT3||MBPT4||MP2||MP3||MP4) default CCSD]
c        [THRESH <double thresh default 1e-6>]
c        [MAXITER <integer maxiter default 100>]
c        [PRINT (none||low||medium||high||debug)]
c        [IO (fortran||eaf||ga||sf||replicated||dra||ga_eaf) default ga]
c        [DIIS <integer diis default 5>]
c        [NROOTS <integer nroots default 0>]
c        [TARGET <integer target default 1>]
c        [TARGETSYM <character targetsym default 'none'>]
c        [SYMMETRY]
c        [DIPOLE]
c        [TILESIZE <no default (automatically adjusted)>]
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
      character*10 module
      double precision thresh
      character*10 ioalgchar
      integer ioalg
      integer reference
      integer diis
      integer nroots
      integer target
      integer tilesize
      character*4 targetsym
      logical symmetry
      logical left
      character*10 perturbative
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
c
c     NROOTS (the number of excited state roots)
c
      nroots=0
      if (.not.rtdb_put(rtdb,'tce:nroots',mt_int,1,nroots))
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
c     DIPOLE (dipole moments & dipole transition moments)
c
      left=.false.
      if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
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
c     MODEL (the name of theory requested)
c
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
c     SYMMETRY (restricts the roots to have the TARGETSYM irrep)
c
      else if (inp_compare(.false.,test,'dipole')) then
        left=.true.
        if (.not.rtdb_put(rtdb,'tce:left',mt_log,1,left))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
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
