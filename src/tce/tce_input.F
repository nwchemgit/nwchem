      subroutine tce_input(rtdb)
c
c $Id: tce_input.F,v 1.9 2002-11-11 18:25:55 sohirata Exp $
c
c Input parser for TCE module for various many-electron theories.
c Also sets default values for input parameters.
c Modified from tddft/tddft_input.F by So Hirata Oct, 2002.
c
c     TCE
c        [FREEZE [[core] (atomic || <integer nfzc default 0>)] \
c                 [virtual <integer nfzv default 0>]]
c        [(LCCD||CCD||CCSD||LCCSD||CCSDT||CCSDTQ|| \ 
c          QCISD||CISD||CISDT||CISDTQ|| \
c          MBPT2||MBPT3||MBPT4) default CCSD]
c        [THRESH <double thresh default 1e-6>]
c        [MAXITER <integer maxiter default 100>]
c        [PRINT (none||low||medium||high||debug)]
c        [IOALGORITHM (0||1||2||3||4) default 2]
c     END
c     
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*20 test
      integer maxiter
      character*10 model
      double precision thresh
      integer ioalg
c
c ------------------
c Set default values
c ------------------
c
c     FREEZE (Frozen cores/virtuals)
c
c     no action is taken
c
c     MODEL (the name of CC model requested)
c
      model='ccsd' 
      if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1  call errquit('tce_input: failed writing to rtdb',0)
c
c     THRESH (convergence threshold for Davidson iteration)
c
      thresh=1.0d-7
      if (.not.rtdb_put(rtdb,'tce:thresh',mt_dbl,1,thresh))
     1  call errquit('tce_input: failed writing to rtdb',0)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      maxiter=100
      if (.not.rtdb_put(rtdb,'tce:maxiter',mt_int,1,maxiter))
     1  call errquit('tce_input: failed writing to rtdb',0)
c
c     IOALGORITHM (I/O method, 0 = Fortran Direct Access, 1 = C Low-Level I/O,
c                              2 = GA Library)
c
      ioalg = 2
      if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1  call errquit('tce_input: failed writing to rtdb',0)
c
c ----------
c Read input
c ----------
c
 10   if (.not. inp_read()) 
     1  call errquit('tce_input: failed reading input',0)
      if (.not. inp_a(test)) 
     1  call errquit('tce_input: failed reading keyword',0)
c
c     FREEZE (Frozen cores/virtuals)
c
      if (inp_compare(.false.,test,'freeze')) then
        call freeze_input(rtdb,'tce')
c
c     MODEL (the name of theory requested)
c
      else if (inp_compare(.false.,test,'ccd')) then
        model='ccd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'lccd')) then
        model='lccd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'ccsd')) then
        model='ccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'lccsd')) then
        model='lccsd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'ccsdt')) then
        model='ccsdt'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'ccsdtq')) then
        model='ccsdtq'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'qcisd')) then
        model='qcisd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'cisd')) then
        model='cisd'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'cisdt')) then
        model='cisdt'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'cisdtq')) then
        model='cisdtq'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'mbpt2')) then
        model='mbpt2'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'mbpt3')) then
        model='mbpt3'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
      else if (inp_compare(.false.,test,'mbpt4')) then
        model='mbpt4'
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0)
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
     1    call errquit('tce_input: failed writing to rtdb',0)
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
     1    call errquit('tce_input: failed writing to rtdb',0)
c
c     IOALGORITHM (I/O method)
c
      else if (inp_compare(.false.,test,'ioalgorithm')) then
        if (.not.inp_i(ioalg)) then
          write(LuOut,*) 'tce_input: ioalgorithm value not found; ',
     1      'default value of 2 will be used'
          ioalg=2
        endif
        if (.not.rtdb_put(rtdb,'tce:ioalg',mt_int,1,ioalg))
     1    call errquit('tce_input: failed writing to rtdb',0)
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
        call errquit('tce_input: unknown directive',0)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
