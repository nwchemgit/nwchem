      subroutine scf_input(rtdb)
C $Id: scf_input.F,v 1.6 1996-10-04 21:48:22 d3g681 Exp $
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      logical status
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 35)
      character*10 dirs(num_dirs)
      character*255 test
      character*5 scftype(1)
      double precision conv, shifts(6), tolnr, tol2e
      integer maxiter, poff, nopen, task_chunks, task_bf, maxsub
      integer memsize, filesize
      integer nriscf, ri_store
      data dirs /'thresh','maxiter','scftype','level','vectors',
     $     'nr', 'profile', 'print', 'noprint', 'tol2e',
     $     'noadapt', 'nopen', 'uhf', 'rohf', 'rhf', 
     $     'singlet', 'doublet', 'triplet', 'quartet', 'quintet',
     $     'chunk', 'nosym', 'ri', 'maxsub', 'lagrang', 'taskbf',
     $     'diis', 'diisbas', 'sym', 'adapt', 'lagrangian', 'lock',
     $     'direct', 'semidirect', 'end'/
c
c     For now, any SCF input results in it being marked non-converged
c
      status = rtdb_delete(rtdb,'scf:converged')
c
 10   if (.not. inp_read()) call errquit
     $     ('scf_input: inp_read failed', 0)
      if (.not. inp_a(test))
     $     call errquit('scf_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('scf_input: unknown directive', 0)
c     
      goto (100,  200,  300,  400,  500,  600,  700,  800,  900, 
     $     1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
     $     1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700,
     $     2800, 2900, 3000, 2500, 3100, 3200, 3300,
     $     9999) ind
      call errquit('scf_input: unimplemented directive', ind)
c     
c     thresh <convergence threshold>
c     
 100  if (.not. inp_f(conv)) call errquit
     $     ('scf_input: failed reading convergence', 0)
      if (.not. rtdb_put(rtdb, 'scf:thresh', mt_dbl, 1, conv))
     $     call errquit('scf-input: rtdb_put failed', 0)
      goto 10
c
c     Maximum iterations
c
 200  if (.not. inp_i(maxiter)) call errquit
     $     ('scf_input: failed reading maxiter', 0)
      if (.not. rtdb_put(rtdb, 'scf:maxiter', mt_int, 1, maxiter))
     $     call errquit('scf-input: rtdb_put failed', 0)
      goto 10
c     
c     Alternate setting for RHF, ROHF, UHF wavefunction selection
c     
 300  if (.not. inp_a(scftype)) call errquit
     $     ('scf_input: failed reading wavefunction type', 0)
      call inp_ucase(scftype)
      if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, scftype ))
     $     call errquit('scf-input: rtdb_put failed', 0)
      goto 10
c
c     Level shift 
c 
c     level [pcg initial [tol final]] [nr initial [tol final]]
c
 400  call dfill(6, -1.0d0, shifts, 1)  ! -1s will be replaced with defaults
 410  if (.not. inp_a(test)) goto 450
      if (inp_compare(.false.,test,'pcg')) then
         poff = 0
      else if (inp_compare(.false.,test,'nr')) then
         poff = 3
      else
         goto 490  ! Unrecognized keyword
      endif
      if (.not. inp_f(shifts(poff+1))) goto 490
      if (.not. inp_f(shifts(poff+2))) then
         if (inp_cur_field().ne.inp_n_field())
     $        call inp_prev_field()
      else
         if (.not. inp_f(shifts(poff+3))) goto 490
      endif
      goto 410
 450  if (.not. rtdb_put(rtdb,'scf:level shift info',MT_DBL,6,shifts))
     $     call errquit('scf_input: rtdb_put of shifts failed',0)
      goto 10
c
 490  call errquit('level [pcg initial [tol final]] [nr initial '//
     $     '[tol final]]',0)
c
c     vectors
c
 500  call vectors_input(rtdb,'scf')
      goto 10
c
c     nr switch
c
 600  if (.not. inp_f(tolnr)) call errquit
     $     ('scf_input: nr switch', 0)
      if (.not. rtdb_put(rtdb,'scf:full hessian switch',MT_DBL,1,tolnr))
     $     call errquit('scf_input:rtdb put of tolnr failed', 0)
      goto 10
c
c     enable profiling
c
 700  if (.not. rtdb_put(rtdb,'scf:pstat', MT_LOG, 1, .true.))
     $     call errquit('scf_input:rtdb_put of pstat failed', 0)
      goto 10
c
c     print and noprint
c
 800  continue
 900  continue
      call util_print_input(rtdb, 'scf')
      goto 10
c
c     tol2e <real tol2e>
c
 1000 if (.not. inp_f(tol2e)) call errquit
     $     ('scf_input: tol2e <thresh>',0)
      if (.not. rtdb_put(rtdb,'scf:tol2e', MT_DBL, 1, tol2e))
     $     call errquit('scf_input: rtdb_put of tol2e failed',0)
      goto 10
c
c     noadapt ... OBSOLETE
c
 1100 write(6,*) ' NOADAPT is superceded by ADAPT [ON|OFF]'
      if(.not. rtdb_put(rtdb,'scf:adapt', MT_LOG, 1, .false.))
     $     call errquit('scf_input: rtdb_put of adapt failed',0)
      goto 10
c
c     nopen <integer nopen>
c
 1200 if (.not. inp_i(nopen)) call errquit('scf_input: nopen <nopen>',0)
      if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, nopen))
     $     call errquit('scf_input: rtdb_put of nopen failed',0)
      goto 10
c
c     RHF/ROHF/UHF
c
 1300 continue
 1400 continue
 1500 continue
      call inp_ucase(test)
      if (.not. rtdb_cput(rtdb,'scf:scftype',1,test)) 
     $     call errquit('scf_input: rtdb put of scftype failed',0)
      goto 10
c
c     singlet/doublet/triplet/quartet/quintet
c
 1600 nopen = 0
      goto 2001
 1700 nopen = 1
      goto 2001
 1800 nopen = 2
      goto 2001
 1900 nopen = 3
      goto 2001
 2000 nopen = 4
 2001 if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, nopen))
     $     call errquit('scf_input: rtdb_put of nopen failed',0)
      call inp_lcase(test)
      if (.not. rtdb_cput(rtdb,'scf:multiplicity', 1, test))
     $     call errquit('scf_input: rtdb_put of multi?',0)
      goto 10
c
c     chunk <integer chunksize>
c
 2100 if (.not. inp_i(task_chunks)) call errquit
     $     ('scf_input: chunk <integer chunk_size>', 0)
      if (task_chunks.lt.1 .or. task_chunks.gt.1024) call errquit
     $     ('scf_input: task_chunks has invalid value',task_chunks)
      if (.not. rtdb_put(rtdb,'fock:task_chunks',MT_INT,1,task_chunks))
     $     call errquit('scf_input: rtdb_put of task_chunks failed',0)
      goto 10
c
c     nosym ... OBSOLETE
c
 2200 write(6,*) ' NOSYM is superceded by SYM [ON|OFF]'
      if (.not. rtdb_put(rtdb,'scf:skeleton',MT_LOG, 1, .false.))
     $     call errquit('scf_input: failed to store skeleton',0)
      goto 10
c
c     RI-SCF
c
 2300 continue
c     set defaults
      nriscf = 2
      ri_store = 1
      if (inp_n_field().gt.1) then
 2310   if (inp_a(test)) then
          if (inp_compare(.false., test, 'hessian')) then
            nriscf = 1
          else if (inp_compare(.false., test, 'full')) then
            nriscf = 2
          else if (inp_compare(.false., test, 'preconverge')) then
            nriscf = 3
          else if (inp_compare(.false., test, 'auto')) then
            ri_store = 1
          else if (inp_compare(.false., test, 'memory')) then
            ri_store = 2
          else if (inp_compare(.false., test, 'disk')) then
            ri_store = 3
          else
            call errquit('scf_input: unknown RI option',0)
          endif
          goto 2310
        endif
      endif
      if (.not. rtdb_put(rtdb,'scf:ri',MT_INT, 1, nriscf))
     $     call errquit('scf_input: failed to store RI option',0)
      if (.not. rtdb_put(rtdb,'scf:ri_store',MT_INT, 1, ri_store))
     $     call errquit('scf_input: failed to store RI option',0)

        goto 10
c
c     maxsub max_subspace_dimension
c
 2400   if (.not. inp_i(maxsub)) call errquit
     $       ('scf_input: maxsub max_subspace_dimension',0)
        if (.not. rtdb_put(rtdb,'scf:maxsub',MT_INT, 1, maxsub))
     $     call errquit('scf_input: failed to store maxsub',0)
        goto 10
c
c     lagrang or lagrangian dump lagrangian
c
 2500 if (.not. rtdb_put(rtdb,'scf:lagrangian',MT_LOG, 1, .true.))
     $     call errquit('scf_input: failed to store lagrangian',0)
      goto 10
c
c     taskbf <integer task_bf>
c
 2600 if (.not. inp_i(task_bf)) call errquit
     $     ('scf_input: taskbf <integer taskbf>', 0)
      if (task_bf.lt.1 .or. task_bf.gt.1024) call errquit
     $     ('scf_input: taskbf has invalid value',task_bf)
      if (.not. rtdb_put(rtdb,'fock:task_bf',MT_INT,1,task_bf))
     $     call errquit('scf_input: rtdb_put of task_bf failed',0)
      goto 10
c
c     toggle DIIS convergence
c
 2700 if(.not. rtdb_put(rtdb,'scf:diis', MT_LOG, 1, .true.))
     $     call errquit('scf_input: rtdb_put of diis failed',0)
      goto 10
c
c     max DIIS basis
c
 2800   if (.not. inp_i(maxsub)) call errquit
     $       ('scf_input: diisbas max_diis_basis',0)
        if (.not. rtdb_put(rtdb,'scf:diisbas',MT_INT, 1, maxsub))
     $     call errquit('scf_input: failed to store diisbas',0)
        goto 10
c
c     sym [off|on] (defaults to on)
c
 2900   if (.not. inp_a(test)) test = 'on'
        status = inp_compare(.false.,test,'on')
        if (.not. rtdb_put(rtdb,'scf:skeleton', MT_LOG, 1, status))
     $       call errquit('scf_input: failed to store skeleton',0)
        goto 10
c
c     adapt [off|on] (defaults to on)
c
 3000   if (.not. inp_a(test)) test = 'on'
        status = inp_compare(.false.,test,'on')
        if (.not. rtdb_put(rtdb,'scf:adapt', MT_LOG, 1, status))
     $       call errquit('scf_input: failed to store adapt',0)
        goto 10
c
c     lock [off|on] (defaults to on)
c
 3100   if (.not. inp_a(test)) test = 'on'
        status = inp_compare(.false.,test,'on')
        if (.not. rtdb_put(rtdb,'scf:lock', MT_LOG, 1, status))
     $       call errquit('scf_input: failed to store lock',0)
        goto 10
c
c     direct
c
 3200   filesize = -1
        memsize  = -1
        if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, filesize))
     $       call errquit('scf_input: failed to store filesize',0)
        if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $       call errquit('scf_input: failed to store filesize',0)
        goto 10
c
c     semidirect [filesize <integer filesize = disksize>] 
c                [memsize  <integer memsize = available>]
c                [filename <string filename>]
c
 3300   filesize = -33
        memsize  = -33
 3301   if (inp_a(test)) then
           if (inp_compare(.false., test, 'filesize')) then
              if (.not. inp_i(filesize)) goto 3302
           else if (inp_compare(.false., test, 'memsize')) then
              if (.not. inp_i(memsize)) goto 3302
           else if (inp_compare(.false., test, 'filename')) then
              if (.not. inp_a(test)) goto 3302
              if (.not. rtdb_cput(rtdb, 'int2e:filename', 1, test))
     $             call errquit('scf_input: bad rtdb : filename?',0)
           else
              goto 3302
           endif
           goto 3301
        endif
        call inp_clear_err()
c
c     User should see a value of zero meaning zero not use default
c
        if (filesize.eq.0) filesize = -1
        if (memsize.eq.0) memsize = -1
        if (filesize.eq.-33) filesize=0
        if (memsize.eq.-33) memsize=0
        if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, filesize))
     $       call errquit('scf_input: failed to store filesize',0)
        if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $       call errquit('scf_input: failed to store filesize',0)
        goto 10
c
 3302   call errquit('semi-direct [filesize <integer filesize'//
     $       ' = disksize>] [memsize <integer memsize = available>]'//
     $       ' [filename <string filename>] ',0)
c
c     end
c     
 9999 continue
c     
      end
