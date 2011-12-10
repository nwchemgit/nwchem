      subroutine scf_input(rtdb)
C     $Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "case.fh"
c     
      integer rtdb
c     
      logical status
      character*255 test
      character*16 multiplicity
      double precision conv, shifts(6), tolnr, tol2e
      integer maxiter, poff, nopen, task_chunks, task_bf, maxsub
      integer memsize, filesize
      integer nriscf, ri_store
      logical frozemb
c     
c     Any SCF input results in it being marked non-converged
c     
      status = rtdb_delete(rtdb,'scf:converged')
c     
      multiplicity = ' '         ! If changed it needs storing at end
      nopen = -1                ! If changed it needs storing at end
      frozemb = .false.
c
c     While loop
c     
 10   if (.not. inp_read()) call errquit('scf_input: inp_read failed',0,
     &       INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('scf_input: failed to read keyword', 0,
     &       INPUT_ERR)
c     
      if (inp_compare(.false.,'thresh', test)) then
c     
c     thresh <convergence threshold>
c     
         if (.not. inp_f(conv)) call errquit
     $        ('scf_input: failed reading convergence', 0, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'scf:thresh', mt_dbl, 1, conv))
     $        call errquit('scf-input: rtdb_put failed', 0, RTDB_ERR)
c     
      else if (inp_compare(.false.,'maxiter', test)) then
c     
c     maxiter <maxiter>
c     
         if (.not. inp_i(maxiter)) call errquit
     $        ('scf_input: failed reading maxiter', 0, INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'scf:maxiter', mt_int, 1, maxiter))
     $        call errquit('scf-input: rtdb_put failed', 0, RTDB_ERR)
c     
      else if (inp_compare(.false.,'level', test)) then
c     
c     Level shift 
c     
c     level [pcg initial [tol final]] [nr initial [tol final]]
c     
         call dfill(6, -1.0d0, shifts, 1) ! -1s will be replaced with defaults
 410     if (inp_a(test)) then
            if (inp_compare(.false.,test,'pcg')) then
               poff = 0
            else if (inp_compare(.false.,test,'nr')) then
               poff = 3
            else
               goto 490         ! Unrecognized keyword
            endif
            if (.not. inp_f(shifts(poff+1))) goto 490
            if (.not. inp_f(shifts(poff+2))) then
               if (inp_cur_field().ne.inp_n_field())
     $              call inp_prev_field()
            else
               if (.not. inp_f(shifts(poff+3))) goto 490
            endif
            goto 410
         endif
         if (.not.rtdb_put(rtdb,'scf:level shift info',MT_DBL,6,shifts))
     $        call errquit('scf_input: rtdb_put of shifts failed',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'vectors', test)) then
c     
         call vectors_input(rtdb,'scf')
c     
      else if (inp_compare(.false.,'nr', test)) then
c     
c     nr <switch>
c     
         if (.not. inp_f(tolnr)) call errquit
     $        ('scf_input: nr switch', 0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'scf:full hessian switch',
     $        MT_DBL,1,tolnr))
     $        call errquit('scf_input:rtdb put of tolnr failed', 0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'profile', test)) then
c     
         if (.not. rtdb_put(rtdb,'scf:pstat', MT_LOG, 1, .true.))
     $        call errquit('scf_input:rtdb_put of pstat failed', 0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'print', test) .or.
     $        inp_compare(.false.,'noprint', test)) then
c     
         call util_print_input(rtdb, 'scf')
c     
      else if (inp_compare(.false.,'tol2e', test)) then
c     
c     tol2e <real tol2e>
c     
         if (.not. inp_f(tol2e)) call errquit
     $        ('scf_input: tol2e <thresh>',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'scf:tol2e', MT_DBL, 1, tol2e))
     $        call errquit('scf_input: rtdb_put of tol2e failed',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'nopen', test)) then
c     
c     nopen <integer nopen>
c     
         if (.not. inp_i(nopen)) call errquit
     $        ('scf_input: nopen <nopen>',0, INPUT_ERR)
c     
      else if (inp_compare(.false.,'rhf', test) .or.
     $        inp_compare(.false.,'uhf', test) .or.
     $        inp_compare(.false.,'rohf', test)) then
c     
         call inp_ucase(test)
         if (.not. rtdb_cput(rtdb,'scf:scftype',1,test)) 
     $        call errquit('scf_input: rtdb put of scftype failed',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'singlet', test)) then
         nopen = 0
         multiplicity = test
c
      else if (inp_compare(.false.,'doublet', test)) then
         nopen = 1
         multiplicity = test
c
      else if (inp_compare(.false.,'triplet', test)) then
         nopen = 2
         multiplicity = test
c
      else if (inp_compare(.false.,'quartet', test)) then
         nopen = 3
         multiplicity = test
c
      else if (inp_compare(.false.,'quintet', test)) then
         nopen = 4
         multiplicity = test
c
      else if (inp_compare(.false.,'sextet', test)) then
         nopen = 5
         multiplicity = test
c
      else if (inp_compare(.false.,'septet', test)) then
         nopen = 6
         multiplicity = test
c
      else if (inp_compare(.false.,'octet', test)) then
         nopen = 7
         multiplicity = test
c
      else if (inp_compare(.false.,'frozemb', test)) then
        frozemb = .true.
        if (.not. rtdb_put(rtdb,'scf:frozemb', MT_LOG, 1, frozemb))
     $     call errquit('scf_input: failed to store frozemb',0,
     &     RTDB_ERR)
c     
      else if (inp_compare(.false.,'chunk', test)) then
c     
c     chunk <integer chunksize>
c     
         if (.not. inp_i(task_chunks)) call errquit
     $        ('scf_input: chunk <integer chunk_size>', 0, INPUT_ERR)
         if (task_chunks.lt.1 .or. task_chunks.gt.1024) call errquit
     $        ('scf_input: task_chunks has invalid value',task_chunks,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'fock:task_chunks',
     $        MT_INT,1,task_chunks)) call errquit
     $        ('scf_input: rtdb_put of task_chunks failed',0, RTDB_ERR)
c     
      else if (inp_compare(.false.,'ri-scf', test)) then
c     
c     RI-SCF
c     
         nriscf = 2             ! set defaults
         ri_store = 1
         if (inp_n_field().gt.1) then
 2610       if (inp_a(test)) then
               if (inp_compare(.false., test, 'hessian')) then
                  nriscf = 1
               else if (inp_compare(.false., test, 'full')) then
                  nriscf = 2
               else if (inp_compare(.false., test, 'preconverge')) then
                  nriscf = 3
               else if (inp_compare(.false., test, 'coulomb')) then
                  nriscf = 10
               else if (inp_compare(.false., test, 'auto')) then
                  ri_store = 1
               else if (inp_compare(.false., test, 'memory')) then
                  ri_store = 2
               else if (inp_compare(.false., test, 'disk')) then
                  ri_store = 3
               else
                  call errquit('scf_input: unknown RI option',0,
     &       INPUT_ERR)
               endif
               goto 2610
            endif
         endif
         if (.not. rtdb_put(rtdb,'scf:ri',MT_INT, 1, nriscf))
     $        call errquit('scf_input: failed to store RI option',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb,'scf:ri_store',MT_INT, 1, ri_store))
     $        call errquit('scf_input: failed to store RI option',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'maxsub', test)) then
c     
c     maxsub max_subspace_dimension
c     
         if (.not. inp_i(maxsub)) call errquit
     $        ('scf_input: maxsub max_subspace_dimension',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'scf:maxsub',MT_INT, 1, maxsub))
     $        call errquit('scf_input: failed to store maxsub',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'lagrangian', test)) then
c     
         if (.not. rtdb_put(rtdb,'scf:lagrangian',MT_LOG, 1, .true.))
     $        call errquit('scf_input: failed to store lagrangian',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'taskbf', test)) then
c     
c     taskbf <integer task_bf>
c     
         if (.not. inp_i(task_bf)) call errquit
     $        ('scf_input: taskbf <integer taskbf>', 0, INPUT_ERR)
         if (task_bf.lt.1 .or. task_bf.gt.1024) call errquit
     $        ('scf_input: taskbf has invalid value',task_bf, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'fock:task_bf',MT_INT,1,task_bf))
     $        call errquit('scf_input: rtdb_put of task_bf failed',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'diis', test)) then
c     
         if(.not. rtdb_put(rtdb,'scf:diis', MT_LOG, 1, .true.))
     $        call errquit('scf_input: rtdb_put of diis failed',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'diisbas', test)) then
c     
c     diisbas <maxsub>
c     
         if (.not. inp_i(maxsub)) call errquit
     $        ('scf_input: diisbas max_diis_basis',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'scf:diisbas',MT_INT, 1, maxsub))
     $        call errquit('scf_input: failed to store diisbas',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'sym', test)) then
c     
c     sym [off|on] (defaults to on)
c     
         if (.not. inp_a(test)) test = 'on'
         status = inp_compare(.false.,test,'on')
         if (.not. rtdb_put(rtdb,'scf:skeleton', MT_LOG, 1, status))
     $        call errquit('scf_input: failed to store skeleton',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'adapt', test)) then
c     
c     adapt [off|on] (defaults to on)
c     
         if (.not. inp_a(test)) test = 'on'
         status = inp_compare(.false.,test,'on')
         if (.not. rtdb_put(rtdb,'scf:adapt', MT_LOG, 1, status))
     $        call errquit('scf_input: failed to store adapt',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'lock', test)) then
c     
c     lock [off|on] (defaults to on)
c     
         if (.not. inp_a(test)) test = 'on'
         status = inp_compare(.false.,test,'on')
         if (.not. rtdb_put(rtdb,'scf:lock', MT_LOG, 1, status))
     $        call errquit('scf_input: failed to store lock',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'direct', test)) then
         filesize = -1
         memsize  = -1
         if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, filesize))
     $        call errquit('scf_input: failed to store filesize',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $        call errquit('scf_input: failed to store filesize',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'semidirect', test)) then
c     
c     semidirect [filesize <integer filesize = disksize>] 
c     [memsize  <integer memsize = available>]
c     [filename <string filename>]
c     
         filesize = -33
         memsize  = -33
 3601    if (inp_a(test)) then
            if (inp_compare(.false., test, 'filesize')) then
               if (.not. inp_i(filesize)) goto 3602
            else if (inp_compare(.false., test, 'memsize')) then
               if (.not. inp_i(memsize)) goto 3602
            else if (inp_compare(.false., test, 'filename')) then
               if (.not. inp_a(test)) goto 3602
               if (.not. rtdb_cput(rtdb, 'int2e:filename', 1, test))
     $              call errquit('scf_input: bad rtdb : filename?',0,
     &       RTDB_ERR)
            else
               goto 3602
            endif
            goto 3601
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
     $        call errquit('scf_input: failed to store filesize',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $        call errquit('scf_input: failed to store filesize',0,
     &       RTDB_ERR)
c     
      else if (inp_compare(.false.,'cam', test)) then
c
c       == preserve cam_exch ==
        cam_exch = .true.
        if (.not.rtdb_put(rtdb,'scf:cam_exch', mt_log, 1, cam_exch))
     &   call errquit('scf_input: rtdb_put failed', 0, RTDB_ERR)
        if (.not.inp_f(cam_omega)) cam_omega = 0.5d0
        if (.not.rtdb_put(rtdb,'scf:cam_omega',mt_dbl,1,cam_omega))
     &   call errquit('scf_input: rtdb_put failed cam_omega', 0,
     &     RTDB_ERR)
c
c       == get cam_alpha ==
        if (.not. inp_a(test)) goto 10
        if (inp_compare(.false.,'cam_alpha', test)) then
          if (.not. inp_f(cam_alpha)) cam_alpha = 1.d0
          if (.not.rtdb_put(rtdb,'scf:cam_alpha', mt_dbl, 1, cam_alpha))
     &      call errquit('scf_input: rtdb_put failed cam_alpha', 0,
     &        RTDB_ERR)
        end if
c
c       == get cam_beta ==
        if (.not. inp_a(test)) goto 10
        if (inp_compare(.false.,'cam_beta', test)) then
          if (.not. inp_f(cam_beta)) cam_beta = 1.d0
          if (.not.rtdb_put(rtdb,'scf:cam_beta', mt_dbl, 1, cam_beta))
     &      call errquit('scf_input: rtdb_put failed cam_beta', 0,
     &        RTDB_ERR)
        end if
c
c       == force texas integrals to false for CAM calculations ==
        if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.))
     &   call errquit('cannot set int:cando_txs',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'noscf', test)) then
c
c       == store the noscf flag ==
        if (.not.rtdb_put(rtdb,'scf:noscf', mt_log, 1, .true.))
     &   call errquit('scf_input: rtdb_put failed', 0, RTDB_ERR)
c
      else if (inp_compare(.false.,'end', test)) then
c     
c     end
c     
      else
         call errquit('scf_input: unrecognized directive', 0, INPUT_ERR)
      endif
c     
c     End of while loop over the input
c     
      if (.not. inp_compare(.false.,'end', test)) goto 10
c     
c     Process the multiplicity of multiplicities
c     
      if (nopen .ne. -1) then
         if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, nopen))
     $        call errquit('scf_input: rtdb_put of nopen failed',0,
     &       RTDB_ERR)
         if (multiplicity .ne. ' ') then
            call inp_lcase(multiplicity)
            if (.not. rtdb_cput(rtdb, 'scf:multiplicity', 1, 
     $           multiplicity)) call errquit
     $           ('scf_input: failed to store multipliciity',0,
     &       RTDB_ERR)
         endif
      endif
c     
      return
c     
 490  call errquit('level [pcg initial [tol final]] [nr initial '//
     $     '[tol final]]',0, UNKNOWN_ERR)
 3602 call errquit('semi-direct [filesize <integer filesize'//
     $     ' = disksize>] [memsize <integer memsize = available>]'//
     $     ' [filename <string filename>] ',0, UNKNOWN_ERR)
c     
c     
      end
