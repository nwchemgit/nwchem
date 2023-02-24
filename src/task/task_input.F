C> \ingroup task
C> @{
C>
C> \brief Reads the task directive from input and stores associated 
C> parameters on the RTDB
C>
C> This routine should be called by the root process only.
C> The objective of this routine is to parse the task directive and
C> check for any inconsistencies. The details of the task directive are
C> stored on the RTDB. 
C> ~~~~
C> task [qmmm] <string theory> [<string operation = energy>] [numerical || analytic] 
C>      [ignore]
C> ~~~~
C> The actual task call itself has to be made elsewhere. In case an
C> inconsistency is found the code will abort with an error.
C>
      subroutine task_input(rtdb)
c
c $Id$
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
      integer rtdb !< [Input] The RTDB handle
c
c     Called by process 0 only.  Does nothing but read the input
c     and set corresponding the database parameters
c
c     Called whenever the input parser encounters a task directive.
c     Terminates by calling errquit if there was any error.
c
c     task [qmmm] <string theory> [<string operation = energy>] [numerical || analytic] 
c          [ignore]
c
      logical numerical, ignore, oqmmm, analytic
      logical orism
      logical omm
      logical opes
c      logical oembed
      logical ofragment
      logical osmd
      character*30 test, operation
      character*32 theory
      character*128 err_string
c
      theory    = ' '
      operation = ' '
      analytic  = .false.
      numerical = .false.
      ignore    = .false.
      oqmmm     = .false.
      omm       = .false.
      orism     = .false.
      opes      = .false.
      ofragment = .false.
      osmd      = .false.
c      oembed    = .false.
c
      call rtdb_clean(rtdb, 'task:') ! Delete ALL previous task parameters
      
c
      if (.not. inp_a(test)) goto 1000 ! task keyword
 10   if (inp_a(test)) then
         if (inp_compare(.false.,test,'energy')) then
            operation = 'energy'
         else if (inp_compare(.false.,test,'abc')) then
            operation = 'abc'
         else if (inp_compare(.false.,test,'fep')) then
            operation = 'fep'
         else if (inp_compare(.false.,test,'property')) then
            operation = 'property'
         else if (inp_compare(.false.,test,'optimize')) then
            operation = 'optimize'
         else if (inp_compare(.false.,test,'saddle')) then
            operation = 'saddle'
         else if (inp_compare(.false.,test,'lst')) then
            operation = 'lst'
         else if (inp_compare(.false.,test,'neb')) then
            operation = 'neb'
         else if (inp_compare(.false.,test,'gradient')) then
            operation = 'gradient'
         else if (inp_compare(.false.,test,'hessian')) then
            operation = 'hessian'
         else if (inp_compare(.false.,test,'dynamics')) then
            operation = 'dynamics'
         else if (inp_compare(.false.,test,'qmd')) then
            operation = 'qmd'
         else if (inp_compare(.false.,test,'rfcf')) then
            operation = 'rfcf'
         else if (inp_compare(.false.,test,'thermodynamics')) then
            operation = 'thermodynamics'
         else if (inp_compare(.false.,test,'mepgs')) then
            operation = 'mepgs'
         else if (inp_compare(.false.,test,'string')) then
            operation = 'string'
         else if (inp_compare(.false.,test,'et')) then
            operation = 'et'
         else if (inp_compare(.false.,test,'rt_tddft')) then
            operation = 'rt_tddft'
         else if (inp_compare(.false.,test,'gw')) then
            operation = 'gw'
         else if (inp_compare(.false.,test,'fcidump')) then
            operation = 'fcidump'
c         else if (inp_compare(.false.,test,'prepare')) then
c            operation = 'prepare'
         else if (inp_compare(.false.,test,'freq') .or.
     $            inp_compare(.false.,test,'frequencies')) then
            operation = 'freq'
         else if (inp_compare(.false.,test,'numerical')) then
            numerical = .true.
         else if (inp_compare(.false.,test,'analytic')) then
            analytic  = .true.
         else if (inp_compare(.false.,test,'vscf')) then
            operation  = 'vscf'
         else if (inp_compare(.false.,test,'raman')) then
            operation  = 'raman'
         else if (inp_compare(.false.,test,'dntmc')) then
            operation  = 'dntmc'
         else if (inp_compare(.false.,test,'ccca')) then
            theory  = 'ccca'
         else if (inp_compare(.false.,test,'shell')) then
            theory = 'shell'
            call task_shell_input(rtdb)
         else if (inp_compare(.false.,test,'debug')) then
            if (.not. rtdb_put(rtdb, 'task:debug', MT_LOG, 1, .true.))
     $           call errquit('task_input: failed putting debug flag',0,
     &       RTDB_ERR)
         else if (inp_compare(.false.,test,'ignore')) then
            ignore = .true.
         else if (inp_compare(.false.,test,'rism')) then
            orism = .true.
         else if (inp_compare(.false.,test,'mm')) then
            omm = .true.
         else if (inp_compare(.false.,test,'qmmm')) then
            oqmmm = .true.
c         else if (inp_compare(.false.,test,'embed')) then
c            oembed = .true.
         else if (inp_compare(.false.,test,'fragment')) then
            ofragment = .true.
         else if (inp_compare(.false.,test,'smd')) then
            osmd = .true.
         else if (inp_compare(.false.,test,'pes')) then
            opes = .true.
         else if (inp_compare(.false.,test,
     >                        'steepest_descent')) then
            operation='steepest_descent'
         else if (inp_compare(.false.,test,
     >                        'espcharge')) then
            operation='espcharge'
         else if (inp_compare(.false.,test,
     >                        'Car-Parrinello')) then
            operation='Car-Parrinello'
         else if (inp_compare(.false.,test,
     >                        'Born-Oppenheimer')) then
            operation='Born-Oppenheimer'
         else if (inp_compare(.false.,test,
     >                        'Metropolis')) then
            operation='Metropolis'
         else if (inp_compare(.false.,test,
     >                        'nwpw-optimize')) then
            operation='nwpw-optimize'
         else if (inp_compare(.false.,test,
     >                        'aimd_properties')) then
            operation='aimd_properties'
         else if (inp_compare(.false.,test,
     >                        'exafs')) then
            operation='exafs'
         else if (inp_compare(.false.,test,
     >                        'conjugate_gradient')) then
            operation='conjugate_gradient'
         else if (inp_compare(.false.,test,
     >                        'wavefunction_initializer')) then
            operation='wavefunction_initializer'
         else if (inp_compare(.false.,test,
     >                        'v_wavefunction_initializer')) then
            operation='v_wavefunction_initializer'
         else if (inp_compare(.false.,test,
     >                        'wavefunction_expander')) then
            operation='wavefunction_expander'
         else if (inp_compare(.false.,test,
     >                        'psp_formatter')) then
            operation='psp_formatter'
         else if (inp_compare(.false.,test,
     >                        'psp_generator')) then
            operation='psp_generator'
         else if (inp_compare(.false.,test,
     >                        'cell_optimize')) then
            operation='cell_optimize'
         else if (inp_compare(.false.,test,
     >                        'stress')) then
            operation='stress'
         else if (inp_compare(.false.,test,
     >                        'noit_energy')) then
            operation='noit_energy'
         else if (inp_compare(.false.,test,
     >                        'pspw_dplot')) then
            operation='pspw_dplot'
         else if (inp_compare(.false.,test,
     >                        'band_dplot')) then
            operation='band_dplot'
         else if (inp_compare(.false.,test,
     >                        'wannier')) then
            operation='wannier'
         else if (inp_compare(.false.,test,
     >                        'pspw_et')) then
            operation='pspw_et'
         else if (inp_compare(.false.,test,
     >                        'expand_cell')) then
            operation='expand_cell'
         else if (inp_compare(.false.,test,
     >                        'ionize')) then
            operation='ionize'
         else if (inp_compare(.false.,test,
     >                        'convert_restricted')) then
            operation='convert_restricted'
         else if (inp_compare(.false.,test,
     >                        'structure')) then
            operation='structure'
         else if (inp_compare(.false.,test,
     >                        'dos')) then
            operation='dos'
         else if (inp_compare(.false.,test,
     >                        'dos_dplot')) then
            operation='dos_dplot'
         else if (inp_compare(.false.,test,
     >                        'rdf')) then
            operation='rdf'
         else if (inp_compare(.false.,test,
     >                        'aoints')) then
            operation='aoints'
         else if (inp_compare(.false.,test,
     >                        'grdtrc')) then
            operation='grdtrc'
         else if (inp_compare(.false.,test,
     >                        'nadtrc')) then
            operation='nadtrc'
         else if (inp_compare(.false.,test,
     >                        'lcao')) then
            operation='lcao'
         else if (inp_compare(.false.,test,
     >                        'pauls_matrix_diag')) then
            operation='pauls_matrix_diag'
         else if (inp_compare(.false.,test,
     >                        'translate')) then
            operation='translate'

         else if (theory .eq. ' ') then
            theory = test
            call inp_lcase(theory)
         else
            goto 1000
         endif
         goto 10
      endif
      call inp_clear_err()
c
      if (numerical.and.analytic) then
        write(luout,*)' ********** fatal input error ************'
        write(luout,*)' both numerical and analytic specified '
        goto 1000
      endif
      if (numerical) then
         if (.not. rtdb_put(rtdb, 'task:numerical', mt_log, 1, 
     $        numerical)) call errquit('task_input: rtdb is corrupt',0,
     &       RTDB_ERR)
      endif
      if (analytic) then
         if (.not. rtdb_put(rtdb, 'task:analytic', mt_log, 1, 
     $        analytic)) call errquit('task_input: rtdb is corrupt',0,
     &       RTDB_ERR)
      endif
      if (osmd) then
        if (.not.rtdb_put(rtdb, 'task:smd', mt_log, 1, .true.))
     &       call errquit('task_input:task:smd rtdb is corrupt',911,
     &       RTDB_ERR)
      end if
      if (ofragment) then
        if (.not.rtdb_put(rtdb, 'task:fragment', mt_log, 1, .true.))
     &       call errquit('task_input:task:pes rtdb is corrupt',911,
     &       RTDB_ERR)
      end if
      if (opes) then
        if (.not.rtdb_put(rtdb, 'task:pes', mt_log, 1, .true.))
     &       call errquit('task_input:task:pes rtdb is corrupt',911,
     &       RTDB_ERR)
      end if
      if (omm) then
        if (.not.rtdb_put(rtdb, 'task:mm', mt_log, 1, .true.))
     &       call errquit('task_input:task:mm rtdb is corrupt',911,
     &       RTDB_ERR)
        if (theory.eq.'md') call errquit
     &      ('task_input:theory must be a QM method for mm',911,
     &       RTDB_ERR)
      endif
       if (orism) then
        if (.not.rtdb_put(rtdb, 'task:rism', mt_log, 1, .true.))
     &       call errquit('task_input:task:QMMM rtdb is corrupt',911,
     &       RTDB_ERR)
        if (theory.eq.'md') call errquit
     &      ('task_input:theory must be a QM method for rism',911,
     &       RTDB_ERR)
        if (operation.ne.'energy') call errquit
     &      ("task_input:"//operation//"not supported with QMMM",911,
     &       RTDB_ERR)
      endif
      if (oqmmm) then
        if (.not.rtdb_put(rtdb, 'task:QMMM', mt_log, 1, .true.))
     &       call errquit('task_input:task:QMMM rtdb is corrupt',911,
     &       RTDB_ERR)
        if (theory.eq.'md') call errquit
     &      ('task_input:theory must be a QM method for QMMM',911,
     &       RTDB_ERR)
        if (operation.ne.'energy' .and.
     >      theory.ne.'dplot' .and.
     >      theory.ne.'esp' .and.
     >      operation.ne.'optimize' .and.
     >      operation.ne.'fep' .and.
     >      operation.ne.'neb' .and.
     >      operation.ne.'saddle' .and.
     >      operation.ne.'frequencies' .and.
     >      operation.ne.'hessian' .and.
     >      operation.ne.'freq' .and.
     >      operation.ne.'mepgs' .and.
     >      operation.ne.'string' .and.
     >      operation.ne.'gradient' .and.
     >      operation.ne.'property' .and.
     >      operation.ne.'dynamics' .and.
     >      operation.ne.'qmd' .and.
     >      operation.ne.'rfcf' .and.
     >      operation.ne.'abc' .and.
     >      operation.ne.'et' ) call errquit
     &      ("task_input:"//operation//"not supported with QMMM",911,
     &       RTDB_ERR)
      endif
c
c      if (oembed) then
c        if (.not.rtdb_put(rtdb, 'task:embed', mt_log, 1, .true.))
c     &       call errquit('task_input:task:embed rtdb is corrupt',911,
c     &       RTDB_ERR)
c        if (operation.ne.'energy' .and.
c     >      operation.ne.'optimize' .and.
c     >      operation.ne.'prepare' .and.
c     >      operation.ne.'neb' .and.
c     >      operation.ne.'saddle' .and.
c     >      operation.ne.'frequencies' .and.
c     >      operation.ne.'hessian' .and.
c     >      operation.ne.'freq' .and.
c     >      operation.ne.'gradient' .and.
c     >      operation.ne.'property' .and.
c     >      operation.ne.'dynamics' .and.
c     >      operation.ne.'abc' .and.
c     >      operation.ne.'et' ) call errquit
c     &      ("task_input:"//operation//"not supported with embed",911,
c     &       RTDB_ERR)
c      endif
c
      if (ignore) then
         if (.not. rtdb_put(rtdb, 'task:ignore', mt_log, 1, 
     $        ignore)) call errquit('task_input: rtdb is corrupt',0,
     &       RTDB_ERR)
      endif
c
      if (theory .eq. ' ') then
        if(osmd) theory="smd"
      end if        
      if (theory .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $        call errquit('task_input: rtdb is corrupt', 0,
     &       RTDB_ERR)
      endif
      if (operation .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'task:operation',1,operation))
     $        call errquit('task_input: rtdb is corrupt', 0,
     &       RTDB_ERR)
      endif
c
      return
c
 1000 err_string =             'usage: task <string theory>'
      err_string = err_string//' [<string operation = energy>] '
      err_string = err_string//'[numerical || analytic] [ignore]'
      call errquit(err_string,911, INPUT_ERR)
c
c     THIS IS NEVER EXECUTED.  IT IS HERE BEACUSE -ltask IS AT
c     THE FRONT OF THE LINK PATH AND THE STATE ROUTINES ARE NOT
c     REFERENCED UNTIL LATER.
c
c     JEM: Using a logical from above as dummy return
c      omm = task_save_state(0,'a')
c      omm = task_restore_state(0,'a')
c      omm = task_delete_state(0,'a')
c
      end
C> @}
