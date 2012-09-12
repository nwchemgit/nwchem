C> \ingroup task
C> @{
      logical function task_ccca(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb
c
c     This is the ccCA main block
c
c     Variables follow
c
      character*256 movecsfile(2,6) 
      character*80 calctype
      character*10 opttype, geom_abel
      integer stage, geom, errcode
      logical tm, nodezero, status
c
c     Var init
      errcode=0
      stage=0
c     geom_abel='abel'
      nodezero=(ga_nodeid().eq.0)
c
c     set the file names (var init)
c
      call util_file_prefix('aVDZ.movecs',movecsfile(1,1))
      call util_file_prefix('aVTZ.movecs',movecsfile(1,2))
      call util_file_prefix('aVQZ.movecs',movecsfile(1,3))
      call util_file_prefix('VTZ.movecs',movecsfile(1,4))
      call util_file_prefix('VTZdk.movecs',movecsfile(1,5))
      call util_file_prefix('aCVTZ.movecs',movecsfile(1,6))
      call util_file_prefix('aVDZdk.movecs',movecsfile(2,1))
      call util_file_prefix('aVTZdk.movecs',movecsfile(2,2))
      call util_file_prefix('aVQZdk.movecs',movecsfile(2,3))
      call util_file_prefix('VTZ.movecs',movecsfile(2,4))
      call util_file_prefix('VTZdk.movecs',movecsfile(2,5))
      call util_file_prefix('aCVDZdk.movecs',movecsfile(2,6))
      if (.not.rtdb_cput(rtdb,'ccca:movecsfile',12,movecsfile))
     +   call errquit('task_ccca: rtdb_cput failed',0,RTDB_ERR)
c
c     Check if energy or optimization is to be done
c
      if (.not.rtdb_cget(rtdb,'ccca:type',1,calctype)) 
     +   call ccca_input_default(rtdb) 
c
      if (calctype.eq.'optimize') then
c
c     The ccCA module has two opt methods
c     1) DFT opt per ccCA protocol (B3LYP/cc-pVTZ)
c     2) direct opt on the ccCA PES
c
        if (nodezero) then
        write(luout,*) 'ccCA Geometry Optimization Procedure'
        write(luout,*) '------------------------------------'
        write(luout,*) ' ' 
        endif
        if (rtdb_cget(rtdb,'ccca:opttype',1,opttype)) then
          if (opttype.eq.'direct') then
            if (nodezero) then 
              write(luout,*) 'Optimization on direct ccCA PES requested'
              write(luout,*) ' ' 
              call errquit('task_ccca: direct optimization not ',
     +                     'available.',0,RTDB_ERR)
            endif
          else if (opttype.eq.'dft') then
            if (nodezero)
     +      write(luout,*) 'DFT optimization conforming to ccCA ',
     +                     'requested'
          else 
            call errquit('task_ccca: optimization method not ',
     +                     'recognized.',0,RTDB_ERR)
          endif
        else
          if (nodezero)
     $    write(luout,*) 'Default ccCA-DFT optimization will be ',
     +                   'performed'
          opttype='dft'
        endif
      else 
        if (nodezero) then
        write(luout,*) 'ccCA Energy Calculation Procedure'
        write(luout,*) '---------------------------------'
        write(luout,*) ' ' 
        endif
      endif
      call util_flush(luout)
c
c     check for restart status 
c     restart is only available for energy calculations 
c
      if (rtdb_get(rtdb,'ccca:stage',MT_INT,1,stage)) then
        if (stage.gt.0) then
          call ccca_err(rtdb_put(rtdb,'ccca:restart',MT_LOG,1,.true.),
     &                  errcode)
          if (nodezero)
     +    write(luout,*)'ccCA restart detected: proceding from stage ',
     +                  stage
        endif
      endif
c
c     Obtain the noble gas frozen core settings
c     Check if there is any transition metal
c
      call ccca_get_ngcore(rtdb,tm)
      call ccca_get_ingcore(rtdb)
c
c     Setup the basis sets
c
      call ccca_setup_basis(rtdb,tm)
c
c     Default single point energy calcuation
c
      call ccca_err(rtdb_cput(rtdb,'ccca:intype',1,'energy'), errcode)
c
c     For ccCA optimizations, call the DFT setup for default optimization
c     No DKH Hamiltonian for DFT optimization
c
      if ((calctype.eq.'optimize').and.(opttype.eq.'dft')) then 
        if (stage.eq.0) then
          if (nodezero) then
            write(luout,*) 
     $      'ccCA: DFT optimization with cc-pVTZ basis'
          endif
          call ccca_dft_opt(rtdb)
        endif
      endif
c
      call ccca_energy(rtdb)
c
c     If we made it this far, we can return successful
      task_ccca=.true.
      call ccca_err(rtdb_put(rtdb,'task:status',MT_LOG,1,.true.),
     +              errcode)
      call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,0),errcode)
c
c     nwchem freezes after ccca: print end of run flag
      if (nodezero) then
        write(luout,*) 'Leaving ccCA module...'
        call util_flush(luout)
      endif
      return
      end
C> @}

c $Id$
