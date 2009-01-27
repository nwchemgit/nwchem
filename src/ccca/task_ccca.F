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
      integer ngcore, ingcore, natoms, errcode
      character*256 avdz_file, avtz_file, avqz_file,
     +              vtz_file, vtzdk_file, acvtz_file
      double precision e_comp
      character*80 proj_vecs(2), calctype
      character*3 mp2type/'tce'/
      character*10 cctype,opttype
      double precision grad(1536)
      integer stage, geom
      logical nodezero,status
c
c     Externals needed
c
      integer ccca_get_ngcore,ccca_get_ingcore
      external ccca_get_ngcore,ccca_get_ingcore
      logical ccca_sane
      external ccca_sane
c
c     Var init
      errcode=0
      stage=0
      nodezero=(ga_nodeid().eq.0)
c
c     Check if energy or optimization is done
c
      if (.not.rtdb_cget(rtdb,'ccca:type',1,calctype))
     +   call errquit('task_ccca: cannot get calc type',0,RTDB_ERR)
c
      if (calctype.eq.'optimize') then
c
c     The ccCA module has two opt methods
c     1) DFT opt per ccCA protocol
c     2) direct opt on the ccCA PES
c
       
        if (nodezero) then
        write(luout,*) 'ccCA Geometry Optimization Procedure'
        write(luout,*) '------------------------------------'
        write(luout,*) ' ' 
        endif
        if (rtdb_cget(rtdb,'ccca:opttype',1,opttype)) then
          if (opttype.eq.'direct'.and.nodezero) then
            write(luout,*) 'Optimization on direct ccCA PES requested'
          else if (opttype.eq.'dft'.and.nodezero) then
            write(luout,*) 'DFT optimization conforming to ccCA ',
     +                     'requested'
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
c
      if (rtdb_get(rtdb,'ccca:stage',MT_INT,1,stage)) then
        if ((stage.gt.0).and.nodezero) then
          write(6,*)'ccCA restart detected: proceding from stage ',stage
        endif
      endif
c
c     set the file names (var init)
c
      call util_file_prefix('aVDZ.movecs',avdz_file)
      call util_file_prefix('aVTZ.movecs',avtz_file)
      call util_file_prefix('aVQZ.movecs',avqz_file)
      call util_file_prefix('VTZ.movecs',vtz_file)
      call util_file_prefix('VTZdk.movecs',vtzdk_file)
      call util_file_prefix('aCVTZ.movecs',acvtz_file)
c
c     Obtain the noble gas frozen core settings
c
      ngcore=ccca_get_ngcore(rtdb)
      ingcore=ccca_get_ingcore(rtdb)
c
c     Setup the basis sets
c
      call ccca_setup_basis(rtdb,calctype)
c
c     For ccCA optimizations, call the DFT setup
c
      if (calctype.eq.'optimize') then 
        if (stage.eq.0) then
           if (ga_nodeid().eq.0) then
              write(luout,*) 
     $        'CCCA: DFT optimization with 6-31G(2df,p) basis'
           endif
          call ccca_dft_opt(rtdb)
        endif
        calctype='energy'
        call ccca_err(rtdb_cput(rtdb,'ccca:type',1,'energy'),errcode)
      endif
c
c     Setup the abelian geometry after optimization
c
      call ccca_force_abelian(rtdb)
c
c     The first calc has no orb proj guess, so it's special
c     stage 1: MP2(fc)/aVDZ
      if (.not.(stage.gt.1)) then
        if (ga_nodeid().eq.0) then
           write(luout,*) 
     $     'CCCA: MP2 calculation with aug-cc-pVDZ basis'
        endif
        call ccca_run_calc(rtdb,'','','aug-cc-pvdz',avdz_file,
     +                     'mp',ngcore,'ccca:mp2(fc)/avdz')
        stage=1
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage)
     +                ,errcode)
      endif
c
c     stage 2: MP2(fc)/aVTZ
      if (stage.lt.2) then
        if (ga_nodeid().eq.0) then
           write(luout,*) 
     $     'CCCA: MP2 calculation with aug-cc-pVTZ basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvdz',avdz_file,
     +     'aug-cc-pvtz',avtz_file,'mp',ngcore,'ccca:mp2(fc)/avtz')
        stage=2
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
      endif
c
c     stage 3: MP2(fc)/aVQZ
      if (stage.lt.3) then
        if (ga_nodeid().eq.0) then
           write(luout,*) 
     $     'CCCA: MP2 calculation with aug-cc-pVQZ basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvtz',avtz_file,
     +     'aug-cc-pvqz',avqz_file,'mp',ngcore,'ccca:mp2(fc)/avqz')
        stage=3
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
      endif
c
c     Currently, I cannot find a reliable way of the getting
c     the mp2 component of the CCSD(T)/VTZ calc
c     So it is run separately below
c
c     stage 4: MP2(fc)/VTZ
      if (stage.lt.4) then
        if (ga_nodeid().eq.0) then
           write(luout,*) 
     $     'CCCA: MP2 calculation with cc-pVTZ basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvdz',avdz_file,
     +     'cc-pvtz',vtz_file,'mp',ngcore,'ccca:mp2(fc)/vtz')
        stage=4
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
      endif
c
c     Need to enable DK integrals for this one, then switch off after
c     stage 5: DK-MP2(fc)/VTZdk
      call ccca_err(rtdb_put(rtdb,'doug_kroll',MT_LOG,1,.true.),errcode)
      call ccca_err(rtdb_put(rtdb,'doug_kroll:type',MT_INT,1,2),errcode)
      if (stage.lt.5) then
        if (ga_nodeid().eq.0) then
           write(luout,*) 
     $     'CCCA: Relativity through Douglas-Kroll MP2 calculation'
           write(luout,*) '      with cc-pVTZ basis'
        endif
        call ccca_run_calc(rtdb,'cc-pvtz',vtz_file,
     +       'cc-pvtz_dk',vtzdk_file,'mp',ngcore,'ccca:mp2(fc)/vtzdk')
        stage=5
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
      endif
      call ccca_err(rtdb_put(rtdb,'doug_kroll',MT_LOG,1,.false.),
     +              errcode)
c
c     Change to the ING core and perform the CV calc
c     stage 6: MP2(fc1)/aCVTZ
      if (stage.lt.6) then
        if (ga_nodeid().eq.0) then
           write(luout,*) 
     $     'CCCA: Calculation of core-valence correction at MP2 level'
           write(luout,*) '      with aug-cc-pcvtz basis'
         endif
        call ccca_run_calc(rtdb,'aug-cc-pvtz',avtz_file,
     +    'aug-cc-pcvtz',acvtz_file,'mp',ingcore,'ccca:mp2(fc1)/acvtz')
        stage=6
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
      endif
c
c     Perform CCSD(T)
c     stage 7: CCSD(T)/VTZ
      if (stage.lt.7) then
        if (ga_nodeid().eq.0) then
           write(luout,*)
     $     'CCCA: CCSD(T) calculation with cc-pVTZ basis'
        endif
        call ccca_run_calc(rtdb,'cc-pvtz',vtz_file,
     +           'cc-pvtz',vtz_file,'cc',ngcore,'ccca:ccsd(t)/vtz')
        stage=7
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
      endif
c
c     We have completed all computational stages
      stage=8
      call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),errcode)
c
c     Restore the level of theory
      call ccca_err(rtdb_cput(rtdb,'task:theory',1,'ccca'),errcode)
c
c     Print the results
      call ccca_print(rtdb)
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

