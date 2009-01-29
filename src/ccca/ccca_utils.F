      subroutine ccca_setup_basis(rtdb,calctype)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "baslibraryP.fh"
#include "geom.fh"
c
      integer rtdb
c
c     Vars
      character*16 bas_data(28),mybas
      character*80 calctype
      data bas_data/
     $ 'aug-cc-pvdz','aug-cc-pvdz','aug-cc-pvdz','aug-cc-pv(d+d)z',
     $ 'aug-cc-pvtz','aug-cc-pvtz','aug-cc-pvtz','aug-cc-pv(t+d)z',
     $ 'aug-cc-pvqz','aug-cc-pvqz','aug-cc-pvqz','aug-cc-pv(q+d)z',
     $ 'cc-pvtz','cc-pvtz','cc-pvtz','cc-pv(t+d)z',
     $ 'cc-pvtz_dk','cc-pvtz_dk','cc-pvtz_dk','cc-pv(t+d)z_dk',
     $ 'aug-cc-pcvtz','aug-cc-pvtz','aug-cc-pcvtz','aug-cc-pcv(t+d)z',
     $ '6-31G(2df,p)','6-31G(2df,p)','6-31G(2df,p)','6-31G(2df,p)'/
c
      character*255 filename
      double precision cart(3),charge
      character*16 tag,element,bas_name
      character*16 atom_lab(nw_max_atom),symbol
      integer errcode, basis, i
      integer l_q,k_q,atn,ibas
      integer nopen, natoms, nunique, idx, idx2, geom
      logical unique,nodezero,ospherical,osegment,orel
c
c     Externs
c
      logical bas_set_spherical,bas_rtdb_store,bas_add_ucnt_tidy
      logical bas_do_destroy,bas_create,bas_add_ucnt_init,bas_print
      external bas_set_spherical,bas_rtdb_store,bas_add_ucnt_tidy
      external bas_do_destroy,bas_create,bas_add_ucnt_init,bas_print
c
c     Var init
c
      errcode=3000
      nodezero=(ga_nodeid().eq.0)
c
c     Load geometry
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Get number of unique centers
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $   ('ccca_setup_basis: geom_ncent failed', 0, GEOM_ERR)
      if(.not.ma_push_get(mt_int,natoms,'atomnr',l_q,k_q)) call
     &  errquit('ccca_setup_basis alloc k_q failed',911,MA_ERR)
c
      nunique=1
      do 10 idx = 1, natoms
        atom_lab(idx)=' '
        unique=.true.
        if (.not.geom_cent_get(geom,idx,tag,cart,charge)) call
     $     errquit('ccca_setup_basis geom_cent_get',idx,GEOM_ERR)
        if (.not.geom_tag_to_element(tag,symbol,element,atn)) call
     $     errquit('ccca_setup_basis tag_to_element',idx,GEOM_ERR)
        if (idx .gt. 1) then
          do 20 idx2 = 1, idx-1
            if (atom_lab(idx2).eq.symbol) unique=.false.
   20     continue
          if (unique) then
             nunique=nunique+1
             atom_lab(nunique)=tag
             int_mb(k_q+nunique-1)=atn
           endif
        else
          atom_lab(1)=tag   
          int_mb(k_q)=atn
        endif
   10 continue
c
      if(.not.ma_pop_stack(l_q)) call
     &  errquit('ccca_setup_basis dealloc k_q failed',911,MA_ERR)
c
      if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_destroy?', 0, GEOM_ERR)
c
c     setup basis set names
c     the common code should be moved to a function sometime
c
      ospherical=.true.
      osegment=.true.
      orel=.false.
      call bas_set_library_name()
      filename = user_library_name
c
c     Set up basis sets and store to RTDB
c
c     Loop over the seven basis set types
c      Create basis
c      Add element exp and coeff from library
c      Store basis
c      Destroy basis
c
      ibas = 5
c
c     If single point, only first 5 are needed, skip 6th
c     If optimize, also setup 6th basis (6-31G...)
c
      if (calctype.eq.'optimize') ibas=6
c
      do i=0,6 ! loop over the seven basis set classes used
c
      bas_name=bas_data(i*4+1)
      if (.not. bas_create(basis, bas_name)) call
     $  errquit('ccca_setup_basis failed to create basis', 0,BASIS_ERR) 
      if (.not.bas_add_ucnt_init(basis)) call errquit
     &    ('ccca_setup_basis: failed to init add_ucnt ',911,
     &       BASIS_ERR)
c
      do 30 idx=1, nunique
        if (int_mb(k_q+idx-1) .lt. 2) then
          mybas=bas_data(i*4+2)
        elseif (int_mb(k_q+idx-1) .lt. 11) then
          mybas=bas_data(i*4+3)
        else
          mybas=bas_data(i*4+4)
        endif
        call bas_tag_lib(basis, osegment, atom_lab(idx), atom_lab(idx),
     $        mybas, filename, orel)
   30 continue
c
      if (.not. bas_set_spherical(basis, ospherical)) call
     &   errquit ('bas_set_spherical failed ',911, BASIS_ERR)
      if (nodezero) then
        if (.not.bas_print(basis)) call errquit('basprint',0,BASIS_ERR)
      endif
      if (.not. bas_rtdb_store(rtdb, bas_name, basis)) call errquit
     $     ('ccca_setup_basis: failed to store basis', 0, BASIS_ERR)
      if (.not. bas_add_ucnt_tidy(basis))
     $     call errquit('bas_input: tidy failed',0, BASIS_ERR)
      if (.not. bas_do_destroy(basis)) call errquit
     $     ('bas_input: bas_destroy failed', 0, BASIS_ERR)
c
      enddo ! generation of basis types of basis sets loop
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Get the noble gas core, return the number of frozen orbs
      integer function ccca_get_ngcore(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
      integer rtdb
      integer errcode, geom
      integer natoms, idx, frozen
      double precision cart(3),charge
      character*16 tag
c
c     Var init
      errcode=3000
      frozen=0
c
c     Load geometry
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Get number of unique centers
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $   ('ccca_setup_basis: geom_ncent failed', 0, GEOM_ERR)
c
c     Loop through the atoms, accum in frozen
c
      do 10 idx = 1, natoms
        if (.not.geom_cent_get(geom,idx,tag,cart,charge)) call
     $     errquit('ccca_setup_basis geom_cent_get',idx,GEOM_ERR)
        if (charge .lt. 3) then
        elseif (charge .lt. 11) then
          frozen=frozen+1
        elseif (charge .lt. 19) then
          frozen=frozen+5
        elseif (charge .lt. 37) then
          frozen=frozen+9
        elseif (charge .lt. 55) then
          frozen=frozen+18
        elseif (charge .lt. 87) then
          frozen=frozen+27
        else
          frozen=frozen+43
        endif
10    continue
c
      if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_destroy?', 0, GEOM_ERR)
c
c     return the result
c
      ccca_get_ngcore=frozen
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Get the inner noble gas core, return the number of frozen orbs
      integer function ccca_get_ingcore(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
      integer rtdb
      integer errcode, geom
      integer natoms, idx, frozen
      double precision cart(3),charge
      character*16 tag
c
c     Var init
      errcode=4000
      frozen=0
c
c     Load geometry
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Get number of unique centers
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $   ('ccca_setup_basis: geom_ncent failed', 0, GEOM_ERR)
c
c     Loop through the atoms, accum in frozen
c
      do 10 idx = 1, natoms
        if (.not.geom_cent_get(geom,idx,tag,cart,charge)) call
     $     errquit('ccca_setup_basis geom_cent_get',idx,GEOM_ERR)
        if (charge .lt. 11) then
        elseif (charge .lt. 19) then
          frozen=frozen+1
        elseif (charge .lt. 37) then
          frozen=frozen+5
        elseif (charge .lt. 55) then
          frozen=frozen+9
        elseif (charge .lt. 87) then
          frozen=frozen+18
        else
          frozen=frozen+27
        endif
10    continue
c
      if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_destroy?', 0, GEOM_ERR)
c
c     return the result
c
      ccca_get_ingcore=frozen
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This function sets up a standard ccCA calculation and routes
c  the energies and gradients appropriately
c
c  Set oldbasis=srcfile='' to turn off orb projection guess
c
      subroutine ccca_run_calc
     +  (rtdb,oldbasis,srcfile,newbasis,destfile,theory,fcore,name)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
c
c     type the parameters
      integer rtdb, fcore, geom
      character*(*) oldbasis,newbasis,srcfile,destfile,name
      character*2 theory
c
c     Externs needed
      logical task_energy,task_gradient,nw_inp_from_file
      external task_energy,task_gradient,nw_inp_from_file
c
c     Vars
c
      character*256 filename,t2file,newname
      character*8 nwctheory,reftype,tceio,calctype,model,perturbation
      character*18 tag,tag2,vectors_in, frozen_key
      character*80 project(2), vectors_out
      logical result,nodezero,fexist,status
      integer errcode,nopen,natoms
      double precision e_comp
      double precision grad(1536)
c
c     Var init
c
      errcode=5000
      nwctheory='        '
      nodezero=(ga_nodeid().eq.0)
c
c     clean up parts of the rtdb
c
      call rtdb_clean(rtdb, 'tce:')
      call rtdb_clean(rtdb, 'mp2:')
      call rtdb_clean(rtdb, 'ccsd:')
      call rtdb_clean(rtdb, 'moints:')
      call rtdb_clean(rtdb, 'fourindex:')
c
c     what is the correlation model?
c
      if (theory.eq.'mp') then
        if (.not.rtdb_cget(rtdb,'ccca:mp2type',1,nwctheory)) then
          nwctheory='tce'
        endif
      else if (theory.eq.'cc') then
        if (.not.rtdb_cget(rtdb,'ccca:cctype',1,nwctheory)) then
          nwctheory='tce'
        endif
      else
        call errquit(
     +  'ccca_run_calc: unknown correlation model (internal error)',
     +   0,INPUT_ERR)
      endif
c
c     Set the basis set
      call ccca_err(rtdb_cput(rtdb,'ao basis',1,newbasis),errcode) 
c
c     if necessary, make the geometry abelian
c     NOTE: ABELIANIZING/RESTORING SHOULD BE DONE 
c           THE FEWEST NUMBER OF TIMES!!
      if ((nwctheory.eq.'tce').or.(nwctheory.eq.'ccsd(t)')) then
        call ccca_err(rtdb_cput(rtdb,'geometry',1,'geometry_abelian'
     +                      ),errcode)
      endif
c
c     If the user requested a specific ref wfn, set it now
c
      if (rtdb_cget(rtdb,'ccca:refwfn',1,reftype)) then
         if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, reftype)) call
     &      errquit('rtdb_cput in run_calc failed',0,RTDB_ERR)
      endif
c
c     set the input/output vectors
c
      vectors_in='atomic'
      vectors_out=destfile(1:inp_strlen(destfile))
c     if (len(oldbasis).gt.0) then
c       vectors_in='project'
c       project(1)=oldbasis(1:inp_strlen(oldbasis))
c       project(2)=srcfile(1:inp_strlen(srcfile))
c       if (.not. rtdb_cput(rtdb, 'vectors:project', 2, project))
c    $     call errquit('vectors_input: failed to store project',0,
c    &       RTDB_ERR)
c     endif
      if (.not. rtdb_cput(rtdb,'scf:input vectors',1,vectors_in)) call
     $  errquit('vectors_input: failed to store vectors_input', 0,
     &       RTDB_ERR)
      if (.not. rtdb_cput(rtdb,'scf:output vectors',1,vectors_out)) call
     $  errquit('vectors_input: failed to store vectors_output', 0,
     &       RTDB_ERR)
c
c     write the correlation model input
c
      if (nwctheory.eq.'tce') then
        write(frozen_key,'(a,'':'',a)') 'tce','frozen core'
        if (theory.eq.'mp') then
          model='mbpt2'
        else
          model='ccsd'
          perturbation='(t)'
          if (.not.rtdb_cput(rtdb,'tce:perturbative',1,perturbation))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        endif
        if (.not.rtdb_cput(rtdb,'tce:model2e',1,'2espin')) call
     $      errquit('tce_input: model2e',0,RTDB_ERR)
        if (.not.rtdb_cput(rtdb,'tce:model',1,model))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:recompf',mt_log,1,.true.))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:fragment',mt_int,1,-1))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c       if (.not.rtdb_put(rtdb,'tce:tilesize',mt_int,1,10))
c    1    call errquit('tce_input: failed writing to rtdb',0,
c    2    RTDB_ERR)
      else if (nwctheory.eq.'mp2') then
         write(frozen_key,'(a,'':'',a)') 'mp2','frozen core'
      else if (nwctheory.eq.'ccsd(t)') then
         write(frozen_key,'(a,'':'',a)') 'ccsd','frozen core'
      endif
      if (.not. rtdb_cput(rtdb,'task:theory',1,nwctheory)) call
     $    errquit('run_calc: failed writing task_theory',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb, frozen_key, mt_int, 1, fcore))
     $    call errquit('run_calc: failed writing rtbd',0,RTDB_ERR)
c
c     here's where we decide to do a gradient or an energy only
c
      call ccca_err(rtdb_cget(rtdb,'ccca:type',1,calctype),errcode)
      if (calctype.eq.'energy') then
        call ccca_err(task_energy(rtdb),errcode)
      else if (calctype.eq.'gradient') then
        call ccca_err(task_gradient(rtdb),errcode)
      else
        call errquit('ccca_run_calc: unknown calc type requested',
     +                0,INPUT_ERR)
      endif
c
c     if a t2 file exists, get rid of it
c
      if (nodezero) then
        call util_file_name('t2', .false.,.false.,t2file)
        inquire(file=t2file,exist=fexist)
        if (fexist) call util_file_unlink(t2file)
      endif
c
c     get the data
c
      call ccca_err(rtdb_get(rtdb,'task:energy',MT_DBL,1,e_comp)
     +              ,errcode)
      if (calctype .eq. 'gradient') then
c
c     Load geometry
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Get number of unique centers
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $   ('ccca_setup_basis: geom_ncent failed', 0, GEOM_ERR)
c
c     Destroy geometry
c
      if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
c
        call ccca_err(rtdb_get(rtdb,'task:gradient',MT_DBL,
     +                natoms*3,grad),errcode)
        newname=' '
        newname=name(1:len(name))//':gradient'
        call ccca_err(rtdb_put(rtdb,newname(1:len(newname)),
     +                                 MT_DBL,natoms*3,grad),errcode)
      endif
      call ccca_err(rtdb_put(rtdb,name(1:len(name)),MT_DBL,1,e_comp)
     +              ,errcode)
c
c     Reset geometry to point to default one
c
      if ((nwctheory.eq.'tce').or.(nwctheory.eq.'ccsd(t)')) then
        call ccca_err(rtdb_cput(rtdb,'geometry',1,'geometry_fullsym'
     +                      ),errcode)
      endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine checks and potentially re-reads the geometry
c     to force an Abelian symmetry
c
      subroutine ccca_force_abelian(rtdb)
#include "nwc_const.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
#include "geomP.fh"
      integer rtdb,ofile,nops
c
c     Vars
      integer symtab(46)
      data symtab/1,2,3,4,1,4,1,4,1,4,9,1,9,1,9,16,2,16,2,16,
     +            21,2,21,2,21,26,4,26,2,26,26,4,4,4,4,4,
     +            4,3,4,4,21,4,4,21,4,21/
      character*3  newsym
      integer natoms, errcode, idx, geom
      logical nodezero
c
c     Var init
      errcode=6000
      nodezero=(ga_nodeid().eq.0)
c
c     Load geometry
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Store with full symmetry
c
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry_fullsym'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Get number of unique centers
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $   ('ccca_setup_basis: geom_ncent failed', 0, GEOM_ERR)
c
c     Adjust to abelian symmetry
c
      group_number(geom)=symtab(group_number(geom))
      call sym_nwc(geom,rtdb,natoms,.false.,1.0d0,1.0d-2,nops)
c
c     Store with abelian symmetry
c
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry_abelian'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c     Destroy geometry
c
      if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute and print the various ccCA energies
c     (and gradients, eventually)
c
      subroutine ccca_print(rtdb)
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
c     Vars
      integer errcode/70/
      double precision E_D, E_T, E_Q, ErefS4, ErefP, zpe, hthermal,
     +                 CV, CVcorr, DK, DKcorr, mp2vtz, CC, CCcorr, sf
      character*3 variant/'S4'/
      character*80 calctype
      logical parsable, therm
c
c     Var init
      errcode=70
      variant='S4'
      parsable=.false.
      therm=.false.
c
c     Get the values into vars
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/avdz',MT_DBL,1,E_D),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/avtz',MT_DBL,1,E_T),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/avqz',MT_DBL,1,E_Q),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/vtzdk',MT_DBL,1,DK),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc1)/acvtz',MT_DBL,1,CV),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/vtz',MT_DBL,1,mp2vtz),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:ccsd(t)/vtz',MT_DBL,1,CC),
     +              errcode)
c
c     Compute the required values
      ErefS4=-0.577163461538*E_T+1.577163461540*E_Q
      ErefP=0.034867106002*E_D-0.711622433526*E_T+1.676755327520*E_Q
      DKcorr=DK-mp2vtz
      CCcorr=CC-mp2vtz
      CVcorr=CV-E_T
c
c     ...and report them
      if (ga_nodeid().eq.0) then
        write(luout,*) ''
        write(luout,*) 'ccCA-S4 reference energy = ',ErefS4
        write(luout,*) 'ccCA-P reference energy  = ',ErefP
        write(luout,*) 'DK correction            = ',DKcorr
        write(luout,*) 'CCSD(T) correction       = ',CCcorr
        write(luout,*) 'CV correction            = ',CVcorr
        write(luout,*) '---------------------------'
        write(luout,*) 'Total ccCA-S4 energy     = ',
     +              ErefS4+DKcorr+CCcorr+CVcorr
        write(luout,*) 'Total ccCA-P energy      = ',
     +              ErefP+DKcorr+CCcorr+CVcorr
        write(luout,*) ''
      endif
c
c     Check for availability of thermochem
      if (rtdb_get(rtdb,'vib:zpe',MT_DBL,1,zpe)
     +   .and.rtdb_get(rtdb,'vib:hthermal',MT_DBL,1,hthermal)
     +   .and.(ga_nodeid().eq.0)) then
        therm=.true.
        write(luout,*)'Thermochemistry available:'
        write(luout,*)'ccCA-S4 E+ZPE    =',ErefS4+DKcorr+CCcorr+CVcorr+
     &                                     zpe
        write(luout,*)'ccCA-P E+ZPE     =',ErefP+DKcorr+CCcorr+CVcorr+
     &                                     zpe
        write(luout,*)'ccCA-S4 Enthalpy =',
     +             ErefS4+DKcorr+CCcorr+CVcorr+hthermal
        write(luout,*)'ccCA-P Enthalpy  =',
     +             ErefP+DKcorr+CCcorr+CVcorr+hthermal
        write(luout,*)''
      endif
c
c     Now set the proper rtdb entries
      if (rtdb_cget(rtdb,'ccca:variant',1,variant).and. 
     &    ga_nodeid().eq.0) then
        write(luout,*) 'User requested ccCA-'//variant//' energy'
      endif
      if ((variant .eq. 'p') .or. (variant .eq. 'P')) then
        call ccca_err(rtdb_put(rtdb,'task:energy',MT_DBL,1,
     +                ErefP+DKcorr+CCcorr+CVcorr),errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:energy',MT_DBL,1,
     +                ErefP+DKcorr+CCcorr+CVcorr),errcode)
      else if ((variant .eq. 's4') .or. (variant .eq. 'S4')) then
        call ccca_err(rtdb_put(rtdb,'task:energy',MT_DBL,1,
     +                ErefS4+DKcorr+CCcorr+CVcorr),errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:energy',MT_DBL,1,
     +                ErefS4+DKcorr+CCcorr+CVcorr),errcode)
      else
        call errquit('ccca: unknown ccCA variant requested',0,INPUT_ERR)
      endif
      if (ga_nodeid().eq.0) 
     &   write(luout,*) 'Wrote ccCA-'//variant//' energy to the RTDB'
c
c     If we are doing gradients, call the gradient-adapted code
      call ccca_err(rtdb_cget(rtdb,'ccca:type',1,calctype),errcode)
      if (calctype .eq. 'gradient') then
        call ccca_print_gradient(rtdb)
      endif
c
c     check for parsable output
c     load appropriate vars before entering node 0 print block
      if (rtdb_get(rtdb,'ccca:parsable',MT_LOG,1,parsable)) then
      sf=1.0
      if (therm) then
        call ccca_err(rtdb_get(rtdb,'vib:scale',MT_DBL,1,sf),errcode)
      endif
      if (ga_nodeid().eq.0) then
        write(luout,*)'CCCA DATA BLOCK BEGINS'
        write(luout,*)'mp2(fc)/avdz=',E_D
        write(luout,*)'mp2(fc)/avtz=',E_T
        write(luout,*)'mp2(fc)/avqz=',E_Q
        write(luout,*)'mp2(fc)/vtzdk=',DK
        write(luout,*)'mp2(fc1)/acvtz=',CV
        write(luout,*)'mp2(fc)/vtz=',mp2vtz
        write(luout,*)'ccsd(t)/vtz=',CC
        write(luout,*)'Eref(S4)=',ErefS4
        write(luout,*)'Eref(P)=',ErefP
        write(luout,*)'DK correction=',DKcorr
        write(luout,*)'CCSD(T) correction=',CCcorr
        write(luout,*)'CV correction=',CVcorr
        write(luout,*)'Etot(S4)=',ErefS4+DKcorr+CCcorr+CVcorr
        write(luout,*)'Etot(P)=',ErefP+DKcorr+CCcorr+CVcorr
        if (therm) then
          write(luout,*)'Scale factor=',sf
          write(luout,*)'ZPE=',zpe
          write(luout,*)'H therm corr=',hthermal
        endif
        write(luout,*)'CCCA DATA BLOCK ENDS'
        call util_flush(luout)
      endif
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This function extrapolates (and prints, eventually) the gradient
      subroutine ccca_print_gradient(rtdb)
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
c     Vars
      integer errcode/80/, idx
      double precision E_D(1536), E_T(1536), E_Q(1536), grad(1536),
     +                 CV(1536), DK(1536), mp2vtz(1536), CC(1536)
      character*3 variant/'S4'/
      integer natoms
c
c     Var init
      errcode=80
      variant='S4'
c
c     Get the number of atoms
      call ccca_err(rtdb_get(rtdb,'geometry:geometry:ncenter',
     +                MT_INT,1,natoms),errcode)
c
c     Get the gradients
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/avdz:gradient',MT_DBL,
     +              natoms*3,E_D),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/avtz:gradient',MT_DBL,
     +              natoms*3,E_T),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/avqz:gradient',MT_DBL,
     +              natoms*3,E_Q),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/vtzdk:gradient',MT_DBL,
     +              natoms*3,DK),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc1)/acvtz:gradient',
     +              MT_DBL,natoms*3,CV),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:mp2(fc)/vtz:gradient',MT_DBL,
     +              natoms*3,mp2vtz),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:ccsd(t)/vtz:gradient',MT_DBL,
     +              natoms*3,CC),errcode)
c
c     Perform the reference extrapolation
      if (rtdb_cget(rtdb,'ccca:variant',1,variant).and.
     &    ga_nodeid().eq.0) then
        write(luout,*) 'User requested ccCA-'//variant//' energy'
      endif
      if ((variant .eq. 'p') .or. (variant .eq. 'P')) then
        do 10 idx=1, natoms*3
          grad(idx)=0.034867106002*E_D(idx)
     +             -0.711622433526*E_T(idx)
     +             +1.676755327520*E_Q(idx)
   10   continue
      else if ((variant .eq. 's4') .or. (variant .eq. 'S4')) then
        do 20 idx=1, natoms*3
          grad(idx)=-0.577163461538*E_T(idx)+1.577163461540*E_Q(idx)
   20   continue
      else
        call errquit('ccca: Invalid variant requested',0,INPUT_ERR)
      endif
c
c     Perform the gradient corrections to the reference
      do 30 idx=1, natoms*3
        grad(idx)=grad(idx)+DK(idx)+CC(idx)+CV(idx)
     +            -2*mp2vtz(idx)-E_T(idx)
   30 continue
c
c     write the result to the rtdb
      call ccca_err(rtdb_put(rtdb,'task:gradient',MT_DBL,natoms*3,
     +              grad),errcode)
      call ccca_err(rtdb_put(rtdb,'ccca:gradient',MT_DBL,natoms*3,
     +              grad),errcode)
c
c     Inform the user
      if (ga_nodeid().eq.0)
     &   write(luout,*) 'Wrote ccCA-'//variant//' gradient to the RTDB'
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Run the DFT optimization (and thermochem)
c
      subroutine ccca_dft_opt(rtdb)
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
c     Vars
      character*255 filename
      integer nopen,iline,errcode/100/
      logical dotherm/.true./,res
      double precision sf
      logical nodezero
      logical  lcfac(40),lxfac(40),nlcfac(40),nlxfac(40),xccomb(40)
      double precision cfac(40),xfac(40)
c
c     Externs
      logical task_optimize,task_freq,nw_inp_from_file
      external task_optimize,task_freq,nw_inp_from_file
c
c     Var init
      errcode=100
      dotherm=.true.
      nodezero=(ga_nodeid().eq.0)
c
c     Set the basis set
      call ccca_err(rtdb_cput(rtdb,'ao basis',1,'6-31G(2df,p)'),errcode)
c
c     Set up fine grid
      call ccca_err(rtdb_cput(rtdb, 'dft:gridtype', 1, 'fine'),errcode)
c
c     Set up B3LYP functional
      do iline = 1, 40 ! numfunc in DFT input
          xfac(iline)=0.d0
          cfac(iline)=0.d0
         lcfac(iline)=.false.
        nlcfac(iline)=.false.
         lxfac(iline)=.false.
        nlxfac(iline)=.false.
        xccomb(iline)=.false.
      enddo
c
      xccomb(2) = .true.
        cfac(2) = 0.81d0
       lcfac(2) = .true.
      nlcfac(2) = .true.
        cfac(7) = 0.19d0
       lcfac(7) = .true.
      nlcfac(7) = .false.
        xfac(1) = 0.20d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(2) = 0.80d0
       lxfac(2) = .true.
      nlxfac(2) = .false.
        xfac(3) = 0.72d0
       lxfac(3) = .false.
      nlxfac(3) = .true.
c
      call ccca_err(rtdb_put(rtdb,'dft:lcfac',mt_log,40,lcfac),errcode)
      call ccca_err(rtdb_put(rtdb,'dft:nlcfac',mt_log,40,nlcfac),
     &              errcode)
      call ccca_err(rtdb_put(rtdb,'dft:lxfac',mt_log,40,lxfac),errcode)
      call ccca_err(rtdb_put(rtdb,'dft:nlxfac',mt_log,40,nlxfac),
     &              errcode)
      call ccca_err(rtdb_put(rtdb,'dft:xccomb',mt_log,40,xccomb),
     &              errcode)
      call ccca_err(rtdb_put(rtdb,'dft:cfac',mt_dbl,40,cfac),errcode)
      call ccca_err(rtdb_put(rtdb,'dft:xfac',mt_dbl,40,xfac),errcode)
c
c     Do the DFT optimization
      call ccca_err(rtdb_cput(rtdb,'task:theory',1,'dft'),errcode)
      call ccca_err(task_optimize(rtdb),errcode)
c
c     Check to see if we should run thermo
      dotherm=.true.
      if(rtdb_get(rtdb,'ccca:opt+therm',MT_LOG,1,dotherm)) then
        if (dotherm) then 
          if (nodezero) 
     &       write(luout,*)'User requested ccCA thermochemistry'
        else
          if (nodezero)
     &       write(luout,*)'User disabled thermochemistry computation'
        endif 
      endif
      if (dotherm) then
c
c       Check for scale factor override
        if (rtdb_get(rtdb,'vib:scale',MT_DBL,1,sf).and.nodezero) then
          write(luout,*)'ccCA Warning: nonstandard scale factor used(',
     +         sf,')'
        else
          sf=0.9854d0
          call ccca_err(rtdb_put(rtdb,'vib:scale',MT_DBL,1,sf),
     +                  errcode)
        endif
        call ccca_err(task_freq(rtdb),errcode)
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     If the first arg is .f., an error is generated with err code
c     corresponding to the int 2nd arg
c     The errcode is automatically incremented to produce meaningful
c     codes that can be back-traced to a specific statement
c
      subroutine ccca_err(result, code)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      logical result
      integer code
c
      if (.not. result) then
        write (0,*) 'ccca error code=',code
        call errquit('ccca module generated an error',0,INPUT_ERR)
      else
        code=code+1
      endif
c
      return
      end

