      subroutine ccca_setup_basis(rtdb,tm)
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
      logical tm
c
c     Vars
      character*18 bas_data(48), mybas, bas_name
      data bas_data/
     $ 'aug-cc-pvdz','aug-cc-pvdz','aug-cc-pvdz','aug-cc-pv(d+d)z',
     $ 'aug-cc-pvtz','aug-cc-pvtz','aug-cc-pvtz','aug-cc-pv(t+d)z',
     $ 'aug-cc-pvqz','aug-cc-pvqz','aug-cc-pvqz','aug-cc-pv(q+d)z',
     $ 'cc-pvtz','cc-pvtz','cc-pvtz','cc-pv(t+d)z',
     $ 'cc-pvtz-dk','cc-pvtz-dk','cc-pvtz-dk','cc-pv(t+d)z-dk',
     $ 'aug-cc-pcvtz','aug-cc-pvtz','aug-cc-pcvtz','aug-cc-pcv(t+d)z',
     $ 'aug-cc-pvdz-dk','aug-cc-pvdz-dk',
     $ 'aug-cc-pvdz-dk','aug-cc-pvdz-dk',
     $ 'aug-cc-pvtz-dk','aug-cc-pvtz-dk',
     $ 'aug-cc-pvtz-dk','aug-cc-pv(t+d)z-dk',
     $ 'aug-cc-pvqz-dk','aug-cc-pvqz-dk',
     $ 'aug-cc-pvqz-dk','aug-cc-pv(q+d)z-dk',
     $ 'cc-pvtz','cc-pvtz','cc-pvtz','cc-pv(t+d)z',
     $ 'cc-pvtz-dk','cc-pvtz-dk','cc-pvtz-dk','cc-pv(t+d)z-dk',
     $ 'aug-cc-pcvdz-dk','aug-cc-pcvdz-dk',
     $ 'aug-cc-pcvdz-dk','aug-cc-pcvdz-dk'/
c
c
      character*256 filename
      double precision cart(3),charge
      character*16 tag,element
      character*16 atom_lab(nw_max_atom),symbol
      character*16 calctype
      character*256 mydir
      integer errcode, basis, i, j, nbas
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
      errcode=2000
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
c     call ccca_err(rtdb_cget(rtdb,'ccca:type',1,calctype),errcode)
c     if (calctype.eq.'energy') then
c       osegment=.false.
c     else
        osegment=.true.
c     endif
      orel=.false.
      mydir=' '
      call bas_set_library_name(mydir)
      filename = user_library_name
c
c     Set up basis sets and store to RTDB
c
c     Loop over the six basis set types
c      Create basis
c      Add element exp and coeff from library
c      Store basis
c      Destroy basis
c
      nbas=0
      if (tm) nbas=6 
      do i=0,5 ! loop over the six basis set classes used
c
      j=i+nbas
      bas_name=bas_data(j*4+1)
      if (.not. bas_create(basis, bas_name)) call
     $  errquit('ccca_setup_basis failed to create basis', 0,BASIS_ERR) 
      if (.not.bas_add_ucnt_init(basis)) call errquit
     &    ('ccca_setup_basis: failed to init add_ucnt ',911,
     &       BASIS_ERR)
c
      do 30 idx=1, nunique
        if (int_mb(k_q+idx-1) .lt. 2) then
          mybas=bas_data(j*4+2)
        elseif (int_mb(k_q+idx-1) .lt. 11) then
          mybas=bas_data(j*4+3)
        elseif (int_mb(k_q+idx-1) .lt. 19) then
          mybas=bas_data(j*4+4)
        else
          mybas=bas_data(j*4+3)
        endif
        call bas_tag_lib(basis, osegment, atom_lab(idx), atom_lab(idx),
     $        mybas, filename, orel, .false.)
c     load bq to keep DFT quiet
        call bas_tag_lib(basis, osegment,
     B       'bq'//atom_lab(idx), atom_lab(idx),
     $        mybas, filename, orel, .false.)
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
c
      subroutine ccca_get_ngcore(rtdb,tm)
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
      logical tm, tmin, nodezero
      integer errcode, geom
      integer natoms, idx, frozen
      double precision cart(3),charge
      character*16 tag
c
c     Var init
      errcode=3000
      frozen=0
      tm=.false.
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
        elseif (charge .lt. 31) then
           call errquit('ccCA for main group only',0,GEOM_ERR)
        elseif (charge .lt. 37) then
          frozen=frozen+14
c       elseif (charge .lt. 55) then
c         frozen=frozen+18
c       elseif (charge .lt. 87) then
c         frozen=frozen+27
c       else
c         frozen=frozen+43
        else
          call errquit(
     +    'ccca_get_ngcore: elements not implemented in ccCA',
     +     0,GEOM_ERR)
        endif
10    continue
c
c
      if (.not.rtdb_put(rtdb,'ccca:tm',MT_LOG,1,tm))
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      if (.not. geom_destroy(geom))
     $    call errquit('task_ccca: geom_destroy?', 0, GEOM_ERR)
c
      call ccca_err(rtdb_put(rtdb,'ccca:ngcore',MT_INT,1,frozen),
     $              errcode)
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Get the inner noble gas core, return the number of frozen orbs
c
      subroutine ccca_get_ingcore(rtdb)
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
c       elseif (charge .lt. 55) then
c         frozen=frozen+9
c       elseif (charge .lt. 87) then
c         frozen=frozen+18
c       else
c         frozen=frozen+27
        else
          call errquit(
     +    'ccca_get_ngcore: elements not implemented in ccCA',
     +     0,GEOM_ERR)
        endif
10    continue
c
      if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_destroy?', 0, GEOM_ERR)
c
      call ccca_err(rtdb_put(rtdb,'ccca:ingcore',MT_INT,1,frozen),
     $              errcode)
c
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
      logical task_energy_doit
      external task_energy_doit
c
c     Vars
c
      character*256 filename,t2file,newname
      character*32 nwctheory
      character*8 reftype,tceio,calctype,model,perturbation
      character*18 tag,tag2,vectors_in, frozen_key
      character*80 project(2), vectors_out
      logical result,nodezero,fexist,status
      integer errcode,nopen,natoms,maxiter/100/
      double precision e_comp 
      double precision grad(1536)
c
c     Var init
c
      errcode=5000
      nopen=0
      nwctheory='        '
      reftype='      '
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
c     closed shell: mp2/mbpt2 and tce-ccsd(t)/ccsd(t) 
c                   mp2 and tce-ccsd(t) is default
c     open shell:   rohf-based mbpt2 and tce-ccsd(t) are default
c                   uhf-based mp2/mbpt2 and tce-ccsd(t) is optional              
c
c     If the user requested a specific ref wfn, set it now
c
      if (rtdb_get(rtdb, 'ccca:nopen',MT_INT,1,nopen)) then
         if (nopen.gt.0) then
           reftype='rohf'
           call ccca_err(rtdb_put(rtdb,'scf:nopen',MT_INT,1,
     &         nopen), errcode)
           call ccca_err(rtdb_cput(rtdb,'scf:scftype',1,
     &         reftype), errcode)
         endif
      endif
      if (rtdb_cget(rtdb,'ccca:refwfn',1,reftype)) then
         if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, reftype)) 
     &     call errquit('rtdb_cput in run_calc failed',0,RTDB_ERR)
      endif
      if (theory.eq.'mp') then
        if (.not.rtdb_cget(rtdb,'ccca:mp2type',1,nwctheory)) then
          nwctheory='mp2'
        endif
        if (reftype.eq.'rohf') then
          nwctheory='tce'
        endif
      else if (theory.eq.'cc') then
        if ((.not.rtdb_cget(rtdb,'ccca:cctype',1,nwctheory)).or.
     &     (nopen.gt.0)) then
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
c
      if ((nwctheory.eq.'tce').or.(nwctheory.eq.'ccsd(t)')) then
        call ccca_err(rtdb_cput(rtdb,'geometry',1,'geometry_abelian'
     &                    ),errcode)
      endif 
c
c     set the input/output vectors
c
      vectors_in='atomic'
      vectors_out=destfile(1:inp_strlen(destfile))
      if (oldbasis.ne.'none') then
        vectors_in='project'
        project(1)=oldbasis(1:inp_strlen(oldbasis))
        project(2)=srcfile(1:inp_strlen(srcfile))
        if (.not. rtdb_cput(rtdb, 'vectors:project', 2, project))
     $     call errquit('vectors_input: failed to store project',0,
     &       RTDB_ERR)
      endif
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
        if (.not.rtdb_put(rtdb,'tce:2eorb',mt_log,1,.true.))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:2emet',mt_int,1,14))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:idiskx',mt_int,1,1))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:recompf',mt_log,1,.true.))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:fragment',mt_int,1,-1))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:tilesize',mt_int,1,10))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        if (.not.rtdb_put(rtdb,'tce:tilesize',mt_int,1,5))
     1    call errquit('tce_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tce:densmat',mt_log,1,.false.))
     1   call errquit('tce_input: failed writing to rtdb',0,
     2   rtdb_err)
      else if (nwctheory.eq.'mp2') then
         write(frozen_key,'(a,'':'',a)') 'mp2','frozen core'
      else if (nwctheory.eq.'ccsd(t)') then
         write(frozen_key,'(a,'':'',a)') 'ccsd','frozen core'
         if (.not. rtdb_put(rtdb, 'ccsd:maxiter', mt_int, 1, maxiter))
     $      call errquit('run_calc: failed writing rtbd',0,RTDB_ERR)
      endif
      if (.not. rtdb_cput(rtdb,'task:theory',1,nwctheory)) call
     $    errquit('run_calc: failed writing task_theory',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb, frozen_key, mt_int, 1, fcore))
     $    call errquit('run_calc: failed writing rtbd',0,RTDB_ERR)
c
c     here's where we decide to do a gradient or an energy only
c
      call ccca_err(rtdb_cget(rtdb,'ccca:intype',1,calctype),errcode)
      if (calctype.eq.'energy') then
        call ccca_err(task_energy_doit(rtdb,nwctheory,e_comp),errcode)
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
      if (calctype .eq. 'gradient') then
c
c       Load geometry  
c
        if (.not. geom_create(geom, 'geometry'))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('task_ccca: no geometry ', 0, RTDB_ERR)
c
c       Get number of unique centers
c
        if (.not. geom_ncent(geom, natoms)) call errquit
     $   ('ccca_setup_basis: geom_ncent failed', 0, GEOM_ERR)
c
c       Destroy geometry
c
        if (.not. geom_destroy(geom))
     $     call errquit('task_ccca: geom_create?', 0, GEOM_ERR)
c
        newname=' '
        newname=name(1:len(name))//':gradient'
        call ccca_err(rtdb_put(rtdb,newname(1:len(newname)),
     +                               MT_DBL,natoms*3,grad),errcode)
      endif
      call ccca_err(rtdb_put(rtdb,name(1:len(name)),MT_DBL,1,e_comp)
     +              ,errcode)
c
c     Reset geometry to point to default one
c
      if (((nwctheory.eq.'tce').or.(nwctheory.eq.'ccsd(t)'))) then
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
      implicit none
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
      integer errcode/7000/
      double precision HF_D, HF_T, HF_Q, E_D, E_T, E_Q, ErefHF, ErefS3, 
     +                 ErefS4, ErefP, ErefPS3, zpe, hthermal, CV, CVFC, 
     +                 CVcorr, DK, DKcorr, mp2vtz, CC, CCcorr, sf
      character*4 variant/'P'/
      character*80 calctype
      logical nodezero, parsable, therm, tm 
      double precision AUKCAL 
      parameter (AUKCAL=627.5093314d0)
c
c     Var init
      errcode=7000
      variant='P'
      parsable=.false.
      therm=.false.
      tm=.false.
      nodezero=(ga_nodeid().eq.0)
c
c     Get the values into vars
      call ccca_err(rtdb_get(rtdb,'ccca:hf/avdz',MT_DBL,1,HF_D),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:hf/avtz',MT_DBL,1,HF_T),
     +              errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:hf/avqz',MT_DBL,1,HF_Q),
     +              errcode)
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
      E_D = E_D - HF_D
      E_T = E_T - HF_T
      E_Q = E_Q - HF_Q
c
c     fit function: E(n) = E(inf) + B*exp(-1.63*n)
      ErefHF=-0.243672155849*HF_T+1.243672155849*HF_Q
c
c     fit function: E(lmax) = E(inf) + B/(lmax)^3
      ErefS3=ErefHF-0.729729729730*E_T+1.729729729730*E_Q
c
c     fit function: E(lmax) = E(inf) + B/(lmax+1/2)^4
      ErefS4=ErefHF-0.577163461538*E_T+1.577163461540*E_Q
c
c     fit function: E(n) = E(inf) + B*exp(-(n-1)) + C*exp(-(n-1)**2)
      ErefP=ErefHF+0.034867106002*E_D-0.711622433526*E_T
     $      +1.676755327520*E_Q
      ErefPS3=(ErefP+ErefS3)*0.5
      DKcorr=DK-mp2vtz
      CVcorr=CV-E_T-HF_T
      CCcorr=CC-mp2vtz
c
c     ...and report them
      if (nodezero) then
        write(luout,*) ''
        write(luout,*) 'ccCA-P  reference energy = ',ErefP
        write(luout,*) 'ccCA-S3 reference energy = ',ErefS3
        write(luout,*) 'ccCA-S4 reference energy = ',ErefS4
        write(luout,*) 'ccCA-PS3 reference energy= ',ErefPS3
        write(luout,*) 'DK correction            = ',DKcorr
        write(luout,*) 'CCSD(T) correction       = ',CCcorr
        write(luout,*) 'CV correction            = ',CVcorr
        write(luout,*) '---------------------------'
        write(luout,*) 'Total ccCA-P   energy    = ',
     +              ErefP+DKcorr+CCcorr+CVcorr
        write(luout,*) 'Total ccCA-S3  energy    = ',
     +              ErefS3+DKcorr+CCcorr+CVcorr
        write(luout,*) 'Total ccCA-S4  energy    = ',
     +              ErefS4+DKcorr+CCcorr+CVcorr
        write(luout,*) 'Total ccCA-PS3 energy    = ',
     +              ErefPS3+DKcorr+CCcorr+CVcorr
        write(luout,*) ''
      endif
c
c     Check for availability of thermochem
      if (rtdb_get(rtdb,'vib:zpe',MT_DBL,1,zpe)
c    +   .and.rtdb_get(rtdb,'vib:hthermal',MT_DBL,1,hthermal)
     +  ) then
        call ccca_err(rtdb_get(rtdb,'vib:scalefreq',MT_DBL,1,sf),
     &           errcode)
        zpe = zpe * sf
c     following line commented . See
c     https://nwchemgit.github.io/Special_AWCforum/st/id2417/index.html
c        zpe = zpe / AUKCAL
c       hthermal = hthermal / AUKCAL
        therm=.true.
        if (nodezero) then
        write(luout,*)'Thermochemistry available:'
        write(luout,*)'           ZPE   =',zpe    
        write(luout,*)'ccCA-P   E+ZPE   =',ErefP+DKcorr+CCcorr+CVcorr+
     &                                     zpe
        write(luout,*)'ccCA-S3  E+ZPE   =',ErefS3+DKcorr+CCcorr+CVcorr+
     &                                     zpe
        write(luout,*)'ccCA-S4  E+ZPE   =',ErefS4+DKcorr+CCcorr+CVcorr+
     &                                     zpe
        write(luout,*)'ccCA-PS3 E+ZPE   =',ErefPS3+DKcorr+CCcorr+CVcorr+
     &                                     zpe
c       write(luout,*)'        Enthalpy =',hthermal    
c       write(luout,*)'ccCA-P  Enthalpy =',
c    +             ErefP+DKcorr+CCcorr+CVcorr+hthermal
c       write(luout,*)'ccCA-S3 Enthalpy =',
c    +             ErefS3+DKcorr+CCcorr+CVcorr+hthermal
c       write(luout,*)'ccCA-S4 Enthalpy =',
c    +             ErefS4+DKcorr+CCcorr+CVcorr+hthermal
c       write(luout,*)'ccCA-PS3 Enthalpy=',
c    +             ErefPS3+DKcorr+CCcorr+CVcorr+hthermal
c       write(luout,*)''
        endif
      endif
c
c     Now set the proper rtdb entries
      if (rtdb_cget(rtdb,'ccca:variant',1,variant)) then
        if (nodezero) 
     +  write(luout,*) 'User requested ccCA-'//variant//' energy'
      endif
      if ((variant .eq. 'p') .or. (variant .eq. 'P')) then
        call ccca_err(rtdb_put(rtdb,'task:energy',MT_DBL,1,
     +                ErefP+DKcorr+CCcorr+CVcorr),errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:energy',MT_DBL,1,
     +                ErefP+DKcorr+CCcorr+CVcorr),errcode)
      else if ((variant .eq. 's3') .or. (variant .eq. 'S3')) then
        call ccca_err(rtdb_put(rtdb,'task:energy',MT_DBL,1,
     +                ErefS3+DKcorr+CCcorr+CVcorr),errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:energy',MT_DBL,1,
     +                ErefS3+DKcorr+CCcorr+CVcorr),errcode)
      else if ((variant .eq. 's4') .or. (variant .eq. 'S4')) then
        call ccca_err(rtdb_put(rtdb,'task:energy',MT_DBL,1,
     +                ErefS4+DKcorr+CCcorr+CVcorr),errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:energy',MT_DBL,1,
     +                ErefS4+DKcorr+CCcorr+CVcorr),errcode)
      else if ((variant .eq. 'ps3') .or. (variant .eq. 'PS3')) then
        call ccca_err(rtdb_put(rtdb,'task:energy',MT_DBL,1,
     +                ErefPS3+DKcorr+CCcorr+CVcorr),errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:energy',MT_DBL,1,
     +                ErefPS3+DKcorr+CCcorr+CVcorr),errcode)
      else
        call errquit('ccca: unknown ccCA variant requested',0,INPUT_ERR)
      endif
      if (nodezero) 
     &   write(luout,*) 'Wrote ccCA-'//variant//' energy to the RTDB'
c
c     check for parsable output
c     load appropriate vars before entering node 0 print block
      if (rtdb_get(rtdb,'ccca:parsable',MT_LOG,1,parsable)
     +   .and.nodezero) then
        write(luout,*) ''
        write(luout,*)'CCCA DATA BLOCK BEGINS'
        write(luout,*)'hf/avdz=',HF_D
        write(luout,*)'hf/avtz=',HF_T
        write(luout,*)'hf/avqz=',HF_Q
        write(luout,*)'mp2(fc)/avdz=',E_D+HF_D
        write(luout,*)'mp2(fc)/avtz=',E_T+HF_T
        write(luout,*)'mp2(fc)/avqz=',E_Q+HF_Q
        write(luout,*)'mp2(fc)/vtzdk=',DK
        write(luout,*)'mp2(fc1)/acvtz=',CV
        write(luout,*)'mp2(fc)/vtz=',mp2vtz
        write(luout,*)'ccsd(t)/vtz=',CC
        write(luout,*)'Eref(S3)=',ErefS3
        write(luout,*)'Eref(S4)=',ErefS4
        write(luout,*)'Eref(P)=',ErefP
        write(luout,*)'Eref(PS3)=',ErefPS3
        write(luout,*)'DK correction=',DKcorr
        write(luout,*)'CCSD(T) correction=',CCcorr
        write(luout,*)'CV correction=',CVcorr
        write(luout,*)'Etot(S3)=',ErefS3+DKcorr+CCcorr+CVcorr
        write(luout,*)'Etot(S4)=',ErefS4+DKcorr+CCcorr+CVcorr
        write(luout,*)'Etot(P)=',ErefP+DKcorr+CCcorr+CVcorr
        write(luout,*)'Etot(PS3)=',ErefPS3+DKcorr+CCcorr+CVcorr
        if (therm) then
          write(luout,*)'Scale factor=',sf
          write(luout,*)'ZPE=',zpe
c         write(luout,*)'H therm corr=',hthermal
        endif
        write(luout,*)'CCCA DATA BLOCK ENDS'
        call util_flush(luout)
      endif
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
      integer nopen,iline,errcode/9000/
      logical dotherm/.true./,res
      double precision sf
      logical nodezero
      logical  lcfac(40),lxfac(40),nlcfac(40),nlxfac(40),xccomb(40)
      double precision cfac(40),xfac(40)
c
c     Externs
      logical task_optimize,task_freq
      external task_optimize,task_freq
c
c     Var init
      errcode=9000
      nopen=0
      nodezero=(ga_nodeid().eq.0)
      if (rtdb_get(rtdb, 'ccca:nopen', MT_INT, 1,nopen)) then 
        call ccca_err(rtdb_put(rtdb,'dft:mult',MT_INT,1,
     &    (nopen+1)), errcode)
      endif
c
c     Set the basis set
      call ccca_err(rtdb_cput(rtdb,'ao basis',1,'cc-pvtz'),errcode)
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
      if(rtdb_get(rtdb,'ccca:therm',MT_LOG,1,dotherm)) then
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
        if (rtdb_get(rtdb,'vib:scalefreq',MT_DBL,1,sf)) then
          if (nodezero) 
     +    write(luout,*)'ccCA Warning: nonstandard scale factor used(',
     +         sf,')'
        else
c       new scale factor for the cc-pVTZ basis set
          sf=0.9889d0
          call ccca_err(rtdb_put(rtdb,'vib:scalefreq',MT_DBL,1,sf),
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
c       write (luerr,*) 'ccca error code=',code
        write (luout,*) 'ccca error code=',code
        call errquit('ccca module generated an error',0,INPUT_ERR)
      else
        code=code+1
      endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Calculations of the ccCA gradients and/or energies (default P)
c
      subroutine ccca_energy(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb
      double precision e_hf,e_mp2
      integer stage, nopen, errcode 
      logical nodezero,tm,readmp2,restart,status
      character*80 calctype, cctheory, geom_abel 
c
c     Externals needed
c
      integer ngcore, ingcore 
      character*256 movecsfile(2,6)
      character*256 avdz_file, avtz_file, avqz_file,
     +              vtz_file, vtzdk_file, acvtz_file,
     +              avdzdk_file, avtzdk_file, avqzdk_file,
     +              acvdzdk_file
      character*256 readbas,readvec_file
c
c     Var init
      errcode=100
      nodezero=(ga_nodeid().eq.0)
c
      stage=0
      nopen=0
      tm=.false.
      readmp2=.false.
      restart=.false.
      cctheory='       '
      readbas=''
      readvec_file=''
c
c
c     restart calculations work only for energy calculation
      if (rtdb_get(rtdb,'ccca:restart',MT_LOG,1,restart)) then
        if (restart)
     $    call ccca_err(rtdb_get(rtdb,'ccca:stage',MT_INT,1,stage),0)
      endif
      call ccca_err(rtdb_cget(rtdb,'ccca:intype',1,calctype),errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:tm',MT_LOG,1,tm),errcode) 
      if (.not.((rtdb_cget(rtdb,'ccca:initbas',1,readbas)).and.
     $    (rtdb_cget(rtdb,'ccca:initmovec',1,readvec_file)))) then
c       readvec_file='none'
        readbas='none'
      endif
c
      call ccca_err(rtdb_cget(rtdb,'ccca:movecsfile',12,movecsfile)
     &     ,errcode)
      avdz_file=movecsfile(1,1)
      avtz_file=movecsfile(1,2)
      avqz_file=movecsfile(1,3)
      vtz_file=movecsfile(1,4)
      vtzdk_file=movecsfile(1,5)
      acvtz_file=movecsfile(1,6)
      call ccca_err(rtdb_get(rtdb,'ccca:ngcore',MT_INT,1,ngcore)
     &     ,errcode)
      call ccca_err(rtdb_get(rtdb,'ccca:ingcore',MT_INT,1,ingcore)
     &     ,errcode)
c
      if (rtdb_get(rtdb,'ccca:nopen',MT_INT,1,nopen).and.nodezero) then
        if (nopen.gt.0)
     &    write(luout,*) 'open shell systems'
      endif
      if (rtdb_cget(rtdb,'ccca:cctype',1,cctheory)) then
        if ((cctheory.eq.'ccsd(t)').and.(calctype.eq.'energy')
     $       .and.(nopen.eq.0)) then
          readmp2=.true. 
          if (nodezero) then
              write(luout,*) 
     $        'MP2/cc-pVTZ(-DK) energy will be read from CCSD(T)' 
          endif
        endif
      endif
c
c     use abelian geometry for ccsd(t) and tce methods
c
      call ccca_force_abelian(rtdb)
c
c     The first calc has no orb proj guess, so it's special
c     stage 1: MP2(fc)/aVDZ
      if (.not.(stage.gt.1)) then
        if (nodezero) then
           write(luout,*) 
     $     'ccCA: MP2 calculation with aug-cc-pVDZ basis'
        endif
        call ccca_run_calc(rtdb,readbas,readvec_file,'aug-cc-pvdz',
     +           avdz_file,'mp',ngcore,'ccca:mp2(fc)/avdz')
        stage=1
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage)
     +                ,errcode)
        call ccca_err(rtdb_get(rtdb,'scf:energy',MT_DBL,1,e_hf)
     +              ,errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:hf/avdz',MT_DBL,1,e_hf)
     +              ,errcode)
      endif
c
c     stage 2: MP2(fc)/aVTZ
      if (stage.lt.2) then
        if (nodezero) then
           write(luout,*) 
     $     'ccCA: MP2 calculation with aug-cc-pVTZ basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvdz',avdz_file,
     +     'aug-cc-pvtz',avtz_file,'mp',ngcore,'ccca:mp2(fc)/avtz')
        stage=2
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
        call ccca_err(rtdb_get(rtdb,'scf:energy',MT_DBL,1,e_hf)
     +              ,errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:hf/avtz',MT_DBL,1,e_hf)
     +              ,errcode)
      endif
c
c     stage 3: MP2(fc)/aVQZ
      if (stage.lt.3) then
        if (nodezero) then
           write(luout,*) 
     $     'ccCA: MP2 calculation with aug-cc-pVQZ basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvtz',avtz_file,
     +     'aug-cc-pvqz',avqz_file,'mp',ngcore,'ccca:mp2(fc)/avqz')
        stage=3
        call ccca_err(rtdb_put(rtdb,'ccca:stage',MT_INT,1,stage),
     +                errcode)
        call ccca_err(rtdb_get(rtdb,'scf:energy',MT_DBL,1,e_hf)
     +              ,errcode)
        call ccca_err(rtdb_put(rtdb,'ccca:hf/avqz',MT_DBL,1,e_hf)
     +              ,errcode)
      endif
c
c     Currently, I cannot find a reliable way of the getting
c     the mp2 component of the CCSD(T)/VTZ calc
c     So it is run separately below
c
c     stage 4: MP2(fc)/VTZ
      if (stage.lt.4) then
        if (.not.readmp2) then
          if (nodezero) then
             write(luout,*) 
     $       'ccCA: MP2 calculation with cc-pVTZ basis'
          endif
          call ccca_run_calc(rtdb,'aug-cc-pvdz',avdz_file,
     +       'cc-pvtz',vtz_file,'mp',ngcore,'ccca:mp2(fc)/vtz')
        endif
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
        if (nodezero) then
           write(luout,*) 
     $     'ccCA: Relativity through Douglas-Kroll MP2 calculation'
           write(luout,*) '      with cc-pVTZ-DK basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvdz',avdz_file,
     +       'cc-pvtz-dk',vtzdk_file,'mp',ngcore,'ccca:mp2(fc)/vtzdk')
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
        if (nodezero) then
           write(luout,*) 
     $     'ccCA: Calculation of core-valence correction at MP2 level'
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
        if (nodezero) then
           write(luout,*)
     $     'ccCA: CCSD(T) calculation with cc-pVTZ basis'
        endif
        call ccca_run_calc(rtdb,'aug-cc-pvdz',avdz_file,
     +           'cc-pvtz',vtz_file,'cc',ngcore,'ccca:ccsd(t)/vtz')
        if (readmp2) then
          call ccca_err(rtdb_get(rtdb,'ccsd:mp2 energy',MT_DBL,1,e_mp2)
     +              ,errcode)
          call ccca_err(rtdb_put(rtdb,'ccca:mp2(fc)/vtz',MT_DBL,1,
     +              e_mp2),errcode)
        endif
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
      if (nodezero) then
        write(luout,*)
     $     'ccCA: calculations done, now printing results'
      endif
      call ccca_print(rtdb)
c
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Calculations of the ccCA gradients and energies 
c
C> \ingroup task
C> @{
      logical function task_ccca_energy(rtdb, theory, energy)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb
      double precision energy 
      character*32 theory 
      integer errcode 
c
c     Externals needed
c
c     Var init
      errcode=300
c
      call ccca_energy(rtdb)
c
c
      call ccca_err(rtdb_get(rtdb, 'ccca:energy', mt_dbl, 1, energy),
     $              errcode)
c
      task_ccca_energy = .true. 
      return
      end
C> @}
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccCA input block 
c
      subroutine ccca_input(rtdb)
*
*
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input] data base handle
c
c
      logical therm, all_print 
      integer nopen
      character*32 opttype, mp2type, cctype, refwfn
      character*32 initbas, initmovec 
      character*32 test
c
      initbas    = ' '
      initmovec  = ' '
c
      call rtdb_clean(rtdb, 'ccca:') ! Delete ALL previous ccca parameters
c
c     CCCA 
c        [(ENERGY||OPTIMIZE) default ENERGY]
c        [(DFT||DIRECT) default DFT]
c        [(MP2||MBPT2) default MP2]
c        [(RHF||ROHF||UHF) default RHF]
c        [(CCSD(T)||TCE) default CCSD(T)]
c        [NOPEN <integer number of unpaired electrons default 0>]
c        [(THERM||NOTHERM) default THERM]
c        [(PRINT||NOPRINT) default NOPRINT]
c        [BASIS <basis name for orbital projection guess>]
c        [MOVEC <file name for orbital projection guess>]
c     END
c
c
c ------------------
c Set default values
c ------------------
c
      call ccca_input_default(rtdb) 
c
c ----------
c Read input
c ----------
c
 10   if (.not. inp_read())
     1  call errquit('ccca_input: failed reading input',0,
     2  RTDB_ERR)
      if (.not. inp_a(test))
     1  call errquit('ccca_input: failed reading keyword',0,
     2  RTDB_ERR)
      if (inp_compare(.false.,test,'energy')) then
        if (.not.rtdb_cput(rtdb,'ccca:type',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'optimize')) then
        if (.not.rtdb_cput(rtdb,'ccca:type',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'dft')) then
        if (.not.rtdb_cput(rtdb,'ccca:opttype',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'direct')) then
        if (.not.rtdb_cput(rtdb,'ccca:opttype',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'mp2')) then
        if (.not.rtdb_cput(rtdb,'ccca:mp2type',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'mbpt2')) then
        mp2type='tce'
        if (.not.rtdb_cput(rtdb,'ccca:mp2type',1,mp2type)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'rhf')) then
        if (.not.rtdb_cput(rtdb,'ccca:refwfn',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'rohf')) then
        if (.not.rtdb_cput(rtdb,'ccca:refwfn',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'uhf')) then
        if (.not.rtdb_cput(rtdb,'ccca:refwfn',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'ccsd(t)')) then
        if (.not.rtdb_cput(rtdb,'ccca:cctype',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'tce')) then
        if (.not.rtdb_cput(rtdb,'ccca:cctype',1,test)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
      else if (inp_compare(.false.,test,'nopen')) then
        if (.not.inp_i(nopen)) 
     1  call errquit('ccca_input: nopen value not found',0,
     2  RTDB_ERR)
        if (.not.rtdb_put(rtdb,'ccca:nopen',mt_int,1,nopen)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2      RTDB_ERR)
      else if (inp_compare(.false.,test,'therm')) then
        therm=.true.
        if (.not.rtdb_put(rtdb,'ccca:therm',mt_log,1,therm)) then
          call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        endif
      else if (inp_compare(.false.,test,'notherm')) then
        therm=.false.
        if (.not.rtdb_put(rtdb,'ccca:therm',mt_log,1,therm)) then
          call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        endif
      else if (inp_compare(.false.,test,'print')) then
        all_print=.true.
        if (.not.rtdb_put(rtdb,'ccca:parsable',mt_log,1,all_print)) then
          call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        endif
      else if (inp_compare(.false.,test,'noprint')) then
        all_print=.false.
        if (.not.rtdb_put(rtdb,'ccca:parsable',mt_log,1,all_print)) then
          call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
        endif
      else if (inp_compare(.false.,test,'basis')) then
        if (.not.inp_a(initbas)) 
     1  call errquit('ccca_input: basis value not found',0,
     2  RTDB_ERR)
        if (.not.rtdb_cput(rtdb,'ccca:initbas',1,initbas)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2      RTDB_ERR)
      else if (inp_compare(.false.,test,'movec')) then
        if (.not.inp_a(initmovec)) 
     1  call errquit('ccca_input: movec value not found',0,
     2  RTDB_ERR)
        if (.not.rtdb_cput(rtdb,'ccca:initmovec',1,initmovec)) 
     1    call errquit('ccca_input: failed writing to rtdb',0,
     2      RTDB_ERR)
c
c     END
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit('ccca_input: unknown directive',0,INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccCA input default 
c
      subroutine ccca_input_default(rtdb)
*
*
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input] data base handle
c
c
      logical therm 
      integer nopen
      character*32 calctype, opttype, mp2type, cctype, refwfn
c
c     CCCA 
c        [(ENERGY||OPTIMIZE) default ENERGY]
c        [(DFT||DIRECT) default DFT]
c        [(MP2||MBPT2) default MP2]
c        [(RHF||ROHF||UHF) default RHF]
c        [(CCSD(T)||TCE) default CCSD(T)]
c        [NOPEN <integer number of unpaired electrons default 0>]
c        [(THERM||NOTHERM) default THERM]
c        [(PRINT||NOPRINT) default NOPRINT]
c     END
c
c
c ------------------
c Set default values
c ------------------
c
c     ENERGY or OPTIMIZE
c
      calctype='energy'
      if (.not.rtdb_cput(rtdb,'ccca:type',1,calctype))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     DFT(B3LYP) or direct geometry optimization
c
      opttype='dft'
      if (.not.rtdb_cput(rtdb,'ccca:opttype',1,opttype))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     MP2 or TCE-MBPT2 
c
      mp2type='mp2'
      if (.not.rtdb_cput(rtdb,'ccca:mp2type',1,mp2type))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     RHF, ROHF or UHF (reference wavefunction) 
c
      refwfn='rhf'
      if (.not.rtdb_cput(rtdb,'ccca:refwfn',1,refwfn))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     CCSD(T) or TCE-CCSD(T) 
c
      cctype='ccsd(t)'
      if (.not.rtdb_cput(rtdb,'ccca:cctype',1,cctype))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c     number of unpaired electrons
c
      nopen=0
      if (.not.rtdb_put(rtdb,'ccca:nopen',mt_int,1,nopen))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2  RTDB_ERR)
c
c     Thermal corrections 
c
      therm=.true.
      if (.not.rtdb_put(rtdb,'ccca:therm',mt_log,1,therm))
     1  call errquit('ccca_input: failed writing to rtdb',0,
     2    RTDB_ERR)
c
c
c ------
c Return
c ------
c
 20   return
      end
c $Id$
