      subroutine  prop_polfromsos(rtdb)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "sym.fh"
#include "util.fh"
c#include "geomP.fh"
      integer rtdb ! [in]
c
      integer geom, basis
      character*32  theory
      integer g_movecs(2)
      character*255 vec_file
      character*255 title_vec,basis_name
      character*20 scftype_vecs,scftype
      integer nset
      integer nbf
      integer nmo(2)
      integer l_eval,k_eval(2),l_ir,k_ir,k_occ,l_occ
      integer iset,iptr
      external ga_create_atom_blocked,movecs_read_header,movecs_read
      integer ga_create_atom_blocked
      logical movecs_read_header,movecs_read
      logical int_normalize
      external int_normalize
      logical oadapt
      character*15 subname
      integer owstart(2),owend(2),nocc(2),i
      character*255 vfilename
      double precision nocf
      parameter (subname='prop_polfromsos')
c
c  load geometry and symmetry info
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit(subname//': geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit(subname//': no geometry ', 0, RTDB_ERR)

c>  load the basis set and get info about it
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit(subname//': bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) then
         if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $        call errquit(subname//': no mo or ao basis set', 0,
     &       RTDB_ERR)
      endif
      if(.not.int_normalize(rtdb,basis))
     &     call errquit(subname//'int_normalize failed ', 5, BASIS_ERR)
      call int_init(rtdb,1,basis)


c theory
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('waimfile: theory not specified',555,
     &       INPUT_ERR)
c movecs
      if (.not. rtdb_cget(rtdb, 'prop:vectors',  1,
     $     vec_file)) then
         if (.not. rtdb_cget(rtdb, theory(1:inp_strlen(theory))//
     A          ':output vectors', 1, vec_file)) then
           call util_file_name('movecs', .false.,.false.,vec_file )
        endif
      endif

      call util_file_name_resolve(vec_file, .false.)
c     > movecs
      if (.not. movecs_read_header(vec_file, title_vec, basis_name,
     &    scftype_vecs, nbf, nset, nmo, 2))
     &  call errquit(subname//': basis set error:', 86, basis_err)
      
      if (.not. ga_create(mt_dbl, nbf, nmo(1),
     +     'movad1', nbf, 0, g_movecs(1))) call errquit(
     +              subname//':gacreate failed ',0, GA_ERR)
      call ga_zero(g_movecs(1))
      if (nset.eq.2) then
         if (.not. ga_create(mt_dbl, nbf, nmo(2),
     +        'movad2', nbf, 0, g_movecs(2))) call errquit(
     +        subname//':gacreate failed ',0, GA_ERR)
        call ga_zero(g_movecs(2))
      endif

c     allocate eigenvalues
      if (.not.MA_Push_Get(MT_Dbl,nbf*nset,'eval',l_eval,k_eval(1)))
     &     call errquit(subname//': cannot allocate eval',0, MA_ERR)
      if(nset.eq.2) k_eval(2)=k_eval(1)+nbf
      if (.not. MA_Push_Get(MT_Dbl,nbf*nset,'mo occ',l_occ,k_occ))
     &   call errquit(subname//': failed to alloc',999, MA_ERR)
c     > read  occupation numbers, evals and eigenvecs
      do iset=1,nset
         if (.not.movecs_read(vec_file, iset,dbl_mb(k_occ+(iset-1)*nbf),
     G        dbl_mb(k_eval(iset)), g_movecs(iset))) call errquit
     $        (subname//': failed read  MOs set ',iset,INPUT_ERR)
         nocf=0d0
         do i=1,nmo(iset)
            nocf=nocf+dbl_mb(k_occ+(iset-1)*nbf+i-1)
         enddo
         if(nset.eq.1) nocf=nocf/2d0
            nocc(iset)=nocf
      enddo

      if (.not.MA_Push_Get(mt_int,nset*nbf,'dft:irreps',l_ir,k_ir))
     &   call errquit('dft_scf: cannot allocate irreps',0, MA_ERR)
      call ifill(nset*nbf,1,int_mb(k_ir),1)
c     get irreps
      if (.not.rtdb_get(rtdb, 'dft:adapt',MT_LOG, 1, oadapt))then
         oadapt = sym_number_ops(geom) .gt. 0
      endif
      
      if(oadapt) then
         do iset=1,nset
            call scf_movecs_sym_adapt(basis, g_movecs(iset), .false.,
     $           nmo(iset), 'junk', .true., int_mb(k_ir))
         enddo
       endif
      if (.not. rtdb_cget(rtdb, theory(1:inp_strlen(theory))//
     A        ':scftype', 1, scftype)) then
c     wild guess here to fix mp2
         if (nset .eq. 1) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
         if(ga_nodeid().eq.0) write(6,*) 'missing scftype set to ',
     c        scftype
      endif
      
      if (.not. rtdb_put(rtdb,'dft:polfromsos',mt_log,1,.true.))
     &     call errquit(subname//': rtdb_put ',12,RTDB_ERR)
      
      call util_file_name('vspec',.false.,.false.,vfilename)
      call dft_vspec(rtdb, basis,vfilename,nset,nmo,nbf,
     &  nocc, g_movecs,k_eval,owstart,owend)



        call int_terminate()
      if (.not.(bas_destroy(basis)))
     &    call errquit
     &    (subname//':error bas_destroy',911,GEOM_ERR)

      if(.not.(geom_destroy(geom)))
     &    call errquit
     &    (subname//':error geom_destroy',911, BASIS_ERR)
      if (.not.ma_chop_stack(l_eval))
     &   call errquit(subname//': cannot chop stack',10, MA_ERR)
      if (.not. ga_destroy(g_movecs(1)))
     &   call errquit(subname//': cannot ga_destroy',1, GA_ERR)
      if (nset.eq.2) then
      if (.not. ga_destroy(g_movecs(2)))
     &   call errquit(subname//': cannot ga_destroy',2, GA_ERR)
      endif

      return
      end
