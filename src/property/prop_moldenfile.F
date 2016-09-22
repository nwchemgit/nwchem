      subroutine  prop_moldenfile(rtdb)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
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
      integer l_eval,k_eval,l_ir,k_ir,k_occ,l_occ
      integer iset,iptr
      external ga_create_atom_blocked,movecs_read_header,movecs_read
      logical ga_create_atom_blocked,movecs_read_header,movecs_read
      character*15 subname
      parameter  (subname='prop_moldenfile')
c
c>  load geometry and symmetry info
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
           write(6,*)' using vec ',vec_file,'@@'
      call util_file_name_resolve(vec_file, .false.)
c> movecs 
      if (.not. movecs_read_header(vec_file, title_vec, basis_name,
     &    scftype_vecs, nbf, nset, nmo, 2))
     &  call errquit(subname//': basis set error:', 86, basis_err)

      g_movecs(1) = ga_create_atom_blocked(geom, basis, 'vecs1')
      call ga_zero(g_movecs(1))
      if (nset.eq.2) then
        g_movecs(2) = ga_create_atom_blocked(geom, basis, 'vecs2')
        call ga_zero(g_movecs(2))
      endif

c     allocate eigenvalues
      if (.not.MA_Push_Get(MT_Dbl,nbf*nset,'eval',l_eval,k_eval))
     &     call errquit(subname//': cannot allocate eval',0, MA_ERR)
      if (.not. MA_Push_Get(MT_Dbl,nbf*nset,'mo occ',l_occ,k_occ))
     &   call errquit('dft_scf: failed to alloc',999, MA_ERR)
c> read  occupation numbers, evals and eigenvecs
      do iset = 1, nset
         iptr=(iset-1)*nbf
        if (.not. movecs_read(vec_file, iset, dbl_mb(k_occ+iptr),
     G        dbl_mb(k_eval+iptr), g_movecs(iset+iptr))) call errquit
     $              (subname//': failed read  MOs set ',iset,
     $               INPUT_ERR)
      enddo



      if (.not.MA_Push_Get(mt_int,nset*nbf,'dft:irreps',l_ir,k_ir))
     &   call errquit('dft_scf: cannot allocate irreps',0, MA_ERR)
      call ifill(nset*nbf,1,int_mb(k_ir),1)
c     get irreps
      call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     $     0, 'junk', .true., int_mb(k_ir))

      write(6,*) ' nbf ',nbf
      write(6,*) ' nset ',nset
      do iset=1,nbf
         write(6,*) ' i = ',iset, ' occ = ',dbl_mb(k_occ+iset-1),
     R        int_mb(k_ir+iset-1)
         enddo


      if (.not. rtdb_cget(rtdb, theory(1:inp_strlen(theory))//
     A        ':scftype', 1, scftype)) then
         call errquit
     &    (subname//':error need scftype',0,0)
      endif
      write(6,*) ' scftype ',scftype,'@'

      call dft_nwmolden(rtdb,geom,
     I     nset,dbl_mb(k_eval),dbl_mb(k_occ),int_mb(k_ir),
     N     nbf, nmo,
     I     g_movecs,scftype)


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
