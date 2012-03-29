      subroutine movecs_project(rtdb, nmo, geom, basis, g_vecs,iset,
     A     calling_module, k_eval)
C$Id$
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
c
      integer rtdb, geom, basis, g_vecs(2), k_eval
      double precision one, zero, mone
c
c     Project the MO vectors in old_basis_name into the current basis.  
c     It is no longer assumed that the old basis is smaller than the current 
c     basis, and it now works when linear dependencies have been removed.
c
c     If nnew > nold
c     --------------
c     For i = 1, nold
c        C_new(u,i) = sum(v,w) S^-1(u,v) Snew-old(v,w) C_old(w,i)
c  
c     For i=nold+1, nnew
c        C_new(u,i) = random number in -0.5,0.5 
c
c     If nold > nnew
c     --------------
c     For i = 1, nnew
c        C_new(u,i) = sum(v,w) S^-1(u,v) Snew-old(v,w) C_old(w,i)
c     
c     S^-1 is the inverse overlap in the new basis
c     Snew-old is the mixed basis overlap
c
c     !! THE RESULTING MOS WILL NOT BE ORTHONORMAL
c
      integer iset ! [in] 1 = alpha 2 =beta
      integer basis_handles(2)
      integer old_nbf, old_nsets, old_nmo(2), file_nbf
      integer g_old_vecs
      integer g_over
      integer g_mix
      integer l_old_occ, k_old_occ, l_old_evals, k_old_evals
      integer l_ran, k_ran
      integer old_basis
      integer nbf, i, j, nmo
      double precision temp
      character*255 old_title,basis_name_from_file,old_basis_name,
     A     old_filename
      character*80 project(2)
      character*20 old_scftype
      character*(*) calling_module
      logical oprint, oprint_info

      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      integer ga_create_atom_blocked,anggg
      external ga_create_atom_blocked
      integer nbases,dft_nbases
      external dft_nbases
      logical int_normalize
      external int_normalize

      data   one,  zero,  mone ! Because of KSR2 bug?
     $     /1.0d0, 0.0d0, -1.0d0/

      oprint = util_print('mo guess', print_default)
      oprint_info = (ga_nodeid() .eq. 0 .and. oprint)
      
      if (oprint_info) then
         write(6,*)
         call util_print_centered(6,'Orbital projection guess',15,
     $        .true.)
         call util_flush(6)
      endif

*     read in old vectors

      if (.not. rtdb_cget(rtdb, 'vectors:project', 2, project))
     $     call errquit('movecs_project: no file/basis', 0, RTDB_ERR)
      old_basis_name = project(1)
      old_filename   = project(2)
      call util_file_name_resolve(old_filename, .false.)
*
      if (.not. bas_create(old_basis,old_basis_name)) call errquit
     $     ('movecs_project: bas_create', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, old_basis, old_basis_name))
     $     call errquit('movecs_project: no old basis?', 0, RTDB_ERR)
      if (.not.int_normalize(rtdb,old_basis))
     $     call errquit('movecs_project: int_normalize failed', 0,
     &       INT_ERR)
      if (.not. movecs_read_header(old_filename, old_title,
     $     basis_name_from_file, old_scftype, 
     $     file_nbf, old_nsets, old_nmo, 2)) 
     $     call errquit
     $     ('movecs_project: failed to read old header', 0, INPUT_ERR)
c
      if (.not. (bas_numbf(basis, nbf))) 
     $     call errquit('movecs project: bas_numbf?', 0, BASIS_ERR)
      if (.not. (bas_numbf(old_basis, old_nbf))) 
     $     call errquit('movecs project: bas_numbf?', 0, BASIS_ERR)
c
      if (nbf .lt. old_nbf) then
          write(6,*) 'The current basis is smaller than the basis',
     $               ' you are projecting from !!!'
c          call errquit
c     $     ('movecs_project: new bas smaller than old bas',0) 
      endif
c
      if (file_nbf .ne. old_nbf) call errquit
     $     ('movecs_project: nbf on file inconsistent with basis', 0,
     &       BASIS_ERR)
c
*ga:1:0
c      if (.not. ga_create(mt_dbl, old_nbf, old_nbf, 'Old vectors',
      if (.not. ga_create(mt_dbl, old_nbf, old_nmo(1), 'Old vectors',
     $     0, 0, g_old_vecs)) call errquit
     $     ('movecs_project: failed creating old vectors', 0, GA_ERR)
      if (.not. ma_push_get(mt_dbl, old_nbf, 'Old occ',
     $     l_old_occ, k_old_occ)) call errquit
     $     ('movecs_project: ma old occ', old_nbf, MA_ERR)
*
      if (.not. movecs_read(old_filename, iset, dbl_mb(k_old_occ),
     $     dbl_mb(k_eval), g_old_vecs)) call errquit
     $     ('movecs_project: failed reading old vectors', 0, INPUT_ERR)
*
      if (oprint_info) then
         write(6,1)
     $        old_filename(1:inp_strlen(old_filename)),
     $        old_title(1:inp_strlen(old_title)),
     $        old_nbf,
     $        old_nsets,
     $        old_nmo(1),
     $        old_basis_name(1:inp_strlen(old_basis_name)),
     $        basis_name_from_file(1:inp_strlen(basis_name_from_file))
 1       format(/
     $        ' Filename           = ',a/
     $        ' Old title          = "',a,'"'/
     $        ' Old basis size     = ', i5/
     $        ' Old no. of sets    = ', i5/
     $        ' Old no. of MOs     = ', i5/
     $        ' Old basis name     = ', a/
     $        ' Basis name on file = ', a/)
         call util_flush(6)
*         write(6,*) ' Old eigen values'
*         call output(dbl_mb(k_old_evals),1,old_nbf,1,1,old_nbf,1,1)
*         write(6,*) ' Old occupations'
*         call output(dbl_mb(k_old_occ),1,old_nbf,1,1,old_nbf,1,1)
      end if
*
      if (.not. ma_pop_stack(l_old_occ)) call errquit
     $     ('movecs_project: ma pop', 0, MA_ERR)

*     Multiply by mixed overlap

      basis_handles(1) = basis
      basis_handles(2) = old_basis

*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, old_nbf, 'Mixed overlap',
     $     0, 0, g_over)) call errquit
     $     ('movecs_project: create mixed overlap', 0, GA_ERR)
*ga:1:0
c      if (.not. ga_create(mt_dbl, nbf, old_nbf, 'Temp',
      if (.not. ga_create(mt_dbl, nbf, old_nmo(1), 'Temp',
     $     0, 0, g_mix)) call errquit
     $     ('movecs_project: create temp', 0, GA_ERR)
      call ga_zero(g_mix)
      call ga_zero(g_over)
      call int_terminate
      call int_init(rtdb,2,basis_handles)
      call int_1e_ga(basis, old_basis, g_over, 'overlap', .false.)
      call int_terminate
      if (.not. bas_high_angular(basis,anggg))
     &   call errquit('rdinput: basis set error:', 0, 0)
      if(calling_module.eq.'dft') then
         nbases=dft_nbases(basis_handles)
         call int_init(rtdb, nbases, basis_handles) !  RIGHT?
      else
         call int_init(rtdb, 1, basis) !  WRONG WRONG WRONG
      endif
      if (.not. bas_destroy(old_basis)) call errquit
     $     ('movecs_project: bas_destroy?',0, BASIS_ERR)

*
*      if (ga_nodeid().eq.0) then
*         call ga_summarize(0)
*         write(6,*) ' nbf, old_nbf, g_over, g_old_vecs, g_mix ',
*     $         nbf, old_nbf, g_over, g_old_vecs, g_mix
*         call util_flush(6)
*      endif
*      call ga_sync
*

c      call ga_dgemm('n','n',nbf,old_nbf,old_nbf,one,g_over,
      call ga_dgemm('n','n',nbf,old_nmo(1),old_nbf,one,g_over,
     $     g_old_vecs,zero,g_mix)
      if (.not. ga_destroy(g_over)) call errquit
     $     ('movecs project: ga_destroy?', 0, GA_ERR)
      if (.not. ga_destroy(g_old_vecs)) call errquit
     $     ('movecs project: ga_destroy?', 0, GA_ERR)

*     Make and invert plain overlap

      g_over = ga_create_atom_blocked(geom, basis, 'Overlap')
*ga:1:0
*      if (.not. ga_create(mt_dbl, nbf, nbf, 'Overlap',
*     $     0, 0, g_over)) call errquit
*     $     ('movecs_project: failed making overlap', 0)
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)

      call ga_matpow(g_over, mone, zero)

*     Multiply by inverse overlap

      call ga_zero(g_vecs(iset))

*
*      if (ga_nodeid().eq.0) then
*         call ga_summarize(0)
*         write(6,*) ' nbf, old_nbf, g_over, g_mix, g_vecs ',
*     $         nbf, old_nbf, g_over, g_mix, g_vecs(iset)
*         call util_flush(6)
*      endif
*      call ga_sync
*
      if (old_nmo(1).gt.nmo) old_nmo(1) = nmo
      call ga_dgemm('n','n',nbf,old_nmo(1),nbf,one,g_over,g_mix,
     $     zero,g_vecs(iset))
      if (.not. ga_destroy(g_over)) call errquit
     $     ('movecs project: ga_destroy?', 0, GA_ERR)
      if (.not. ga_destroy(g_mix)) call errquit
     $     ('movecs project: ga_destroy?', 0, GA_ERR)

*     Fill in missing columns with random numbers
      if (nbf.gt.old_nbf) then
      if (.not. ma_push_get(mt_dbl, nbf, 'ran #',
     $     l_ran, k_ran)) call errquit
     $     ('movecs_project: ma ran #', nbf, MA_ERR)
      temp = util_random(ga_nodeid()*314159+1)
c      do i = old_nbf+1+ga_nodeid(), nbf, ga_nnodes()
cbert      do i = old_nbf+1+ga_nodeid(), nmo, ga_nnodes()
       do i = old_nmo(1)+1+ga_nodeid(), nmo, ga_nnodes()
*         write(6,*) ' Randomizing ', i
*         call util_flush(6)
         do j = 1, nbf
            dbl_mb(k_ran+j-1) = (util_random(0) - 0.5d0)
         end do
         call ga_put(g_vecs(iset), 1, nbf, i, i, dbl_mb(k_ran), 1)
      end do
      if (.not. ma_pop_stack(l_ran)) call errquit
     $     ('movecs_project: ma pop?', 0, MA_ERR)
      endif

      end
