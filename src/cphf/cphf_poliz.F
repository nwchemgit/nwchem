c
c     Generalized to handle DFT as well as HF
c
c     BGJ - 8/98
c
      logical function cphf_poliz(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
c !!! BGJ
#include "bgj.fh"
c !!! BGJ
      integer rtdb              ! [input] database handle
c
      logical movecs_read,movecs_read_header, cphf, file_write_ga,
     $     file_read_ga
      external movecs_read,movecs_read_header, cphf, file_write_ga,
     $     file_read_ga
      logical int_normalize
      external int_normalize
      character*(nw_max_path_len) movecs, fname
      character*10 scftype
      character*20 scftype_vecs
      character*255 title_vecs, basis_vecs, title
      integer nbf, nbf_vecs, nsets, nmos(2), nmo
      integer noa, nob, nva, nvb, nvirt, nclosed, nopen
      integer vlen
      integer k_occ_a, l_occ_a, k_eval_a, l_eval_a
      integer k_occ_b, l_occ_b, k_eval_b, l_eval_b
      integer geom, basis
      integer g_vecs_a, g_vecs_b, g_dip(3),ix, g_tmp, g_rhs, g_scratch
      integer g_sol
      integer noc(2)
      logical hf_job
c
c     Compute the SCF polarizability by solving the CPHF equations
c
c     Get MO vectors etc.
c
      hf_job = .not. bgj_have_xc()
      if(.not.rtdb_cget(rtdb, 'title', 1, title)) title = ' '
      if(.not.geom_create(geom, 'geometry'))
     $     call errquit('cphf_poliz: geom_create?', 0, GEOM_ERR)
      if(.not.geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('cphf_poliz: no geometry ', 0, RTDB_ERR)
      if(.not.bas_create(basis, 'ao basis'))
     $     call errquit('cphf_poliz: bas_create?', 0, BASIS_ERR)
      if(.not.bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('cphf_poliz: no ao basis ', 0, RTDB_ERR)
      if (.not.int_normalize(rtdb,basis))
     $     call errquit('cphf_poliz: int_normalize failed ', 0,
     &       UNKNOWN_ERR)
      if (hf_job) then
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs))
     $        call errquit('cphf_poliz: no SCF MO vectors',0, RTDB_ERR)
c         write(*,*)'*** cphf_poliz: hf movecs: ',movecs
         if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype))
     $        call errquit('cphf_poliz: scftype', 0, RTDB_ERR)
      else
c         write(*,*)'*** cphf_poliz: dft case:'
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs))
     $        call errquit('cphf_poliz: no DFT MO vectors',0, RTDB_ERR)
c         write(*,*)'*** cphf_poliz: dft movecs: ',movecs
c !!! Hack to set up scftype properly for later code
         if (bgj_restricted()) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
      endif
c      write(*,*)'cphf_poliz: scftype: ',scftype
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('cphf_poliz: could not get nbf',0, BASIS_ERR)
c
      if(.not.movecs_read_header(movecs,title_vecs,basis_vecs,
     $     scftype_vecs, nbf_vecs,nsets,nmos,2))
     $     call errquit('cphf_poliz: problem with movecs_read_header',1,
     &       INPUT_ERR)
      if(nbf.ne.nbf_vecs)call errquit('cphf: MO vectors mismatch?',0,
     &       INPUT_ERR)
      nmo = nmos(1)
c
      if(.not.ma_push_get(mt_dbl,nbf,'occ alpha',
     $     l_occ_a,k_occ_a))      
     $     call errquit('cphf_poliz ma problem with alpha occ',1,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbf,'alpha evals',
     $     l_eval_a,k_eval_a))
     $     call errquit('cphf_poliz: problem with ma for a evals',1,
     &       MA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nbf,nmo,'alpha_vecs',nbf,0,
     $     g_vecs_a))call errquit('cphf_poliz: ga_create failed',0,
     &       GA_ERR)
      if (scftype .eq. 'UHF') then
*ga:1:0
         if(.not.ga_create(mt_dbl,nbf,nmo,'beta_vecs',nbf,0,
     $        g_vecs_b))call errquit('cphf_poliz: ga_create failed',0,
     &       GA_ERR)
         if(.not.ma_push_get(mt_dbl,nbf,'occ beta',
     $        l_occ_b,k_occ_b))
     $        call errquit('cphf_poliz ma problem with beta occ',1,
     &       MA_ERR)
         if(.not.ma_push_get(mt_dbl,nbf,'beta evals',
     $        l_eval_b,k_eval_b))
     $        call errquit('cphf_poliz: problem with ma for b evals',1,
     &       MA_ERR)
      endif
      if (hf_job) then
         if(scftype.eq.'UHF')then
            if(.not.rtdb_get(rtdb,'scf:nalpha',mt_int,1,noa))
     $           call errquit('no alphas',1, RTDB_ERR)
            if(.not.rtdb_get(rtdb,'scf:nbeta',mt_int,1,nob))
     $           call errquit('no betas',1, RTDB_ERR)
         else if(scftype.eq.'RHF')then
            if(.not.rtdb_get(rtdb,'scf:nclosed',mt_int,1,nclosed))
     $           call errquit('no occ',1, RTDB_ERR)
            noa=nclosed
            nob=nclosed
            nopen = 0
         else
            call errquit('mp2_grad: unknown SCFTYPE ', 0, INPUT_ERR)
         end if
      else
         if (.not. rtdb_get(rtdb,'dft:noc',mt_int,2,noc))
     &        call errquit('cphf_poliz: rtdb_get of noc failed',0,
     &       RTDB_ERR)
c         write(*,*)'cphf_poliz: dft noc:',noc
         noa=noc(1)
         if (scftype.eq.'UHF') then
            nob=noc(2)
         else
            nclosed=noa
            nob=noa
         endif
      endif
      if (bgj_print() .gt. 0)
     &     write(*,*)'*** cphf_poliz: noa nob',noa,nob
      nopen = 0
      nva=nmo-noa
      nvb=nmo-nob
      if (scftype.eq.'UHF') then 
         vlen = noa*nva + nob*nvb
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         nvirt = nmo - nclosed - nopen
         vlen = nclosed*(nopen+nvirt) + nopen*nvirt
         if (bgj_print() .gt. 0)
     &   write(*,*)'nvirt,nmo,nclosed,nopen',nvirt,nmo,nclosed,nopen
      else
         call errquit('cphf: unknown SCF type',0,
     &       INPUT_ERR)
      endif
c
      if(.not.movecs_read(movecs,1,dbl_mb(k_occ_a),
     $     dbl_mb(k_eval_a),g_vecs_a))
     $     call errquit('cphf_poliz: could not read alpha vecs',1,
     &       INPUT_ERR)
      if(scftype.eq.'UHF')then
         if(.not.movecs_read(movecs,2,dbl_mb(k_occ_b),
     $        dbl_mb(k_eval_b),g_vecs_b))
     $        call errquit('cphf_poliz: could not read beta vecs',1,
     &       INPUT_ERR)
      endif
c     
c     Generate dipole integrals
c
      do ix = 1, 3
*ga:1:0
         if(.not. ga_create(mt_dbl,nbf,nbf,' ao integs x ',
     $        0,0,g_dip(ix))) call errquit('cphf_poliz: ga',ix, GA_ERR)
      enddo
c
      call int_init(rtdb,1,basis)
      call int_dip_ga(basis, basis, g_dip(1), g_dip(2), g_dip(3))
      call int_terminate()
c
      if (.not. bas_destroy(basis)) call errquit('cphf_poliz: bas?',0,
     &       BASIS_ERR)
      if (.not. geom_destroy(geom)) call errquit('cphf_poliz: geom?',0,
     &       GEOM_ERR)
c
c     Assemble the RHS vectors ... currently just solve one at a  time
c
*ga:1:0
      if(.not.ga_create(mt_dbl, nbf, nbf,'scratch', 0, 0,g_scratch)) 
     $     call errquit('cphf_poliz: ga create of scratch', nbf, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl, nmo, nmo,'tmp', 0, 0,g_tmp)) 
     $     call errquit('cphf_poliz: ga create of tmp', nbf, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl, vlen, 1,'rhs', 0, 0,g_rhs)) 
     $     call errquit('cphf_poliz: ga create of rhs', vlen, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl, vlen, 1,'sol', 0, 0,g_sol)) 
     $     call errquit('cphf_poliz: ga create of sol', vlen, GA_ERR)
c
      call ga_print(g_vecs_a)
      if (scftype .eq. 'UHF') call ga_print(g_vecs_b)
      do ix = 1, 3
         call two_index_transf(g_dip(ix), g_vecs_a, g_vecs_a, 
     $        g_scratch, g_tmp)
         call ga_print(g_tmp)
         call ga_copy_patch('n', g_tmp, noa+1, nmo, 1, noa,
     $        g_rhs, 1, noa*nva, 1, 1)
         if (scftype .eq. 'UHF') then
            call two_index_transf(g_dip(ix), g_vecs_b, g_vecs_b, 
     $           g_scratch, g_tmp)
            call ga_copy_patch('n', g_tmp, nob+1, nmo, 1, nob,
     $           g_rhs, noa*nva+1, vlen, 1, 1)
         endif
         call cphf_fname('cphf_rhs',fname)
         call ga_print(g_rhs)
         if(.not.file_write_ga(fname,g_rhs))
     $        call errquit('cphf_poliz: could not write cphf_rhs',1,
     &       UNKNOWN_ERR)
         if (.not. cphf(rtdb)) call errquit('cphf_poliz: cphf?',0,
     &       RTDB_ERR)
         call cphf_fname('cphf_sol',fname)
         if(.not.file_read_ga(fname,g_sol))
     $        call errquit('cphf_poliz: could not read cphf_rhs',1,
     &       INPUT_ERR)
         call ga_print(g_sol)
c
         if (scftype .eq. 'UHF') then
            write(6,*) ' POLIZ ', ix, ga_ddot(g_rhs, g_sol)*4.0d0
         else
            write(6,*) ' POLIZ ', ix, ga_ddot(g_rhs, g_sol)*16.0d0
         endif
c
      enddo
c
      if (bgj_print() .gt. 0)
     &     write(*,*)'*** cphf_poliz: done with x y z loop'
      if (.not. ga_destroy(g_vecs_a)) call errquit('cphf_poliz: ga',1,
     &       GA_ERR)
      if (scftype .eq. 'UHF') then
         if (.not. ga_destroy(g_vecs_b))call errquit('cphf_poliz: ga',2,
     &       GA_ERR)
      endif
      if (.not. ga_destroy(g_rhs)) call errquit('cphf_poliz: ga',3,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp)) call errquit('cphf_poliz: ga',4,
     &       GA_ERR)
      if (.not. ga_destroy(g_scratch)) call errquit('cphf_poliz: ga',5,
     &       GA_ERR)
      do ix = 1, 3
         if (.not. ga_destroy(g_dip(ix))) 
     $        call errquit('cphf_poliz: ga',5+ix, GA_ERR)
      enddo
      if (.not. ga_destroy(g_sol)) call errquit('cphf_poliz: ga',6,
     &       GA_ERR)
c
      if (scftype .eq. 'UHF') then
         if (.not.ma_pop_stack(l_eval_b))call errquit('cphf_poliz:ma',1,
     &       MA_ERR)
         if (.not.ma_pop_stack(l_occ_b))call errquit('cphf_poliz:ma',2,
     &       MA_ERR)
      endif
      if (.not.ma_pop_stack(l_eval_a))call errquit('cphf_poliz:ma',3,
     &       MA_ERR)
      if (.not.ma_pop_stack(l_occ_a))call errquit('cphf_poliz:ma',4,
     &       MA_ERR)
c
      cphf_poliz = .true.
c
      end
