C
C $Id: hess_cphf.F,v 1.4 2000-05-30 20:56:43 d3e129 Exp $
C
      subroutine hess_cphf(rtdb)
c
c  Solves the CPHF equations to get the responses
c
      implicit none
c
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "nwc_const.fh"
c
      integer rtdb 
      integer i, ij, jnd
      integer nocc, nc, nxyz, num, num2, norb, nrot, nc3 
      integer num3, norb2, nocc2, nbf2, ind, ndns, in, nfock
      integer l_evals, k_evals, l_occ, k_occ
      integer g_vecs(2), g_modx
*not passed cannot be dynamic      integer g_wax(ncent*3) 
*not passed cannot be dynamic      integer g_aa(ncent*3)
      integer g_wax(nw_max_atom*3)
      integer g_aa(nw_max_atom*3)
      integer g_tmp, g_sol
      character*256 cphf_rhs, cphf_sol
      character*255 movecs
      logical hf_job, dbug, out
      double precision jfact, kfact, tol2e
      double precision pt5, one, two, four, eight
      data pt5,one,two         /0.5d+00, 1.0d+00, 2.0d+00/
      data four, eight         /4.0d+00, 8.0d+00/
c
      logical  cphf, file_write_ga, file_read_ga, movecs_read
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      external cphf, file_write_ga, file_read_ga, movecs_read
c
c call terminate since cphf will call it's own integral initialization
c
      call schwarz_tidy()    ! need to get rid of temp buffers
      call intdd_terminate()
c
      hf_job = theory .eq. 'HF'
c
      write(6,*) 'k_exy at entry ',k_exy
      call flush(6)
      dbug     = .false.
      out      = .false.
      out      = out.or.dbug
      if(out) then
         write(6,*) 'in hnd_cphf'
      endif

c     ----- get nocc -----
c     The code that depends on nocc will have to be changed to use nalpha and
c     nbeta to handle UHF case (nalpha = nclosed in RHF case)
c 

      nocc = nclosed

c     ----- define parameters -----
 
      nc   = n3xyz
      nxyz = nc              ! substitue nc or nxyz for one or the other (later)
      num  = nbf
      norb = num
      nrot = nocc * nvirt
      nc3  = nc * nc
      num2 =(num*(num+1))/2
      num3 = num* num
      norb2= norb*norb
      nocc2= nocc*nocc
c
c may need to pick up hessian off of disk, but for now assuming it is 
c in memory
c

c
c Get MO vectors, first figure out name of movecs file, get memory, and
c then get actual data
c
      if (theory.eq.'dft') then
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs))
     *    call errquit('onedd_cont: DFT MO vectors not defined',555)
      else
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs))
     $    call errquit('gradients: SCF MO vectors not defined',0)
      endif
      call util_file_name_resolve(movecs,.false.)
c
c  Get memory, note that the occ space is reused for
c  UHF since it isn't used in the actual calculation
c
      nbf2 = nbf
      if (scftype .eq. 'UHF') nbf2 = nbf + nbf
      g_vecs(1) = ga_create_atom_blocked(geom, basis, 'mo vectors')
c     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('hess_cphf: could not allocate l_evals',nbf)
c     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf,'occ. numbers', l_occ, k_occ))
     $     call errquit('hess_cphf: could not allocate l_occ',nbf)
c
      if (.not. movecs_read (movecs,1,dbl_mb(k_occ),dbl_mb(k_evals), 
     $     g_vecs(1)))
     $     call errquit('hess_cphf: could not read mo vectors', 110)
 
      if (scftype .eq. 'UHF') then
         g_vecs(2) = ga_create_atom_blocked(geom, basis, 
     $       'uhf mo vectors')
         if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), 
     $        dbl_mb(k_evals+nbf+1), g_vecs(2))) call errquit
     $        ('hess_cphf: could not read beta UHF MOs',0)
      end if
c
c print out some debug info, note that occ and evals are for the
c beta orbitals if UHF
c
       if (out.and.ga_nodeid().eq.0) then
           call hnd_prsq(dbl_mb(k_occ),1,num,num)
           call hnd_prsq(dbl_mb(k_evals),1,num,num)
           call ga_print(g_vecs(1))
           if (scftype .eq. 'UHF') call ga_print(g_vecs(2))
       endif
c
c give back unnecessary memory
c
      if (.not.ma_pop_stack(l_occ)) call errquit
     &  ('hess_cphf: cannot pop stack l_occ',555)
c
c  Transform the current AO based Fock derivatives to the MO basis
c
      g_modx =ga_create_atom_blocked ( geom, basis, 'mo df' )
c
      ind = 0
      do ndns = 1, ndens
       do in =1,nc
        ind = ind + 1
        call ga_zero(g_modx)
        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(ndns), 1, num,   1, num,
     $                         g_rhs(ind),   1, num,   1, num,
     $                         g_modx,       1, num,   1, num)
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modx,       1, num,  1, num,
     $                         g_vecs(ndns), 1, num,  1, num,
     $                         g_rhs(ind),   1, num,  1, num )
c
        if (out) then
          write(6,*) 'Transformed RHS before symmetrization'
          call ga_print(g_rhs(ind))
        endif
        call ga_symmetrize(g_rhs(ind))
        if (out) then
          write(6,*) 'Transformed RHS after symmetrization'
          call ga_print(g_rhs(ind))
        endif
c
        call ga_dscal(g_rhs(ind),four)
c
c  Pull out contributions to wax and aa, first getting memory for each.
c
        if(.not.ga_create(mt_dbl,nocc,nvirt,'g_wax',0,0,g_wax(ind)))
     $     call errquit('cphf_hess:could not alloc wax array',555) 
        if(.not.ga_create(mt_dbl,nocc,nocc,'g_aa',0,0,g_aa(ind)))
     $     call errquit('cphf_hess:could not alloc aa array',555) 
c
c Need to look up the exact calls, but should be something like
c
        write(6,*) 'before the wax copy_patch ',ind
        call flush(6)
        call ga_copy_patch('N',g_rhs(ind),1,nocc,nocc+1,num,
     $     g_wax(ind),1,nocc,1,nvirt)
        write(6,*) 'before the aa copy_patch ',ind
        call flush(6)
        call ga_copy_patch('N',g_rhs(ind),1,nocc,1,nocc,
     $     g_aa(ind),1,nocc,1,nocc)
        write(6,*) 'after the aa copy_patch ',ind
        call flush(6)
        
       enddo
      enddo
c
c This loop and the one before can probably be put together in one loop.
c Will reuse g_rhs for the ds/dx term.
c Transform ds/dx term into MO basis.
c First read in ds/dx to g_rhs(ind)
c
      call sx_read(g_rhs,nc,num,sx_dra_handle)
c
c This is not correct for UHF (need to figure out which orbitals to
c use for this case)
c
      ind = 0
      do ndns = 1, ndens
       do in =1,nc
         ind = ind + 1
         call ga_zero(g_modx)

         call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(ndns), 1, num,   1, num,
     $                         g_rhs(ind),   1, num,   1, num,
     $                         g_modx,       1, num,   1, num)
         call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modx,       1, num,  1, num,
     $                         g_vecs(ndns), 1, num,  1, num,
     $                         g_rhs(ind),   1, num,  1, num )

         call ga_symmetrize(g_rhs(ind))
         call ga_dscal(g_rhs(ind),-pt5)
c
c get contribution to wax, need to put into ga language 
c try to do this a row at a time as a ga_add_patch
c
         do i = 1, nocc
           write(6,*) 'index ',(ndns-1)*nbf+i-1
           call ga_add_patch(one,g_wax(ind),i,i,1,nvirt,
     $       eight*dbl_mb(k_evals+(ndns-1)*nbf+i-1),g_rhs(ind),i,i,
     $       nocc+1, num,g_wax(ind),i,i,1,nvirt)
         enddo
c
c        g_wax(ind,i,a) = g_wax(ind,i,a) + g_rhs(ind,i,a)*eight
c    $        *dbl_mb(k_evals,i)
       enddo
      enddo
c
c
c Do wxy = wxy + 8*sum(i)sum(j)[(ei+ej)*Ta((i-1*nocc)+j,x))*Ta((i-1*nocc)+j,y)]
c This is the hondo equivalent of cpwab1 [ Ta is part of g_rhs, e is 
c dbl_mb(k_evals) and wxy is dbl_mb(k_exy)
c   Benny is doing this routine.
c

      if (.not.ma_pop_stack(l_evals)) call errquit
     &  ('hess_cphf: cannot pop stack l_evals',555)
c
c Do modified fock build that contracts integrals with Ta and then
c adds to g_aa and g_wax.  Ta needs to be transformed back into the AO
c basis, combined with the integrals, and then transformed back into the
c MO basis.
c
      do ind = 1, nc
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_vecs(1),    1, num,   1, nocc,
     $                         g_rhs(ind),   1, nocc,  1, nocc,
     $                         g_modx,       1, num,   1, nocc)
       call ga_matmul_patch ('n','t', 1.D0, 0.D0,
     $                         g_modx,       1, num,  1, nocc,
     $                         g_vecs(1),    1, nocc, 1, num,
     $                         g_rhs(ind),   1, num,  1, num )
c
c Set up some parameters for the fock_2e call.  Note that I am not
c sure if jfact and kfact should stay the same for DFT!
c
        nfock = 1
        jfact =  8.0d+00
        kfact = -4.0d+00
        if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) then
          tol2e = 1.0d-9
        endif 

        call ga_copy(g_rhs(ind),g_modx)  ! copy Ta to g_modx
        call ga_zero(g_rhs(ind))   ! zeroed to get the "new" fock matrices

        call int_init(rtdb,1,basis)
        call schwarz_init(geom,basis)
        call scf_get_fock_param(rtdb,tol2e)
        call fock_2e(geom,basis,nfock,jfact,kfact,tol2e,skelsym,
     &              g_modx,g_rhs(ind), .false.)
        call fock_2e_tidy(rtdb)
        call int_terminate()
        call schwarz_tidy()

        call ga_zero(g_modx) ! prepare for transformation

        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs(ind),   1, num,    1, num,
     $                         g_modx,       1, num,    1, num )
 
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modx,       1, num,    1, num,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs(ind),   1, num,    1, num ) 

c
c Now form the Wax and Aa contributions
c
        call ga_add_patch(one,g_wax(ind),1,nocc,1,nvirt,
     $                    two,g_rhs(ind),1,nocc,nocc+1,nocc+nvirt,
     $                        g_wax(ind),1,nocc,1,nvirt)

        call ga_add_patch(one,g_aa(ind), 1,nocc,1,nocc,
     $                    one,g_rhs(ind),1,nocc,1,nocc,
     $                        g_aa(ind), 1,nocc,1,nocc)
      enddo
c
c Do wxy = wxy + sum(i)sum(j){Ta([i-1]*nocc+j,x)*AA([i-1]*nocc+j,y) +
c Ta([i-1]*nocc+j,y)*AA([i-1]*nocc+j,x)
c This is the hondo equivalent of cpwab2
c
      write(6,*) 'k_exy at Ta*AA sum ',k_exy
      call flush(6)
      ij = 0
      do jnd = 1, nc
        do ind = 1, nc
          ij = ij + 1
          dbl_mb(k_exy+ij) = 
     $      ga_ddot_patch(g_rhs(ind),'N',1,nocc,1,nocc,
     $                    g_aa(jnd), 'N',1,nocc,1,nocc) +
     $      ga_ddot_patch(g_rhs(jnd),'N',1,nocc,1,nocc,
     $                    g_aa(ind), 'N',1,nocc,1,nocc)
        enddo
      enddo
      write(6,*) 'k_exy after Ta*AA sum ',k_exy
      call flush(6)
c
c get rid of some extra memory
c
      if (.not.ga_destroy(g_modx))
     $  call errquit('hess_cphf: cannot destroy g_modf',555)
      do ind = 1, nc*ndens
        if (.not.ga_destroy(g_aa(ind))) call errquit
     $    ('cannot destroy g_aa ',ind)
        if (.not.ga_destroy(g_rhs(ind))) call errquit
     $    ('cannot destroy g_rhs ',ind)
      enddo
      write(6,*) 'k_exy after ga destroys ',k_exy
      call flush(6)

c
c Now actually set up the CPHF matrices to run through solutions one at a
c time.
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
 
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)

      write(6,*) 'k_exy after getting file names ',k_exy
      call flush(6)
      if (.not.ga_create(mt_dbl,nvirt*nocc,1,'cphf_tmp',0,0,g_tmp))
     $  call errquit('hess_cphf: could not create g_tmp',555)
        if(.not.ga_create(mt_dbl,nocc,nvirt,'g_sol',0,0,g_sol))
     $     call errquit('hess_cphf:could not alloc sol array',555) 

      do ind = 1, nc
        write(6,*) 'before the tmp copy_patch ',ind
        call flush(6)
        call ga_copy_patch('T',g_wax(ind),1,nocc,1,nvirt,
     $    g_tmp,1,nvirt*nocc,1,1)   ! put into format for current cphf

        if(.not.file_write_ga(cphf_rhs,g_tmp)) call errquit
     $    ('hess_cphf:could not write cphf_rhs',ind)
c
c Actually run the CPHF on one rhs
c
      write(6,*) 'k_exy before cphf ',k_exy,ind
      call flush(6)
        if (.not.cphf(rtdb)) call errquit
     $    ('hess_cphf: failure in cphf for index ',ind)
      write(6,*) 'k_exy after cphf ',k_exy,ind
      call flush(6)
c
c Get the solution off of disk, put it into appropriate shape and then 
c combine it with the wax matrices to get final contribution to hessian.
c
        if(.not.file_read_ga(cphf_sol,g_tmp)) call errquit
     $    ('hess_cphf:could not read cphf_rhs',ind)
        write(6,*) 'before the sol copy_patch ',ind
        call flush(6)
        call ga_copy_patch('T',g_tmp,1,nvirt*nocc,1,1,
     $    g_sol,1,nocc,1,nvirt)
      write(6,*) 'k_exy before hess_pxfy ',k_exy
      call flush(6)
        call hess_pxfy(dbl_mb(k_exy),g_sol,g_wax,log_mb(k_act),ncent,
     $    ind, hf_job)
      enddo
c
c Get rid of some extra memory
c
      do ind = 1, ndens
        if (.not.ga_destroy(g_vecs(ind))) call errquit
     $    ('cannot destroy g_vecs ',ind)
      enddo
      do ind = 1, ndens*nc
        if (.not.ga_destroy(g_wax(ind))) call errquit
     $    ('cannot destroy g_wax ',ind)
      enddo
      if (.not.ga_destroy(g_tmp)) call errquit
     $  ('cannot destroy g_tmp ',555)
      if (.not.ga_destroy(g_sol)) call errquit
     $  ('cannot destroy g_sol ',555)
c
      return
      end
