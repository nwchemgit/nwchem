C
C $Id: hess_cphf.F,v 1.19 2000-09-19 21:59:52 bjohnson Exp $
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
#include "msgids.fh"
c
      integer rtdb 
      integer i, ij, jnd
      integer nocc(2), nc, nxyz, num, num2, norb, nc3 
      integer num3, norb2, nbf2, ind, ndns, in, nfock
      integer l_evals, k_evals, l_occ, k_occ
      integer g_vecs(2), g_modx
      integer g_wax(nw_max_atom*3)
      integer g_aa(nw_max_atom*3)
      integer g_tmp, g_tmp2, g_sol
      character*256 cphf_rhs, cphf_sol
      character*255 movecs
      logical hf_job, dbug, out
      double precision jfact, kfact, tol2e, dmax
      double precision pt5, one, two, four, eight
      data pt5,one,two         /0.5d+00, 1.0d+00, 2.0d+00/
      data four, eight         /4.0d+00, 8.0d+00/
c
      logical  cphf, file_write_ga, file_read_ga, movecs_read
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      external cphf, file_write_ga, file_read_ga, movecs_read
c
      hf_job = theory .eq. 'hf'
c
      write(6,*) 'k_exy at entry ',k_exy
      call flush(6)
      dbug     = .false.
      out      = .true.
      out      = out.or.dbug
      if(out) then
         write(6,*) 'in hess_cphf'
      endif

c     ----- get nocc -----
c 

      if (ndens.eq.1) then  ! RHF case
       nocc(1) = nclosed
      else                  ! UHF case
       nocc(1) = nalpha
       nocc(2) = nbeta
      endif

c     ----- define parameters -----
 
      nc   = n3xyz
      nxyz = nc          ! substitue nc or nxyz for one or the other (later)
      num  = nbf
      norb = num
      nc3  = nc * nc
      num2 =(num*(num+1))/2
      num3 = num* num
      norb2= norb*norb
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
     $        dbl_mb(k_evals+nbf), g_vecs(2))) call errquit
     $        ('hess_cphf: could not read beta UHF MOs',0)
      end if
c
c print out some debug info, note that occ and evals are for the
c beta orbitals if UHF
c
       if (out.and.ga_nodeid().eq.0) 
     .     write(6,*) 'occupations, eigenvalues and vectors'
       call hnd_prsq(dbl_mb(k_occ),1,num,num)
       call hnd_prsq(dbl_mb(k_evals),1,num,num)
       call ga_print(g_vecs(1))
       if (scftype .eq. 'UHF') call ga_print(g_vecs(2))
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
        if (out) then
          write(6,*) 'RHS before MO transformation'
          call ga_print(g_rhs(ind))
        endif
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
        write(6,*) 'nocc and nvirt ',nocc(ndns),nvirt(ndns)
        call flush(6)
        if(.not.ga_create(MT_DBL,nocc(ndns),nvirt(ndns),'g_wax',0,0,
     $     g_wax(ind)))
     $     call errquit('cphf_hess:could not alloc wax array',555) 
        if(.not.ga_create(MT_DBL,nocc(ndns),nocc(ndns),'g_aa',0,0,
     $     g_aa(ind)))
     $     call errquit('cphf_hess:could not alloc aa array',555) 
c
c Need to look up the exact calls, but should be something like
c
        write(6,*) 'before the wax copy_patch ',ind
        call flush(6)
        call ga_copy_patch('N',g_rhs(ind),1,nocc(ndns),nocc(ndns)+1,
     $     num, g_wax(ind),1,nocc(ndns),1,nvirt(ndns))
        write(6,*) 'before the aa copy_patch ',ind
        call flush(6)
        call ga_copy_patch('N',g_rhs(ind),1,nocc(ndns),1,nocc(ndns),
     $     g_aa(ind),1,nocc(ndns),1,nocc(ndns))
        write(6,*) 'after the aa copy_patch ',ind
        call flush(6)
        if (out) then
          write(6,*) 'g_wax(ind) ',ind
          call ga_print(g_wax(ind))
          write(6,*) 'g_aa(ind) ',ind
          call ga_print(g_aa(ind))
        endif
        
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
         write(6,*) 'Ta after transform'
         call ga_print(g_rhs(ind))
         call ga_dscal(g_rhs(ind),-pt5)
         write(6,*) 'Ta after scaling'
         call ga_print(g_rhs(ind))
c
c get contribution to wax, need to put into ga language 
c try to do this a row at a time as a ga_add_patch
c
         do i = 1, nocc(ndns)
           write(6,*) 'index ',(ndns-1)*nbf+i
           call ga_add_patch(one,g_wax(ind),i,i,1,nvirt(ndns),
     $       eight*dbl_mb(k_evals+(ndns-1)*nbf+i-1),g_rhs(ind),i,i,
     $       nocc(ndns)+1, num,g_wax(ind),i,i,1,nvirt(ndns))
         enddo
         if (out) then
           write(6,*) 'g_wax after T addition'
           call ga_print(g_wax(ind))
         endif
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
c
      do ndns = 1, ndens
      write(6,*) 'nocc and ncent before hess_sxsy ',nocc(ndns),ncent
      call flush(6)
      call hess_sxsy(dbl_mb(k_exy),dbl_mb(k_evals),g_rhs,
     &    log_mb(k_act),ncent,nocc(ndns),hf_job)
      enddo
c
c Print intermediate result
c
      write(6,*) 'Hessian after Ta*Tb contribution'
      call flush(6)
      call hnd_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
c
c Do final summation of the hessian to get result
c
      call ga_dgop(msg_hess_exy,dbl_mb(k_exy),n3xyz2,'+') 
c
c  Get rid of some memory
c
      if (.not.ma_pop_stack(l_evals)) call errquit
     &  ('hess_cphf: cannot pop stack l_evals',555)
c
c Do modified fock build that contracts integrals with Ta and then
c adds to g_aa and g_wax.  Ta needs to be transformed back into the AO
c basis, combined with the integrals, and then transformed back into the
c MO basis.
c
      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) then
        tol2e = 1.0d-9
      endif 
      g_tmp =ga_create_atom_blocked ( geom, basis, 'temp' )
      do ndns = 1, ndens
      do ind = 1, nc
        call ga_copy(g_rhs(ind),g_tmp)  ! copy Ta to temp
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                g_vecs(1),    1, num,         1, nocc(ndns),
     $                g_rhs(ind),   1, nocc(ndns),  1, nocc(ndns),
     $                g_modx,       1, num,         1, nocc(ndns))
       call ga_matmul_patch ('n','t', 1.D0, 0.D0,
     $                g_modx,       1, num,         1, nocc(ndns),
     $                g_vecs(1),    1, nocc(ndns),  1, num,
     $                g_rhs(ind),   1, num,         1, num )
c
c Set up some parameters for the fock_2e call.  Note that I am not
c sure if jfact and kfact should stay the same for DFT!  This needs to
c be dependent on the jfac and kfac in hess_info.fh.
c
        call ga_maxelt(g_rhs(ind),dmax)
        if (dmax.gt.tol2e) then  ! the value of tol2e may not be the best here
        nfock = 1
        jfact =  8.0d+00*jfac(ndns)
        kfact = -4.0d+00*kfac(ndns)

        call ga_copy(g_rhs(ind),g_modx)  ! copy Ta to g_modx
        call ga_zero(g_rhs(ind))   ! zeroed to get the "new" fock matrices

        if (out) then
          write(6,*) 'Before the fock_2e call for ind ',ind
          write(6,*) 'nfock, jfact, kfact ',nfock, jfact, kfact
          write(6,*) 'density like matrix before fock_2e'
          call flush(6)
          call ga_print(g_modx)
        endif
        call int_init(rtdb,1,basis)
        call schwarz_init(geom,basis)
        call scf_get_fock_param(rtdb,tol2e) ! need to change for DFT?
c     Activate XC contribution to OO derivative Fock matrix
        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $       call errquit('hess_cphf: rtdb_put of xc_active failed',0)
        if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $       call errquit('hess_cphf: rtdb_put of calc_type failed',0)
        if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $       call errquit('hess_cphf: rtdb_put of j_derfit failed',0)
        call fock_2e(geom,basis,nfock,jfact,kfact,tol2e,skelsym,
     &              g_modx,g_rhs(ind), .false.)
        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('hess_cphf: rtdb_put of xc_active failed',1)
        if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $       call errquit('hess_cphf: rtdb_put of calc_type failed',0)
        call fock_2e_tidy(rtdb)
        call int_terminate()
        call schwarz_tidy()

       write(6,*) 'fock matrix after fock_2e'
       call flush(6)
       call ga_print(g_rhs(ind))
        call ga_zero(g_modx) ! prepare for transformation

        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs(ind),   1, num,    1, num,
     $                         g_modx,       1, num,    1, num )
 
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modx,       1, num,    1, num,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs(ind),   1, num,    1, num ) 
        if (out) then
          write(6,*) 'after transformation back to MO basis'
          call ga_print(g_rhs(ind))
        endif

c
c Now form the Wax and Aa contributions
c
        call ga_add_patch(one,g_wax(ind),1,nocc(ndns),1,nvirt(ndns),
     $                    two,g_rhs(ind),1,nocc(ndns),nocc(ndns)+1,
     $                        nocc(ndns)+nvirt(ndns),
     $                        g_wax(ind),1,nocc(ndns),1,nvirt(ndns))

        call ga_add_patch(one,g_aa(ind), 1,nocc(ndns),1,nocc(ndns),
     $                    one,g_rhs(ind),1,nocc(ndns),1,nocc(ndns),
     $                        g_aa(ind), 1,nocc(ndns),1,nocc(ndns))
        endif
        call ga_copy(g_tmp,g_rhs(ind))  ! copy temp to Ta
      enddo
      enddo
      if (.not.ga_destroy(g_tmp))
     $  call errquit('hess_cphf: cannot destroy g_tmp',555)
c
c Print out intermediate results
c
      write(6,*) 'WAX after fock contribution'
      call flush(6)
      do ind = 1, nc
        call ga_print(g_wax(ind))
      enddo
      write(6,*) 'AA after fock contribution'
      call flush(6)
      do ind = 1, nc
        call ga_print(g_aa(ind))
      enddo
      write(6,*) 'g_rhs after fock contribution'
      call flush(6)
      do ind = 1, nc
        call ga_print(g_rhs(ind))
      enddo
c
c Do wxy = wxy + sum(i)sum(j){Ta([i-1]*nocc+j,x)*AA([i-1]*nocc+j,y) +
c Ta([i-1]*nocc+j,y)*AA([i-1]*nocc+j,x)
c This is the hondo equivalent of cpwab2
c
      write(6,*) 'k_exy at Ta*AA sum ',k_exy
      call flush(6)
      do ndns = 1, ndens
      ij = 0
      do jnd = 1, nc
        do ind = 1, nc
          dbl_mb(k_exy+ij) = 
     $      dbl_mb(k_exy+ij) +
     $      ga_ddot_patch(g_rhs(ind),'N',1,nocc(ndns),1,nocc(ndns),
     $                    g_aa(jnd), 'N',1,nocc(ndns),1,nocc(ndns)) +
     $      ga_ddot_patch(g_rhs(jnd),'N',1,nocc(ndns),1,nocc(ndns),
     $                    g_aa(ind), 'N',1,nocc(ndns),1,nocc(ndns))
          ij = ij + 1
        enddo
      enddo
      enddo
      write(6,*) 'k_exy after Ta*AA sum ',k_exy
      call flush(6)
c
c Print intermediate result
c
      write(6,*) 'Hessian after Ta*AA contribution'
      call flush(6)
      call hnd_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
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

      do ndns = 1, ndens
        if (.not.ga_create(MT_DBL,nvirt(ndns)*nocc(ndns),1,'cphf_tmp',
     $    0,0,g_tmp))
     $    call errquit('hess_cphf: could not create g_tmp',555)
        if(.not.ga_create(MT_DBL,nvirt(ndns),nocc(ndns),'g_tmp2',
     $     0,0,g_tmp2))
     $     call errquit('hess_cphf:could not alloc tmp2 array',555) 
        if(.not.ga_create(MT_DBL,nocc(ndns),nvirt(ndns),'g_sol',
     $     0,0,g_sol))
     $     call errquit('hess_cphf:could not alloc sol array',555) 
      do ind = 1, nc
        write(6,*) 'before the tmp copy_patch ',ind
        call flush(6)
        call ga_copy_patch('T',g_wax(ind),1,nocc(ndns),1,nvirt(ndns),
     $    g_tmp,1,nvirt(ndns)*nocc(ndns),1,1)   ! put into format for current cphf

        write(6,*) 'g_wax_tr ',ind
        call flush(6)
        call ga_print(g_tmp)
        call ga_maxelt(g_tmp,dmax)
        if (dmax.gt.tol2e) then  ! the value of tol2e may not be the best here
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
        call ga_print(g_tmp)
        call ga_copy_patch('N',g_tmp,1,nvirt(ndns)*nocc(ndns),1,1,
     $    g_tmp2,1,nvirt(ndns),1,nocc(ndns))
      write(6,*) 'g_tmp2 ',ind
      call flush(6)
      call ga_print(g_tmp2)
        call ga_transpose(g_tmp2,g_sol)
      write(6,*) 'k_exy before hess_pxfy ',k_exy
      write(6,*) 'g_sol after CPHF for ind ',ind
      call flush(6)
      call ga_print(g_sol)
        call hess_pxfy(dbl_mb(k_exy),g_sol,g_wax,log_mb(k_act),ncent,
     $    ind, hf_job)
        endif
      enddo
        if (.not.ga_destroy(g_tmp)) call errquit
     $    ('cannot destroy g_tmp ',555)
        if (.not.ga_destroy(g_tmp2)) call errquit
     $    ('cannot destroy g_tmp2 ',555)
        if (.not.ga_destroy(g_sol)) call errquit
     $    ('cannot destroy g_sol ',555)
      enddo
c
c Print intermediate result
c
      if(ga_nodeid().eq.0) then
         write(6,*) 'Hessian after CPHF contribution'
         call flush(6)
         call hnd_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
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
c
      return
      end
