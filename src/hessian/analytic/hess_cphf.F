C
C $Id: hess_cphf.F,v 1.29 2001-11-07 21:29:04 windus Exp $
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
#include "nwc_const.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer rtdb 
      integer i, j, ij, jnd, index
      integer nocc(2), nc, nxyz, num, num2, norb, nc3 
      integer num3, norb2, nbf2, ind, ndns, in, nfock
      integer l_evals, k_evals, l_occ, k_occ
      integer g_vecs(2), g_modx
      integer g_wax(nw_max_atom*3)
      integer g_aa(nw_max_atom*3)
      integer g_tmp, g_tmp2, g_sol
c temp variables
      integer g_rhs_x(200), alo(3), ahi(3), blo(2), bhi(2)
c
      character*256 cphf_rhs, cphf_sol
      character*255 movecs
      logical hf_job, oprint, ocphfprint
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
c Set up printing
c
      oprint = util_print('hess_follow',print_high)
      oprint = oprint .and. (ga_nodeid().eq.0)
      ocphfprint = util_print('cphf_cont',print_debug)
      if (oprint) write(LuOut,*) 'HESS: setting up cphf'
c
c     ----- get nocc -----
c 

      if (ndens.eq.1) then  ! RHF case
       nocc(1) = nclosed
      else                  ! UHF case
       nocc(1) = nalpha
       nocc(2) = nbeta
      endif
c
c     This is temporary until I can get the code working with multiple
c     dimensional GAs: copy g_rhs into g_rhs_x which is an array of two-
c     dimensional GAs.
c
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      j = 0
      do i = 1, ncent
         j = j + 1
         alo(1) = j
         ahi(1) = j
         g_rhs_x(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a x')
         call nga_copy_patch('N',g_rhs,alo,ahi,g_rhs_x(j),blo,bhi)
         j = j + 1
         alo(1) = j
         ahi(1) = j
         g_rhs_x(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a y')
         call nga_copy_patch('N',g_rhs,alo,ahi,g_rhs_x(j),blo,bhi)
         j = j + 1
         alo(1) = j
         ahi(1) = j
         g_rhs_x(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test a z')
         call nga_copy_patch('N',g_rhs,alo,ahi,g_rhs_x(j),blo,bhi)
      enddo
      if (ndens.eq.2) then
      do i = 1, ncent
         j = j + 1
         alo(1) = j
         ahi(1) = j
         g_rhs_x(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test b x')
         call nga_copy_patch('N',g_rhs,alo,ahi,g_rhs_x(j),blo,bhi)
         j = j + 1
         alo(1) = j
         ahi(1) = j
         g_rhs_x(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test b y')
         call nga_copy_patch('N',g_rhs,alo,ahi,g_rhs_x(j),blo,bhi)
         j = j + 1
         alo(1) = j
         ahi(1) = j
         g_rhs_x(j) = ga_create_atom_blocked
     &        (geom, basis, 'CPKS RHS test b z')
         call nga_copy_patch('N',g_rhs,alo,ahi,g_rhs_x(j),blo,bhi)
      enddo
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
c print out some debug info, note that occ is for the
c beta orbitals if UHF and the evals are for the alpha orbitals
c
      if (ocphfprint) then
         if(ga_nodeid().eq.0) then
         write(LuOut,*)
         write(LuOut,*) 'occupations, eigenvalues and vectors'
         call hnd_prsq(dbl_mb(k_occ),1,num,num)
         write(LuOut,*)
         call hnd_prsq(dbl_mb(k_evals),1,num,num)
         write(LuOut,*)
         endif
         call ga_print(g_vecs(1))
         if (scftype .eq. 'UHF') call ga_print(g_vecs(2))
         if(ga_nodeid().eq.0) write(LuOut,*)
      endif
c
c give back unnecessary memory
c
      if (.not.ma_pop_stack(l_occ)) call errquit
     &  ('hess_cphf: cannot pop stack l_occ',555)
c
      g_modx =ga_create_atom_blocked ( geom, basis, 'mo df' )
      ind = 0
      do ndns = 1, ndens
       do in = 1, nc
        ind = ind + 1
        if(.not.ga_create(MT_DBL,nocc(ndns),nvirt(ndns),'g_wax',0,0,
     $     g_wax(ind)))
     $     call errquit('cphf_hess:could not alloc wax array',555) 
        if(.not.ga_create(MT_DBL,nocc(ndns),nocc(ndns),'g_aa',0,0,
     $     g_aa(ind)))
     $     call errquit('cphf_hess:could not alloc aa array',555) 
       enddo
      enddo
c
c  Transform the current AO based Fock derivatives to the MO basis,
c  scale by four, and populate g_aa and g_wax 
c
      call hess_aotomo(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      ind = 0
      do ndns = 1, ndens
       do in = 1, nc
        ind = ind + 1
        call ga_dscal(g_rhs_x(ind),four)
       enddo
      enddo
      call hess_getaawax(g_rhs_x, g_aa, g_wax, nocc, nvirt,
     $          ndens, n3xyz, num, ocphfprint)
c
c We are now done with the Fock derivatives.
c Will reuse g_rhs for the ds/dx term.
c Transform ds/dx term into MO basis and scale by -0.5.
c
      call sx_read_orig(g_rhs_x,nc,num,sx_dra_handle)
c
      if (oprint) write(LuOut,*) 'HESS: Ta contribution to Wax'
      call hess_aotomo(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      ind = 0
      do ndns = 1, ndens
       do in =1,nc
         ind = ind + 1
         call ga_dscal(g_rhs_x(ind),-pt5)
         if (ocphfprint) then
           if(ga_nodeid().eq.0)write(LuOut,*)
           if(ga_nodeid().eq.0)write(LuOut,*) 'Ta after scaling',ind
           call ga_print(g_rhs_x(ind))
         endif
       enddo
      enddo
c
c get contribution to wax
c
      call hess_tuwax(g_rhs_x, g_wax, dbl_mb(k_evals), nocc, nvirt,
     $          ocphfprint, ndens, n3xyz, num)
c
c
c Do wxy = wxy + 8*sum(i)sum(j)[(ei+ej)*Ta((i-1*nocc)+j,x))*Ta((i-1*nocc)+j,y)]
c This is the hondo equivalent of cpwab1 [ Ta is part of g_rhs, e is 
c dbl_mb(k_evals) and wxy is dbl_mb(k_exy)
c
      if (oprint) write(LuOut,*) 'HESS: Ta contribution to Hessian'
      do ndns = 1, ndens
c      write(6,*) 'nocc and ncent before hess_sxsy ',nocc(ndns),ncent
      call hess_sxsy(dbl_mb(k_exy),dbl_mb(k_evals),g_rhs_x,
     &    log_mb(k_act),ncent,nocc(ndns),hf_job)
      enddo
c
c Do final summation of the hessian to get result
c
      call ga_dgop(msg_hess_exy,dbl_mb(k_exy),n3xyz2,'+') 
c
c Print intermediate result
c
      if (ocphfprint.and.ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*) 'Hessian after Ta*Tb contribution'
        call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
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
      if (oprint) write(LuOut,*) 'HESS: Modified Fock builds'
      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) 
     $  tol2e = 1.0d-9
c
      call hess_moootoao(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      call hess_modfock(g_rhs_x, g_vecs, g_modx, ndens, n3xyz,
     $                  tol2e, rtdb, geom, basis, jfac, kfac,
     $                  num, ocphfprint)
c
c Now form the Wax and Aa contributions
c
      call hess_fxwaxaa(g_rhs_x, g_wax, g_aa, nocc, nvirt,
     $                  ndens, n3xyz)
c
c Transform ds/dx term into MO basis and scale by -0.5.
c
      write(LuOut,*) 'Getting Ta again'
      call sx_read_orig(g_rhs_x,nc,num,sx_dra_handle)
c
      call hess_aotomo(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      ind = 0
      do ndns = 1, ndens
       do in =1,nc
         ind = ind + 1
         call ga_dscal(g_rhs_x(ind),-pt5)
         if (ocphfprint) then
           if(ga_nodeid().eq.0)write(LuOut,*)
           if(ga_nodeid().eq.0)write(LuOut,*) 'Ta after scaling',ind
           call ga_print(g_rhs_x(ind))
         endif
       enddo
      enddo
c
c Print out intermediate results
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'WAX after fock contribution'
          endif
        do ind = 1, nc
          call ga_print(g_wax(ind))
        enddo
        if(ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*) 'AA after fock contribution'
        endif
        do ind = 1, nc
          call ga_print(g_aa(ind))
        enddo
        if(ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*) 'g_rhs after fock contribution'
        endif
        do ind = 1, nc
          call ga_print(g_rhs_x(ind))
        enddo
      endif
c
c Do wxy = wxy + sum(i)sum(j){Ta([i-1]*nocc+j,x)*AA([i-1]*nocc+j,y) +
c Ta([i-1]*nocc+j,y)*AA([i-1]*nocc+j,x)
c
      if (oprint) write(LuOut,*) 'HESS: Ta and AA contributions'
      call hess_taaa(dbl_mb(k_exy),g_rhs_x, g_aa, nocc, 
     $               ndens, n3xyz)
c
c Print intermediate result
c
      if (ocphfprint.and.(ga_nodeid().eq.0)) then
        write(LuOut,*)
        write(LuOut,*) 'Hessian after Ta*AA contribution'
        call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
c
c get rid of some extra memory
c
      if (.not.ga_destroy(g_modx))
     $  call errquit('hess_cphf: cannot destroy g_modf',555)
      do ind = 1, nc*ndens
        if (.not.ga_destroy(g_aa(ind))) call errquit
     $    ('cannot destroy g_aa ',ind)
        if (.not.ga_destroy(g_rhs_x(ind))) call errquit
     $    ('cannot destroy g_rhs ',ind)
      enddo
c
c Now actually set up the CPHF matrices 
c
      if (oprint) write(LuOut,*) 'HESS: CPHF contributions'
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
 
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)

c
c The next line obviously won't work for UHF!
c
       if (.not.ga_create(MT_DBL,nvirt(1)*nocc(1),nc*ndens,'cphf_tmp',
     $    0,0,g_tmp))
     $    call errquit('hess_cphf: could not create g_tmp',555)
c
c Put the CPHF matrices into one GA in the correct format
c
      index = 1
      do ndns = 1, ndens
      do ind = 1, nc
        call ga_copy_patch('T',g_wax(ind),1,nocc(ndns),1,nvirt(ndns),
     $    g_tmp,1,nvirt(ndns)*nocc(ndns),index,index)   ! put into format for current cphf
        index = index + 1
      enddo
      enddo
      if (ocphfprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        call ga_print(g_tmp)
      endif
      if(.not.file_write_ga(cphf_rhs,g_tmp)) call errquit
     $  ('hess_cphf:could not write cphf_rhs',ind)
c
c Actually run the CPHF on all rhs
c
      if (.not.cphf(rtdb)) call errquit
     $  ('hess_cphf: failure in cphf for index ',ind)
      if(.not.file_read_ga(cphf_sol,g_tmp)) call errquit
     $  ('hess_cphf:could not read cphf_rhs',ind)
c
c Get the solution off of disk, put it into appropriate shape and then 
c combine it with the wax matrices to get final contribution to hessian.
c
      write(6,*) 'going into the g_tmp2 create'
      call flush(6)
      if(.not.ga_create(MT_DBL,nvirt(1),nocc(1),'g_tmp2',
     $   0,0,g_tmp2))
     $   call errquit('hess_cphf:could not alloc tmp2 array',555) 
      write(6,*) 'going into the g_sol create'
      call flush(6)
      if(.not.ga_create(MT_DBL,nocc(1),nvirt(1),'g_sol',
     $   0,0,g_sol))
     $   call errquit('hess_cphf:could not alloc sol array',555) 
      index = 1
      do ndns = 1, ndens
      do ind = 1, nc
c       if (ocphfprint) then
c         if(ga_nodeid().eq.0)write(LuOut,*)
c         if(ga_nodeid().eq.0)
c    $         write(LuOut,*) 'before the sol copy_patch ',ind
c         call ga_print(g_tmp)
c       endif
        call ga_copy_patch('N',g_tmp,1,nvirt(ndns)*nocc(ndns),index,
     $    index,g_tmp2,1,nvirt(ndns),1,nocc(ndns))
        if (ocphfprint) then
          if(ga_nodeid().eq.0)write(LuOut,*)
          if(ga_nodeid().eq.0)write(LuOut,*) 'g_tmp2 ',ind
          call ga_print(g_tmp2)
        endif
        call ga_transpose(g_tmp2,g_sol)
        if (ocphfprint) then
          if(ga_nodeid().eq.0)write(LuOut,*)
          if(ga_nodeid().eq.0)
     $         write(LuOut,*) 'g_sol after CPHF for ind ',ind
          call ga_print(g_sol)
        endif
        call hess_pxfy(dbl_mb(k_exy),g_sol,g_wax,log_mb(k_act),ncent,
     $    ind, hf_job)
        index = index + 1
      enddo
      enddo
      if (.not.ga_destroy(g_tmp)) call errquit
     $  ('cannot destroy g_tmp ',555)
      if (.not.ga_destroy(g_tmp2)) call errquit
     $  ('cannot destroy g_tmp2 ',555)
      if (.not.ga_destroy(g_sol)) call errquit
     $  ('cannot destroy g_sol ',555)
c
c Print intermediate result
c
      if(ocphfprint.and.ga_nodeid().eq.0) then
         write(LuOut,*)
         write(LuOut,*) 'Hessian after CPHF contribution'
         call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
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
c
      subroutine hess_aotomo(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphf,
     $          nocc, nvirt, ndens, n3xyz, num)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical oprint, ocphf
      integer ndens, n3xyz, num, nocc(ndens), nvirt(ndens)
      integer g_rhs_x(200), g_vecs(ndens)
      integer g_modx, ind, ndns, in
c
c  Transform the current AO based Fock derivatives to the MO basis
c
      if (oprint) write(LuOut,*) 'HESS: setting up Wax and AA'
      ind = 0
      do ndns = 1, ndens
       do in =1,n3xyz
        ind = ind + 1
        if (ocphf) then
          if(ga_nodeid().eq.0) 
     $          write(LuOut,*) 'RHS before MO transformation', ind
          call ga_print(g_rhs_x(ind))
        endif
        call ga_zero(g_modx)
        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(ndns), 1, num,   1, num,
     $                         g_rhs_x(ind),   1, num,   1, num,
     $                         g_modx,       1, num,   1, num)
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modx,       1, num,  1, num,
     $                         g_vecs(ndns), 1, num,  1, num,
     $                         g_rhs_x(ind),   1, num,  1, num )
c
        if (ocphf) then
          if(ga_nodeid().eq.0)  write(LuOut,*) 
     $          'Transformed RHS before symmetrization',ind
          call ga_print(g_rhs_x(ind))
        endif
c        call ga_symmetrize(g_rhs_x(ind))
        if (ocphf) then
          if(ga_nodeid().eq.0) write(LuOut,*) 
     $          'Transformed RHS after symmetrization',ind
          call ga_print(g_rhs_x(ind))
        endif
c
       enddo
      enddo
c
      return
      end
c
      subroutine hess_moootoao(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphf,
     $          nocc, nvirt, ndens, n3xyz, num)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical oprint, ocphf
      integer ndens, n3xyz, num, nocc(ndens), nvirt(ndens)
      integer g_rhs_x(200), g_vecs(ndens)
      integer g_modx, ind, ndns, in
c
c  Transform the current MO based matrices to the AO basis
c
      ind = 0
      do ndns = 1, ndens
       do in =1,n3xyz
        ind = ind + 1
        if (ocphf) then
          if(ga_nodeid().eq.0) 
     $          write(LuOut,*) 'RHS before AO transformation', ind
          call ga_print(g_rhs_x(ind))
        endif
        call ga_zero(g_modx)
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                g_vecs(ndns), 1, num,         1, nocc(ndns),
     $                g_rhs_x(ind), 1, nocc(ndns),  1, nocc(ndns),
     $                g_modx,       1, num,         1, nocc(ndns))
       call ga_matmul_patch ('n','t', 1.D0, 0.D0,
     $                g_modx,       1, num,         1, nocc(ndns),
     $                g_vecs(ndns), 1, nocc(ndns),  1, num,
     $                g_rhs_x(ind),   1, num,         1, num )
c
        if (ocphf) then
          if(ga_nodeid().eq.0) write(LuOut,*) 
     $          'Transformed RHS ',ind
          call ga_print(g_rhs_x(ind))
        endif
c
       enddo
      enddo
c
      return
      end
c
      subroutine hess_getaawax(g_rhs_x, g_aa, g_wax, nocc, nvirt,
     $          ndens, n3xyz, num, ocphf)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical ocphf
      integer ndens, n3xyz, nocc(ndens), nvirt(ndens), num
      integer g_rhs_x(200), g_aa(200), g_wax(200)
      integer ind, ndns, in
c
      ind = 0
      do ndns = 1, ndens
       do in =1,n3xyz
        ind = ind + 1
        call ga_copy_patch('N',g_rhs_x(ind),1,nocc(ndns),nocc(ndns)+1,
     $     num, g_wax(ind),1,nocc(ndns),1,nvirt(ndns))
        call ga_copy_patch('N',g_rhs_x(ind),1,nocc(ndns),1,nocc(ndns),
     $     g_aa(ind),1,nocc(ndns),1,nocc(ndns))
        if (ocphf) then
          if(ga_nodeid().eq.0) write(LuOut,*) 'g_wax(ind) ',ind
          call ga_print(g_wax(ind))
          if(ga_nodeid().eq.0) write(LuOut,*) 'g_aa(ind) ',ind
          call ga_print(g_aa(ind))
        endif
       enddo
      enddo
c
      return
      end
c
      subroutine hess_tuwax(g_rhs_x, g_wax, evals, nocc, nvirt,
     $          ocphfprint, ndens, n3xyz, num)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical ocphfprint
      integer ndens, n3xyz, nocc(ndens), nvirt(ndens), num
      integer g_rhs_x(200), g_wax(200)
      integer ind, ndns, in, i
      double precision evals(*)
      double precision one, eight
      data one, eight /1.0d+00, 8.0d+00/
c
      write(6,*) 'in tuwax: ndens,n3xyz,nocc(1)',ndens,n3xyz,nocc(1)
      ind = 0
      do ndns = 1, ndens
       do in =1,n3xyz
         ind = ind + 1
         do i = 1, nocc(ndns)
           call ga_add_patch(one,g_wax(ind),i,i,1,nvirt(ndns),
     $       eight*evals((ndns-1)*num+i),g_rhs_x(ind),i,i,
     $       nocc(ndns)+1, num,g_wax(ind),i,i,1,nvirt(ndns))
         enddo
         if (ocphfprint) then
           if(ga_nodeid().eq.0)write(LuOut,*) 
     $           'g_wax after T addition',ind
           call ga_print(g_wax(ind))
         endif
c
       enddo
      enddo
c
      return
      end
c
      subroutine hess_taaa(exy,g_rhs_x, g_aa, nocc, 
     $                     ndens, n3xyz)
c
      implicit none
c
#include "global.fh"
c
c Do exy = exy + sum(i)sum(j){Ta([i-1]*nocc+j,x)*AA([i-1]*nocc+j,y) +
c Ta([i-1]*nocc+j,y)*AA([i-1]*nocc+j,x)
c
      integer ndens, n3xyz, nocc(ndens)
      integer g_rhs_x(200), g_aa(200)
      double precision exy(n3xyz,n3xyz)
      integer ndns, jnd, ind
c
      do ndns = 1, ndens
      do jnd = 1, n3xyz
        do ind = 1, n3xyz
          exy(ind, jnd) = 
     $      exy(ind, jnd) +
     $      ga_ddot_patch(g_rhs_x(ind),'N',1,nocc(ndns),1,nocc(ndns),
     $                    g_aa(jnd), 'N',1,nocc(ndns),1,nocc(ndns)) +
     $      ga_ddot_patch(g_rhs_x(jnd),'N',1,nocc(ndns),1,nocc(ndns),
     $                    g_aa(ind), 'N',1,nocc(ndns),1,nocc(ndns))
        enddo
      enddo
      enddo
c
      return
      end
c
      subroutine hess_modfock(g_rhs_x, g_vecs, g_modx, ndens, n3xyz,
     $                        tol2e, rtdb, geom, basis, jfac, kfac,
     $                        num, ocphfprint)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer ndens, n3xyz, rtdb, geom, basis, num
      integer g_rhs_x(200), g_vecs(ndens), g_modx
      integer ind, ndns, in, nfock
      double precision tol2e, jfac(ndens), kfac(ndens)
      double precision jfact, kfact, dmax
      logical ocphfprint
c
      ind = 0
      do ndns = 1, ndens
c      do ind = 1, n3xyz
      do in = 1, n3xyz
         ind = ind+1
c
c Set up some parameters for the fock_2e call.  Note that I am not
c sure if jfact and kfact should stay the same for DFT!  This needs to
c be dependent on the jfac and kfac in hess_info.fh.
c
        call ga_maxelt(g_rhs_x(ind),dmax)
        if (dmax.gt.tol2e) then  ! the value of tol2e may not be the best here
        nfock = 1
        jfact =  8.0d+00*jfac(ndns)
        kfact = -4.0d+00*kfac(ndns)

        call ga_copy(g_rhs_x(ind),g_modx)  ! copy Ta to g_modx
        call ga_zero(g_rhs_x(ind))   ! zeroed to get the "new" fock matrices

        if (ocphfprint) then
           if(ga_nodeid().eq.0) then
          write(LuOut,*) 'Before the fock_2e call for ind ',ind
          write(LuOut,*) 'nfock, jfact, kfact ',nfock, jfact, kfact
          write(LuOut,*) 'density like matrix before fock_2e'
          endif
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
!        call fock_2e(geom,basis,nfock,jfact,kfact,tol2e,skelsym,
        call fock_2e(geom,basis,nfock,jfact,kfact,tol2e,.false.,
     &              g_modx,g_rhs_x(ind), .false.)
        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('hess_cphf: rtdb_put of xc_active failed',1)
        if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $       call errquit('hess_cphf: rtdb_put of calc_type failed',0)
        call int_terminate()
        call schwarz_tidy()
c       if(skelsym) call 
c    $       sym_symmetrize(geom,basis,.false.,g_rhs_x(ind))
        if (ocphfprint) then
           if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'fock matrix after fock_2e',ind
          endif
          call ga_print(g_rhs_x(ind))
        endif
        call ga_zero(g_modx) ! prepare for transformation

        write(LuOut,*) 'before first ga_matmul_patch'
        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs_x(ind),   1, num,    1, num,
     $                         g_modx,       1, num,    1, num )
 
        write(LuOut,*) 'before second ga_matmul_patch'
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_modx,       1, num,    1, num,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs_x(ind),   1, num,    1, num ) 
        else
          call ga_zero(g_rhs_x(ind))
        endif
        if (ocphfprint) then
           if(ga_nodeid().eq.0)
     $   write(LuOut,*) 'after transformation back to MO basis',ind
          call ga_print(g_rhs_x(ind))
        endif
      enddo
      enddo
c
      return
      end
c
      subroutine hess_modfock2(g_rhs_x, g_vecs, g_tmp, ndens, n3xyz,
     $                        tol2e, rtdb, geom, basis, jfac, kfac,
     $                        num, ocphfprint)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer ndens, n3xyz, rtdb, geom, basis, num
      integer g_rhs_x, g_vecs(ndens), g_tmp
      integer ind, ndns, in, nfock
      double precision tol2e, jfac(ndens), kfac(ndens)
      double precision jfact, kfact, dmax
      logical ocphfprint
c
c
c Set up some parameters for the fock_2e call.  Note that I am not
c sure if jfact and kfact should stay the same for DFT!  This needs to
c be dependent on the jfac and kfac in hess_info.fh.
c
      call ga_maxelt(g_rhs_x,dmax)
c     if (dmax.gt.tol2e) then  ! the value of tol2e may not be the best here
      nfock = 1
c     jfact =  8.0d+00*jfac(ndns)
c     kfact = -4.0d+00*kfac(ndns)
      jfact =  8.0d+00*jfac(1)
      kfact = -4.0d+00*kfac(1)

      call ga_copy(g_rhs_x,g_tmp)  ! copy Ta to g_tmp
      call ga_zero(g_rhs_x)   ! zeroed to get the "new" fock matrices

      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'Before the fock_2e call for ind ',ind
          write(LuOut,*) 'nfock, jfact, kfact ',nfock, jfact, kfact
          write(LuOut,*) 'density like matrix before fock_2e'
        endif
        call ga_print(g_tmp)
      endif
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call scf_get_fock_param(rtdb,tol2e) ! need to change for DFT?
c     Activate XC contribution to OO derivative Fock matrix
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of j_derfit failed',0)
c
      call fock_2e(geom,basis,nfock,jfact,kfact,tol2e,.false.,
     &            g_tmp,g_rhs_x, .false.)
c
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',1)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0)
      call int_terminate()
      call schwarz_tidy()
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'fock matrix after fock_2e',ind
        endif
        call ga_print(g_rhs_x)
      endif

      ind = 0
      do ndns = 1, ndens
      do in = 1, n3xyz
        ind = ind+1
        call ga_zero(g_tmp) ! prepare for transformation
        write(LuOut,*) 'before first ga_matmul_patch'
        call ga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs_x,   1, num,    1, num,
     $                         g_tmp,       1, num,    1, num )
 
        write(LuOut,*) 'before second ga_matmul_patch'
        call ga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                         g_tmp,       1, num,    1, num,
     $                         g_vecs(1),    1, num,    1, num,
     $                         g_rhs_x,   1, num,    1, num ) 
      enddo
      enddo
      if (ocphfprint) then
        if(ga_nodeid().eq.0)
     $    write(LuOut,*) 'after transformation back to MO basis',ind
        call ga_print(g_rhs_x)
      endif
c
      return
      end
c
      subroutine hess_fxwaxaa(g_rhs_x, g_wax, g_aa, nocc, nvirt,
     $                        ndens, n3xyz)
c
      implicit none
c
#include "global.fh"
c
c Now form the Wax and Aa contributions
c
      integer ndens, n3xyz
      integer g_rhs_x(200), g_wax(200), g_aa(200) 
      integer nocc(ndens), nvirt(ndens)
      integer ind, ndns, in
      double precision one, two
      data one,two         /1.0d+00, 2.0d+00/
c
      ind = 0
      do ndns = 1, ndens
      do in = 1, n3xyz
        ind = ind+1
        call ga_add_patch(one,g_wax(ind),1,nocc(ndns),1,nvirt(ndns),
     $                    two,g_rhs_x(ind),1,nocc(ndns),nocc(ndns)+1,
     $                        nocc(ndns)+nvirt(ndns),
     $                        g_wax(ind),1,nocc(ndns),1,nvirt(ndns))

        call ga_add_patch(one,g_aa(ind), 1,nocc(ndns),1,nocc(ndns),
     $                    one,g_rhs_x(ind),1,nocc(ndns),1,nocc(ndns),
     $                        g_aa(ind), 1,nocc(ndns),1,nocc(ndns))
      enddo
      enddo
c
      return
      end
