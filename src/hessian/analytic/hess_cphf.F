C
C $Id: hess_cphf.F,v 1.37 2001-12-28 18:45:36 windus Exp $
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
      integer ind, index
      integer nocc(2), nc, nxyz, num, num2, norb, nc3 
      integer num3, norb2, nbf2, ndns
      integer l_evals, k_evals, l_occ, k_occ
      integer g_vecs(2), g_modx
      integer g_wax
      integer g_aa
      integer g_tmp, g_tmp2, g_sol
c temp variables
      integer alo(3), ahi(3), blo(2), bhi(2)
      integer dims(3), chunk(3)
c
      character*256 cphf_rhs, cphf_sol
      character*255 movecs
      logical hf_job, oprint, ocphfprint
      double precision tol2e
      double precision pt5, four
      data pt5, four      /0.5d+00, 4.0d+00/
c
      logical  cphf2, file_write_ga, file_read_ga, movecs_read
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      external cphf2, file_write_ga, file_read_ga, movecs_read
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
c
      dims(1) = ndens*n3xyz
      dims(2) = nocc(1)   ! Note that this is correct only for ndens=1
      dims(3) = nocc(1)
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'AA matrix',chunk,
     &     g_aa))
     &   call errquit('hess_cphf: could not allocate g_aa',555)
c
      dims(3) = nvirt(1) ! Note that this is correct only for ndens=1
      if (.not. nga_create (MT_DBL, 3, dims, 'WAX matrix',chunk,
     &     g_wax))
     &   call errquit('hess_cphf: could not allocate g_wax',555)
c
c  Transform the current AO based Fock derivatives to the MO basis,
c  scale by four, and populate g_aa and g_wax 
c
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      call ga_dscal(g_rhs,four)
      call hess_getaawax(g_rhs, g_aa, g_wax, nocc, nvirt,
     $          ndens, n3xyz, num, ocphfprint)
c
c We are now done with the Fock derivatives.
c Will reuse g_rhs for the ds/dx term.
c Transform ds/dx term into MO basis and scale by -0.5.
c
      call sx_read(g_rhs,nc,num,sx_dra_handle)
c
      if (oprint) write(LuOut,*) 'HESS: Ta contribution to Wax'
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      call ga_dscal(g_rhs,-pt5)
      if (ocphfprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        if(ga_nodeid().eq.0)write(LuOut,*) 'Ta after scaling'
        call ga_print(g_rhs)
      endif
c
c get contribution to wax
c
      call hess_tuwax(g_rhs, g_wax, dbl_mb(k_evals), nocc, nvirt,
     $          ocphfprint, ndens, n3xyz, num)
c
c Do wxy = wxy + 8*sum(i)sum(j)[(ei+ej)*Ta((i-1*nocc)+j,x))*Ta((i-1*nocc)+j,y)]
c This is the hondo equivalent of cpwab1 [ Ta is part of g_rhs, e is 
c dbl_mb(k_evals) and wxy is dbl_mb(k_exy)
c
      if (oprint) write(LuOut,*) 'HESS: Ta contribution to Hessian'
      do ndns = 1, ndens
c      write(6,*) 'nocc and ncent before hess_sxsy ',nocc(ndns),ncent
       call hess_sxsy2(dbl_mb(k_exy),dbl_mb(k_evals),g_rhs,
     &    log_mb(k_act),ncent,nocc(ndns),hf_job)
      enddo
c
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
      call hess_moootoao(g_rhs, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      if (.not. ga_duplicate (g_rhs, g_tmp,'Temporary g_rhs'))
     $  call errquit('hess_cphf: could not duplicate g_rhs',555)
      call hess_modfock(g_rhs, g_vecs, g_tmp, ndens, n3xyz,
     $                  tol2e, rtdb, geom, basis, jfac, kfac,
     $                  num, ocphfprint)
      if (.not.ga_destroy(g_tmp))
     $  call errquit('hess_cphf: cannot destroy g_tmp',555)
c
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
c
c Now form the Wax and Aa contributions
c
      call hess_fxwaxaa(g_rhs, g_wax, g_aa, nocc, nvirt,
     $                  ndens, n3xyz)
c
c Transform ds/dx term into MO basis and scale by -0.5.
c
      write(LuOut,*) 'Getting Ta again'
      call sx_read(g_rhs,nc,num,sx_dra_handle)
c
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          oprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, num)
      call ga_dscal(g_rhs,-pt5)
      if (ocphfprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        if(ga_nodeid().eq.0)write(LuOut,*) 'Ta after scaling'
        call ga_print(g_rhs)
      endif
c
c Print out intermediate results
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'WAX after fock contribution'
        endif
        call ga_print(g_wax)
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'AA after fock contribution'
        endif
        call ga_print(g_aa)
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'g_rhs after fock contribution'
        endif
        call ga_print(g_rhs)
      endif
c
c Do wxy = wxy + sum(i)sum(j){Ta([i-1]*nocc+j,x)*AA([i-1]*nocc+j,y) +
c Ta([i-1]*nocc+j,y)*AA([i-1]*nocc+j,x)
c
      if (oprint) write(LuOut,*) 'HESS: Ta and AA contributions'
      call hess_taaa(dbl_mb(k_exy),g_rhs, g_aa, nocc, 
     $               ndens, n3xyz)
c
c We can now get rid of the RHS memory
c
      if (.not.ga_destroy(g_rhs))
     $  call errquit('hess_cphf: problem destroying g_rhs',555)
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
      if (.not.ga_destroy(g_aa)) call errquit
     $  ('cannot destroy g_aa ',555)
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
      alo(2) = 1
      alo(3) = 1
      blo(1) = 1
      index = 1
      do ndns = 1, ndens
        ahi(2) = nocc(ndns)
        ahi(3) = nvirt(ndns)
        bhi(1) = nocc(ndns)*nvirt(ndns)
        do ind = 1, nc
          alo(1) = index
          ahi(1) = index
          blo(2) = index
          bhi(2) = index
          call nga_copy_patch('T',g_wax,alo,ahi,
     $      g_tmp,blo,bhi)   ! put into format for current cphf
          index = index + 1
        enddo
      enddo
      if (ocphfprint) then
        if(ga_nodeid().eq.0)write(LuOut,*)
        call ga_print(g_tmp)
      endif
c
      if(.not.file_write_ga(cphf_rhs,g_tmp)) call errquit
     $  ('hess_cphf:could not write cphf_rhs',ind)
c
c Actually run the CPHF on all rhs
c
      if (.not.cphf2(rtdb)) call errquit
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
        call hess_pxfy2(dbl_mb(k_exy),g_sol,g_wax,log_mb(k_act),
     $    ncent,
     $    ind, hf_job, nocc(ndns),nvirt(ndns),ndns)
        index = index + 1
      enddo
      enddo
      if (.not.ga_destroy(g_tmp)) call errquit
     $  ('cannot destroy g_tmp ',555)
      if (.not.ga_destroy(g_tmp2)) call errquit
     $  ('cannot destroy g_tmp2 ',555)
      if (.not.ga_destroy(g_wax)) call errquit
     $  ('cannot destroy g_wax ',555)
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
      integer g_rhs_x, g_vecs(ndens)
      integer g_modx, ind, ndns, in
      integer alo(2), ahi(2), blo(3), bhi(3)
c
c  Transform the current AO based Fock derivatives to the MO basis
c
      if (oprint) write(LuOut,*) 'HESS: setting up Wax and AA'
      if (ocphf) then
        if(ga_nodeid().eq.0) 
     $        write(LuOut,*) 'RHS before MO transformation', ind
        call ga_print(g_rhs_x)
      endif
c
      alo(1) = 1
      ahi(1) = num
      alo(2) = 1
      ahi(2) = num  ! g_vecs and g_modx dimensions
      blo(2) = 1
      bhi(2) = num
      blo(3) = 1
      bhi(3) = num
      ind = 0
      do ndns = 1, ndens
        do in = 1, n3xyz
          ind = ind + 1
          blo(1) = ind
          bhi(1) = ind
          call ga_zero(g_modx)
          call nga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                           g_vecs(ndns), alo, ahi,
     $                           g_rhs_x     , blo, bhi,
     $                           g_modx      , alo, ahi)
          call nga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                           g_modx      , alo, ahi,
     $                           g_vecs(ndns), alo, ahi,
     $                           g_rhs_x     , blo, bhi)
c
        enddo
      enddo
c
      if (ocphf) then
        if(ga_nodeid().eq.0) write(LuOut,*) 
     $        'Transformed RHS after symmetrization',ind
        call ga_print(g_rhs_x)
      endif
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
      integer g_rhs_x, g_vecs(ndens)
      integer g_modx, ind, ndns, in
      integer alo1(2), ahi1(2), blo1(3), bhi1(3)
      integer alo2(2), ahi2(2), blo2(3), bhi2(3)
c
      if(ga_nodeid().eq.0) 
     $    write(LuOut,*) 'RHS before AO transformation'
      call ga_print(g_rhs_x)
c
c  Set up some arrays
c
      alo1(1) = 1
      ahi1(1) = num
      alo1(2) = 1
      alo2(1) = 1
      alo2(2) = 1
      ahi2(2) = num
      blo1(2) = 1
      blo1(3) = 1
      blo2(2) = 1
      bhi2(2) = num
      blo2(3) = 1
      bhi2(3) = num
c
c  Transform the current MO based matrices to the AO basis
c
      ind = 0
      do ndns = 1, ndens
       ahi1(2) = nocc(ndns)
       ahi2(1) = nocc(ndns)
       bhi1(2) = nocc(ndns)
       bhi1(3) = nocc(ndns)
       do in =1,n3xyz
        ind = ind + 1
        blo1(1) = ind
        bhi1(1) = ind
        blo2(1) = ind
        bhi2(1) = ind
        call ga_zero(g_modx)
        call nga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                g_vecs(ndns), alo1, ahi1,
     $                g_rhs_x     , blo1, bhi1,
     $                g_modx      , alo1, ahi1)
        call nga_matmul_patch ('n','t', 1.D0, 0.D0,
     $                g_modx      , alo1, ahi1,
     $                g_vecs(ndns), alo2, ahi2,
     $                g_rhs_x     , blo2, bhi2)
c
       enddo
      enddo
c
      if (ocphf) then
        if(ga_nodeid().eq.0) write(LuOut,*) 
     $        'Transformed RHS '
        call ga_print(g_rhs_x)
      endif
c
      return
      end
c
      subroutine hess_getaawax(g_rhs, g_aa, g_wax, nocc, nvirt,
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
      integer g_rhs, g_aa, g_wax
      integer ind, ndns, in
      integer alo1(3), ahi1(3), alo2(3), ahi2(3)
      integer blo(3), bhi(3), clo(3), chi(3)
c
      write(6,*) 'num and ndens ',num,ndens
      call flush(6)
      ind = 0
      do ndns = 1, ndens
       write(6,*) 'nocc and nvirt ',nocc(ndns),nvirt(ndns)
       alo1(2) = 1
       ahi1(2) = nocc(ndns)
       alo1(3) = nocc(ndns) + 1
       ahi1(3) = num
       blo(2)  = 1
       bhi(2)  = nocc(ndns)
       blo(3)  = 1
       bhi(3)  = nvirt(ndns)
       alo2(2) = 1
       ahi2(2) = nocc(ndns)
       alo2(3) = 1
       ahi2(3) = nocc(ndns)
       clo(2)  = 1
       chi(2)  = nocc(ndns)
       clo(3)  = 1
       chi(3)  = nocc(ndns)
       do in =1,n3xyz
        ind = ind + 1
        alo1(1) = ind
        ahi1(1) = ind
        blo(1)  = ind
        bhi(1)  = ind
        alo2(1) = ind
        ahi2(1) = ind
        clo(1)  = ind
        chi(1)  = ind
        call nga_copy_patch('N',g_rhs,alo1, ahi1,
     $     g_wax,blo, bhi)
        call nga_copy_patch('N',g_rhs,alo2, ahi2,
     $     g_aa, clo, chi)
       enddo
      enddo
c
      if (ocphf) then
        if(ga_nodeid().eq.0) write(LuOut,*) 'g_wax'
        call ga_print(g_wax)
        if(ga_nodeid().eq.0) write(LuOut,*) 'g_aa'
        call ga_print(g_aa)
      endif
c
      return
      end
c
      subroutine hess_tuwax(g_rhs, g_wax, evals, nocc, nvirt,
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
      integer g_rhs, g_wax
      integer ind, ndns, in, i
      integer alo(3), ahi(3), blo(3), bhi(3)
      double precision evals(*)
      double precision one, eight
      data one, eight /1.0d+00, 8.0d+00/
c
      ind = 0
      do ndns = 1, ndens
       alo(3) = nocc(ndns) + 1
       ahi(3) = num
       blo(3) = 1
       bhi(3) = nvirt(ndns)
       do in =1,n3xyz
         ind = ind + 1
         alo(1) = ind
         ahi(1) = ind
         blo(1) = ind
         bhi(1) = ind
         do i = 1, nocc(ndns)
           alo(2) = i
           ahi(2) = i
           blo(2) = i
           bhi(2) = i
           call nga_add_patch(one,g_wax,blo,bhi,
     $       eight*evals((ndns-1)*num+i),g_rhs,alo,ahi,
     $       g_wax,blo,bhi)
         enddo
       enddo
      enddo
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0)write(LuOut,*) 
     $        'g_wax after T addition'
        call ga_print(g_wax)
      endif
c
      return
      end
c
      subroutine hess_taaa(exy,g_rhs, g_aa, nocc, 
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
      integer g_rhs, g_aa
      double precision exy(n3xyz,n3xyz)
      integer ndns, jnd, ind
      integer alo1(3), ahi1(3), alo2(3), ahi2(3)
      integer blo1(3), bhi1(3), blo2(3), bhi2(3)
c
      alo1(2) = 1
      alo1(3) = 1
      alo2(2) = 1
      alo2(3) = 1
      blo1(2) = 1
      blo1(3) = 1
      blo2(2) = 1
      blo2(3) = 1
      do ndns = 1, ndens
      ahi1(2) = nocc(ndns)
      ahi1(3) = nocc(ndns)
      ahi2(2) = nocc(ndns)
      ahi2(3) = nocc(ndns)
      bhi1(2) = nocc(ndns)
      bhi1(3) = nocc(ndns)
      bhi2(2) = nocc(ndns)
      bhi2(3) = nocc(ndns)
      do jnd = 1, n3xyz
        alo2(1) = jnd
        ahi2(1) = jnd
        blo2(1) = jnd
        bhi2(1) = jnd
        do ind = 1, n3xyz
          alo1(1) = ind
          ahi1(1) = ind
          blo1(1) = ind
          bhi1(1) = ind
          exy(ind, jnd) = 
     $      exy(ind, jnd) +
     $      nga_ddot_patch(g_rhs,'N',blo1,bhi1,
     $                    g_aa, 'N',alo2,ahi2) +
     $      nga_ddot_patch(g_rhs,'N',blo2,bhi2,
     $                    g_aa, 'N',alo1,ahi1)
        enddo
      enddo
      enddo
c
      return
      end
c
      subroutine hess_modfock(g_rhs_x, g_vecs, g_tmp, ndens, n3xyz,
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
      integer ind, ndns, in
      integer l_jfac, k_jfac, l_kfac, k_kfac
      double precision tol2e, jfac(ndens), kfac(ndens)
      logical ocphfprint
c
c
c Set up some parameters for the fock_2e call.  Note that I am not
c sure if jfact and kfact should stay the same for DFT!  This needs to
c be dependent on the jfac and kfac in hess_info.fh.
c
c     call ga_maxelt(g_rhs_x,dmax)
c     if (dmax.gt.tol2e) then  ! the value of tol2e may not be the best here
      if (.not. ma_push_get(MT_DBL, ndens*n3xyz,'j factors',
     $     l_jfac, k_jfac))
     $     call errquit
     $     ('hess_modfock: could not allocate l_jfac',ndens*n3xyz)
      if (.not. ma_push_get(MT_DBL, ndens*n3xyz,'k factors',
     $     l_kfac, k_kfac))
     $     call errquit
     $     ('hess_modfock: could not allocate l_kfac',ndens*n3xyz)
      ind = 0
      do ndns = 1, ndens
        do in = 1, n3xyz
          dbl_mb(k_jfac+ind) =  8.0d+00*jfac(ndns)
          dbl_mb(k_kfac+ind) =  8.0d+00*kfac(ndns)
          ind = ind + 1
        enddo
      enddo

      call ga_copy(g_rhs_x,g_tmp)  ! copy Ta to g_tmp
      call ga_zero(g_rhs_x)   ! zeroed to get the "new" fock matrices

      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
c         write(LuOut,*) 'nfock ',nfock
          write(LuOut,*) 'density like matrix before fock_2e'
        endif
        call ga_print(g_tmp)
      endif
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call scf_get_fock_param(rtdb,tol2e) ! need to change for DFT?
c
c     Activate XC contribution to OO derivative Fock matrix
c
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of j_derfit failed',0)
c
c     call shell_fock_build(basis,0,ndens,jfac,kfac,tol2e,
c    &            g_tmp,g_rhs_x)
      call shell_fock_build(basis,0,ndens*n3xyz,dbl_mb(k_jfac),
     &            dbl_mb(k_kfac),tol2e,
     &            g_tmp,g_rhs_x,.true.)
c
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',1)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0)
      call schwarz_tidy()
      call int_terminate()
c
      if (.not.ma_pop_stack(l_kfac))
     $  call errquit('hess_cphf: pop problem with l_jfac',555)
      if (.not.ma_pop_stack(l_jfac))
     $  call errquit('hess_cphf: pop problem with l_kfac',555)
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'fock matrix after fock_2e'
        endif
        call ga_print(g_rhs_x)
      endif

c     ind = 0
c     do ndns = 1, ndens
c     do in = 1, n3xyz
c       ind = ind+1
c       call ga_zero(g_tmp) ! prepare for transformation
c       write(LuOut,*) 'before first ga_matmul_patch'
c       call ga_matmul_patch ('t','n', 1.D0, 0.D0,
c    $                         g_vecs(1),    1, num,    1, num,
c    $                         g_rhs_x,   1, num,    1, num,
c    $                         g_tmp,       1, num,    1, num )
c
c       write(LuOut,*) 'before second ga_matmul_patch'
c       call ga_matmul_patch ('n','n', 1.D0, 0.D0,
c    $                         g_tmp,       1, num,    1, num,
c    $                         g_vecs(1),    1, num,    1, num,
c    $                         g_rhs_x,   1, num,    1, num ) 
c     enddo
c     enddo

c
c Put jfac and kfac back to what they were before entering the routine
c
c     do ndns = 1, ndens
c       jfac(ndns) = jfac(ndns)/8.0d+00
c       kfac(ndns) = kfac(ndns)/-4.0d+00
c     enddo
c
c     if (ocphfprint) then
c       if(ga_nodeid().eq.0)
c    $    write(LuOut,*) 'after transformation back to MO basis',ind
c       call ga_print(g_rhs_x)
c     endif
c
      return
      end
c
      subroutine hess_fxwaxaa(g_rhs, g_wax, g_aa, nocc, nvirt,
     $                        ndens, n3xyz)
c
      implicit none
c
#include "global.fh"
c
c Now form the Wax and Aa contributions
c
      integer ndens, n3xyz
      integer g_rhs, g_wax, g_aa
      integer nocc(ndens), nvirt(ndens)
      integer ind, ndns, in
      integer alo1(3), ahi1(3), blo1(3), bhi1(3)
      integer alo2(3), ahi2(3), blo2(3), bhi2(3)
      double precision one, two
      data one,two         /1.0d+00, 2.0d+00/
c
c  Set up some arrays
c
      alo1(2) = 1
      alo1(3) = 1
      alo2(2) = 1
      alo2(3) = 1
      blo1(2) = 1
      blo2(2) = 1
      blo2(3) = 1
      ind = 0
      do ndns = 1, ndens
       ahi1(2) = nocc(ndns)
       ahi1(3) = nvirt(ndns)
       ahi2(2) = nocc(ndns)
       ahi2(3) = nocc(ndns)
       blo1(3) = nocc(ndns) + 1
       bhi1(2) = nocc(ndns)
       bhi1(3) = nocc(ndns) + nvirt(ndns)
       bhi2(2) = nocc(ndns)
       bhi2(3) = nocc(ndns)
       do in = 1, n3xyz
        ind = ind+1
        alo1(1) = ind
        ahi1(1) = ind
        blo1(1) = ind
        bhi1(1) = ind 
        blo2(1) = ind
        bhi2(1) = ind
        alo2(1) = ind
        ahi2(1) = ind
        call nga_add_patch(one,g_wax,alo1,ahi1,
     $                     two,g_rhs,blo1,bhi1,
     $                         g_wax,alo1,ahi1)

        call nga_add_patch(one,g_aa,alo2,ahi2,
     $                     one,g_rhs,blo2,bhi2,
     $                         g_aa,alo2,ahi2)
       enddo
      enddo
c
      return
      end
