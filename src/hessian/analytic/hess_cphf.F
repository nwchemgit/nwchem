C
C $Id$
C
      subroutine hess_cphf(rtdb,restr)
c
c  Solves the CPHF equations to get the responses
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "util.fh"
#include "geom.fh"
c
      integer rtdb 
      integer ind, index, index1
      integer nc
      integer nbf2, ndns, ndip
      integer l_evals, k_evals, l_occ, k_occ
      integer g_vecs(2), g_modx, g_dip
      integer g_wax
      integer g_aa
      integer g_tmp, g_sol
c
      integer ilo(3), ihi(3)
      integer alo(3), ahi(3), blo(2), bhi(2)
      integer dims(3), chunk(3)
      integer i, idens
c
      character*256 cphf_rhs, cphf_sol
      character*255 movecs
      character*255 tmpname
      logical hf_job, oprint, olprint, ocphfprint, ouhf
      double precision tol2e
      double precision pt5, four, zero
      double precision dnrm, enuc, enuc_in
      data pt5, four, zero/0.5d+00, 4.0d+00, 0.0d+00/
c
      integer  restr
      logical  cphf2, file_write_ga, file_read_ga, movecs_read
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      external cphf2, file_write_ga, file_read_ga, movecs_read
c
      hf_job = theory .eq. 'hf'
      ouhf   = scftype .eq. 'UHF'
c
c Set up printing
c
      oprint  = util_print('hess_follow',print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
      ocphfprint = util_print('cphf_cont',print_debug)
      if (olprint) write(LuOut,*) 'HESS: setting up cphf'
c 
c     ----- define parameters -----
c
      nc   = n3xyz
      call cphf_fname('cphf_rhs',cphf_rhs)
      call cphf_fname('cphf_sol',cphf_sol)
c
c Get MO vectors, first figure out name of movecs file, get memory, and
c then get actual data
c
      if (theory.eq.'dft') then
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs))
     *    call errquit('onedd_cont: DFT MO vectors not defined',555,
     &       RTDB_ERR)
      else
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs))
     $    call errquit('gradients: SCF MO vectors not defined',0,
     &       RTDB_ERR)
      endif
      call util_file_name_resolve(movecs,.false.)
c
c  Get memory, note that the occ space is reused for
c  UHF since it isn't used in the actual calculation
c
      nbf2 = nbf
      if (ouhf) nbf2 = nbf + nbf
      g_vecs(1) = ga_create_atom_blocked(geom, basis, 'mo vectors')
c     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('hess_cphf: could not allocate l_evals',nbf2,
     &       MA_ERR)
c     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf,'occ. numbers', l_occ, k_occ))
     $     call errquit('hess_cphf: could not allocate l_occ',nbf,
     &       MA_ERR)
c
      call ga_zero(g_vecs(1))
      call dfill(nbf,zero,dbl_mb(k_occ),1)
      call dfill(nbf,zero,dbl_mb(k_evals),1)
c
      if (.not. movecs_read (movecs,1,dbl_mb(k_occ),dbl_mb(k_evals), 
     $     g_vecs(1)))
     $     call errquit('hess_cphf: could not read mo vectors', 110,
     &       DISK_ERR)

      if (ouhf) then
         g_vecs(2) = ga_create_atom_blocked(geom, basis, 
     $       'uhf mo vectors')
         call ga_zero(g_vecs(2))
         call dfill(nbf,zero,dbl_mb(k_evals+nbf),1)
         if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), 
     $        dbl_mb(k_evals+nbf), g_vecs(2))) call errquit
     $        ('hess_cphf: could not read beta UHF MOs',0,
     &       DISK_ERR)
      end if
c
c print out some debug info, note that occ is for the
c beta orbitals if UHF and the evals are for the alpha orbitals
c
      if (ocphfprint) then
         if(ga_nodeid().eq.0) then
           write(LuOut,*)
           if (ouhf) write(LuOut,*) 'Please note that for the rest ',
     $       'of the CPHF that g_wax and g_aa have the maximum ',
     $       'alpha or beta dimension'
           write(LuOut,*) 'occupations, eigenvalues and vectors'
           write(LuOut,*) 'occupations are for the beta values'
           call hnd_prsq(dbl_mb(k_occ),1,nmo,nbf)
           write(LuOut,*)
           call hnd_prsq(dbl_mb(k_evals),1,nmo,nbf)
           write(LuOut,*)
           if (ouhf) then
             write(LuOut,*) 'Beta eigenvalues'
             call hnd_prsq(dbl_mb(k_evals+nbf),1,nmo,nbf)
             write(LuOut,*)
           endif
           call util_flush(LuOut)
         endif
         call ga_print(g_vecs(1))
         if (ouhf) call ga_print(g_vecs(2))
         if(ga_nodeid().eq.0) write(LuOut,*)
      endif
      if (oprint) then
        call ga_normf(g_vecs(1),dnrm)
        if (olprint) then
          write(LuOut,*)'hess_cphf',': norm of vectors = ',dnrm
        endif
        if (ouhf) then
          call ga_normf(g_vecs(2),dnrm)
          if (olprint) then
            write(LuOut,*)'hess_cphf',': norm of beta vectors = ',dnrm
          endif
        endif
      endif
c
c give back unnecessary memory
c
      if (.not.ma_pop_stack(l_occ)) call errquit
     &  ('hess_cphf: cannot pop stack l_occ',555, MA_ERR)
c
      g_modx =ga_create_atom_blocked ( geom, basis, 'mo df' )
c
c Note that the maximum size of the occupied and virtual spaces are used so that
c we need only one multidimensioned array.
c
      dims(1) = ndens*n3xyz
      dims(2) = nocc(1)
      if (ouhf) dims(2) = max(nocc(1), nocc(2))
      dims(3) = dims(2)
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'AA matrix',chunk,
     &     g_aa))
     &   call errquit('hess_cphf: could not allocate g_aa',555, GA_ERR)
c
      dims(3) = nvirt(1)
      if (ouhf) dims(3) = max(nvirt(1),nvirt(2))
      if (.not. nga_create (MT_DBL, 3, dims, 'WAX matrix',chunk,
     &     g_wax))
     &   call errquit('hess_cphf: could not allocate g_wax',555, GA_ERR)
c
c  Transform the current AO based Fock derivatives to the MO basis,
c  scale by four, and populate g_aa and g_wax 
c
      call hess_copymat(ndens,n3xyz,nbf,nbf,g_rhs_fx,g_rhs)
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          olprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, nbf)
      call ga_dscal(g_rhs,four)
      call ga_zero(g_wax)
      call ga_zero(g_aa)
      call hess_getaawax(g_rhs, g_aa, g_wax, nocc, nvirt,
     $          ndens, n3xyz, nmo, ocphfprint)
c
c We are now done with the Fock derivatives.
c Will reuse g_rhs for the ds/dx term.
c Transform ds/dx term into MO basis and scale by -0.5.
c
      call hess_copymat(ndens,n3xyz,nbf,nbf,g_rhs_sx,g_rhs)
c
      if (olprint) write(LuOut,*) 'HESS: Ta contribution to Wax'
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          olprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, nbf)
      call ga_dscal(g_rhs,-pt5)
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'Ta after scaling'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_rhs,ilo,ihi,1)
          enddo
        enddo
      endif
      if (oprint) then
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_normf_patch(g_rhs,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": Ta =",i2,i4,f24.8)')
     +        'hess_cphf',
     +        idens,i,dnrm
            endif
          enddo
        enddo
      endif
c
c get contribution to wax
c
      call hess_tuwax(g_rhs, g_wax, dbl_mb(k_evals), nocc, nvirt,
     $          ocphfprint, ndens, n3xyz, nbf)
c
c Do wxy = wxy + 8*sum(i)sum(j)[(ei+ej)*Ta((i-1*nocc)+j,x))*Ta((i-1*nocc)+j,y)]
c This is the hondo equivalent of cpwab1 [ Ta is part of g_rhs, e is 
c dbl_mb(k_evals) and wxy is dbl_mb(k_exy)
c
      if (olprint) then
        write(LuOut,*) 'HESS: Ta contribution to Hessian: before'
        call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
      do ndns = 1, ndens
        if (nocc(ndns).gt.0) then
           call hess_sxsy(dbl_mb(k_exy),
     &        dbl_mb(k_evals+(ndns-1)*nbf),g_rhs,
     &        log_mb(k_act),ncent,nocc(ndns),(ndns-1)*n3xyz,ouhf)
        endif
      enddo
c
c
c Do final summation of the hessian to get result
c
      call ga_dgop(msg_hess_exy,dbl_mb(k_exy),n3xyz2,'+') 
c
c Print intermediate result
c
      if (olprint.or.(ocphfprint.and.ga_nodeid().eq.0)) then
        write(LuOut,*)
        write(LuOut,*) 'Hessian after Ta*Tb contribution'
        call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
c
c  Get rid of some memory
c
      if (.not.ma_pop_stack(l_evals)) call errquit
     &  ('hess_cphf: cannot pop stack l_evals',555, MA_ERR)
c
c Do modified fock build that contracts integrals with Ta and then
c adds to g_aa and g_wax.  Ta needs to be transformed back into the AO
c basis, combined with the integrals, and then transformed back into the
c MO basis.
c
      if (olprint) write(LuOut,*) 'HESS: Modified Fock builds'
      if (.not. rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)) 
     $  tol2e = 1.0d-9
c
      call hess_moootoao(g_rhs, g_vecs, g_modx,
     $          olprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, nbf)
      if (.not. ga_duplicate (g_rhs, g_tmp,'Temporary g_rhs')) then
        if (ga_nodeid().eq.0) then
          write(LuOut,*) 'Could not duplicate g_rhs'
          write(LuOut,*) 'The most likely reason is not enough memory'
        endif
        call errquit('hess_cphf: could not duplicate g_rhs',555, GA_ERR)
      endif
      call hess_modfock(g_rhs, g_vecs, g_tmp, ndens, n3xyz,
     $                  tol2e, rtdb, geom, basis, jfac, kfac,
     $                  nbf, ocphfprint)
      if (.not.ga_destroy(g_tmp))
     $  call errquit('hess_cphf: cannot destroy g_tmp',555, GA_ERR)
c
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          olprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, nbf)
c
c Now form the Wax and Aa contributions
c
      call hess_fxwaxaa(g_rhs, g_wax, g_aa, nocc, nvirt,
     $                  ndens, n3xyz)
c
c Transform ds/dx term into MO basis and scale by -0.5.
c
      if (ocphfprint) then
        write(LuOut,*) 'Getting Ta again'
        call util_flush(LuOut)
      endif
      call hess_copymat(ndens,n3xyz,nbf,nbf,g_rhs_sx,g_rhs)
c
      if (.not.ga_destroy(g_rhs_sx))
     $  call errquit('hess_cphf: problem destroying g_rhs',555, GA_ERR)

      if (.not.ga_destroy(g_rhs_fx))
     $  call errquit('hess_cphf: problem destroying g_rhs',555, GA_ERR)
c
      call hess_aotomo(g_rhs, g_vecs, g_modx,
     $          olprint, ocphfprint,
     $          nocc, nvirt, ndens, n3xyz, nbf)
c
c Before doing the scaling, we will get the dipole integrals into 3
c (or 6 for UHF) matrices (x, y, z), transform them into the MO 
c basis and then compute the overlap derivative times the dipole 
c integral contribution to the dipole derivative term.
c
      dims(1) = ndens*3
      dims(2) = nbf
      dims(3) = dims(2)
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'dipole matrix',chunk,
     &     g_dip))
     &   call errquit('hess_cphf: could not allocate g_dip',555, GA_ERR)
      call ga_zero(g_dip)
      call hess_dipint(rtdb, geom, basis, g_dip, ndens, nbf)
      ndip = 3
      call hess_aotomo(g_dip, g_vecs, g_modx, olprint,
     $          ocphfprint, nocc, nvirt, ndens, ndip, nbf)
      call hess_sxdip(g_rhs, dbl_mb(k_dder), g_dip, nocc, ndens, 
     $   n3xyz, nbf)
      if (ocphfprint.and.ga_nodeid().eq.0) then
        write(luout,*) 'after dipole times overlap derivative'
        write(tmpname,*) ' derivative dipole '
        call stpr_print_ddipole(dbl_mb(k_dder),tmpname,ncent,-1.0d00)
      endif
c
c Now scale
c
      call ga_dscal(g_rhs,-pt5)
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'Ta after scaling'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_rhs,ilo,ihi,1)
          enddo
        enddo
      endif
      if (oprint) then
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_normf_patch(g_rhs,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": Ta after scaling =",i2,i4,f24.8)')
     &        'hess_cphf',
     &        idens,i,dnrm
            endif
          enddo
        enddo
      endif
c
c Print out intermediate results
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'WAX after fock contribution'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)
            ihi(3) = nvirt(idens)
            call nga_print_patch(g_wax,ilo,ihi,1)
          enddo
        enddo
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'AA after fock contribution'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)
            ihi(3) = nocc(idens)
            call nga_print_patch(g_aa,ilo,ihi,1)
          enddo
        enddo
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'g_rhs after fock contribution'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)
            ihi(3) = nvirt(idens)
            call nga_print_patch(g_wax,ilo,ihi,1)
          enddo
        enddo
      endif
      if (oprint) then
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)
            ihi(3) = nvirt(idens)
            call nga_normf_patch(g_wax,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": WAX after fock =",i2,i4,f24.8)')
     &        'hess_cphf',
     &        idens,i,dnrm
            endif
          enddo
        enddo
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)
            ihi(3) = nocc(idens)
            call nga_normf_patch(g_aa,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": AA after fock =",i2,i4,f24.8)')
     &        'hess_cphf',
     &        idens,i,dnrm
            endif
          enddo
        enddo
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_normf_patch(g_rhs,ilo,ihi,dnrm)
            if (olprint) then
              write(LuOut,*)'hess_cphf',": RHS after fock = ",
     &                      idens,i,dnrm
            endif
          enddo
        enddo
      endif
c
c Do wxy = wxy + sum(i)sum(j){Ta([i-1]*nocc+j,x)*AA([i-1]*nocc+j,y) +
c Ta([i-1]*nocc+j,y)*AA([i-1]*nocc+j,x)
c
      if (olprint) then
         write(LuOut,*) 'HESS: Ta and AA contributions'
         call util_flush(LuOut)
      endif
      if (ndens.gt.1) call ga_dscal(g_aa,0.5d0)
      call hess_taaa(dbl_mb(k_exy),g_rhs, g_aa, nocc, 
     $               ndens, n3xyz)
c
c We can now get rid of the RHS memory
c
      if (.not.ga_destroy(g_rhs))
     $  call errquit('hess_cphf: problem destroying g_rhs',555, GA_ERR)

c      if (.not.ga_destroy(g_rhs_sx))
c     $  call errquit('hess_cphf: problem destroying g_rhs',555, GA_ERR)
c
c      if (.not.ga_destroy(g_rhs_fx))
c     $  call errquit('hess_cphf: problem destroying g_rhs',555, GA_ERR)

c
c Print intermediate result
c
      if (olprint.or.(ocphfprint.and.(ga_nodeid().eq.0))) then
        write(LuOut,*)
        write(LuOut,*) 'Hessian after Ta*AA contribution'
        call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
c
c get rid of some extra memory
c
      if (.not.ga_destroy(g_modx))
     $  call errquit('hess_cphf: cannot destroy g_modf',555, GA_ERR)
      if (.not.ga_destroy(g_aa)) call errquit
     $  ('cannot destroy g_aa ',555, GA_ERR)
c
c Now actually set up the CPHF matrices 
c
      if (olprint) write(LuOut,*) 'HESS: CPHF contributions'
c
c quick and dirty scale of g_wax for UHF
c
       if (ndens.gt.1) call ga_dscal(g_wax,0.5d0)
c
c Create a temporary array to reshape the CPHF RHS and send to the CPHF
c
       if (.not.ga_create(MT_DBL,vlen,nc,'cphf_sol',
     $    0,0,g_sol))
     $    call errquit('hess_cphf: could not create g_sol',555, GA_ERR)
c
c Try reading RHS if it is available
c
      if (restr.gt.0) then
      if(.not.geom_nuc_rep_energy(geom, enuc))
     c     call errquit('hess_cphf: cannot get enuc',0,0)
      if(.not. rtdb_get(rtdb,'hess:cphf_enuc', MT_DBL, 1, enuc_in))
     t     then
         if(ga_nodeid().eq.0) write(luout,*)
     A        ' cannot reuse RHS because of no enuc in rtdb '
         restr=0
      else
         if(abs(enuc-enuc_in).gt.1d-3) then
            if(ga_nodeid().eq.0) write(luout,*)
     A           ' cannot reuse RHS because of enuc mismatch ',
     E           abs(enuc-enuc_in)
            restr=0
         endif
      endif
      if(restr.eq.0) then
      if (.not.rtdb_put(rtdb,'solver:restart',mt_int,1,restr))
     &        call errquit('hess_cphf: cannot store restr ',0,RTDB_ERR)
         goto 50
      endif
       if(.not.file_read_ga(cphf_rhs,g_sol)) goto 50
       if (ga_nodeid().eq.0) 
     &      write(luout,*) "Restarting rhs from: ", cphf_rhs
       goto 100  ! jump to cphf call
      end if ! restr 
c
 50   continue
c
c Put the CPHF matrices into one GA in the correct format
c Note that the alpha and beta matrices are concatenated
c (not one after the other)
c
      alo(2) = 1
      alo(3) = 1
      bhi(1) = 0
      index1 = 1
      blo(1) = 1
      do ndns = 1, ndens
        ahi(2) = nocc(ndns)
        ahi(3) = nvirt(ndns)
        blo(1) = bhi(1)+ 1
        bhi(1) = bhi(1) + nocc(ndns)*nvirt(ndns)
        do ind = 1, nc
          alo(1) = index1
          ahi(1) = index1
          blo(2) = ind
          bhi(2) = ind
          call ga_mask_sync(.true.,.false.)
          call nga_copy_patch('T',g_wax,alo,ahi,
     $      g_sol,blo,bhi)   ! put into format for current cphf
          index1 = index1 + 1
        enddo
      enddo
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0)then
          write(LuOut,*)
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            blo(1) = ind
            blo(2) = 1
            bhi(1) = ind
            bhi(2) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_sol,blo,bhi,1)
          enddo
        enddo
      endif
      if (oprint) then
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            blo(1) = ind
            blo(2) = 1
            bhi(1) = ind
            bhi(2) = nocc(idens)+nvirt(idens)
            call nga_normf_patch(g_sol,blo,bhi,dnrm)
            if (olprint) then
              write(LuOut,'(a,": solution =",i2,i4,f24.8)')
     &        'hess_cphf',
     &        idens,i,dnrm
            endif
          enddo
        enddo
      endif
c
 100  continue
c
c store nucl. rep. energy to verify that restart rhs corresponds to the same geom
c
      if(.not.geom_nuc_rep_energy(geom, enuc))
     c     call errquit('hess_cphf: cannot get enuc',0,0)
      if(.not. rtdb_put(rtdb,'hess:cphf_enuc', MT_DBL, 1, enuc))
     $     call errquit('hess_cphf: rtdb_put of enuc failed',1,
     &       RTDB_ERR)
c
      if(.not.file_write_ga(cphf_rhs,g_sol)) call errquit
     $  ('hess_cphf:could not write cphf_rhs',ind, DISK_ERR)
      if (.not.ga_destroy(g_sol)) call errquit
     $  ('cannot destroy g_sol ',555, GA_ERR)
      if (.not.ga_destroy(g_wax)) call errquit
     $  ('cannot destroy g_wax ',555, GA_ERR)
c
c Actually run the CPHF on all rhs
c
      if (.not.cphf2(rtdb)) call errquit
     $  ('hess_cphf: failure in cphf ',555, RTDB_ERR)
c
c Get the RHS and solution off of disk.
c
      if (.not.ga_create(MT_DBL,vlen,nc,'cphf_rhs',
     $   0,0,g_wax))
     $   call errquit('hess_cphf: could not create g_wax',555, GA_ERR)
      if(.not.file_read_ga(cphf_rhs,g_wax)) call errquit
     $  ('hess_cphf:could not read cphf_wax',555, DISK_ERR)
      if (.not.ga_create(MT_DBL,vlen,nc,'cphf_sol',
     $   0,0,g_sol))
     $   call errquit('hess_cphf: could not create g_sol',555, GA_ERR)
      if(.not.file_read_ga(cphf_sol,g_sol)) call errquit
     $  ('hess_cphf:could not read cphf_rhs',555, DISK_ERR)
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'g_wax after disk'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_wax,ilo,ihi,1)
          enddo
        enddo
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'g_sol after disk'
          call util_flush(LuOut)
        endif
        ind = 0
        do idens = 1, ndens
          do i = 1, n3xyz
            ind = ind + 1
            ilo(1) = ind
            ilo(2) = 1
            ilo(3) = 1
            ihi(1) = ind
            ihi(2) = nocc(idens)+nvirt(idens)
            ihi(3) = nocc(idens)+nvirt(idens)
            call nga_print_patch(g_sol,ilo,ihi,1)
          enddo
        enddo
      endif
c
c Combine the RHS and CPHF solution to get final contribution to hessian.
c
      call hess_pxfy(dbl_mb(k_exy),g_sol,g_wax,log_mb(k_act),
     $  ncent, hf_job, vlen)
c
c Get rid of temp memory
c
      if (.not.ga_destroy(g_wax)) call errquit
     $  ('cannot destroy g_wax ',555, GA_ERR)
c
c Now get the dipole derivative term
c
      call hess_resdip(g_sol,dbl_mb(k_dder),g_dip,nocc,nvirt,
     $  ndens, n3xyz, nmo)
      if (ocphfprint.and.ga_nodeid().eq.0) then
        write(tmpname,*) ' derivative dipole '
        call stpr_print_ddipole(dbl_mb(k_dder),tmpname,ncent,-1.0d00)
      endif
c
c Get rid of temp memory
c
      if (.not.ga_destroy(g_sol)) call errquit
     $  ('cannot destroy g_sol ',555, GA_ERR)
      if (.not.ga_destroy(g_dip)) call errquit
     $  ('cannot destroy g_dip ',555, GA_ERR)
c
c Print intermediate result
c
      if(olprint.or.(ocphfprint.and.ga_nodeid().eq.0)) then
         write(LuOut,*)
         write(LuOut,*) 'Hessian after CPHF contribution'
         call hess_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
      endif
c
c Get rid of some extra memory
c
      do ind = 1, ndens
        if (.not.ga_destroy(g_vecs(ind))) call errquit
     $    ('cannot destroy g_vecs ',ind, GA_ERR)
      enddo
c
      return
      end
c
      subroutine hess_aotomo(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphf,
     $          nocc, nvirt, ndens, n3xyz, nbf)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical oprint, ocphf
      integer ndens, n3xyz, nbf, nocc(ndens), nvirt(ndens)
      integer g_rhs_x, g_vecs(ndens)
      integer g_modx, ind, ndns, in, nmo
      integer alo(2), ahi(2), blo(3), bhi(3)
      integer clo(3), chi(3), dlo(2), dhi(2)
c
c  Transform the current AO based Fock derivatives to the MO basis
c
      if (ocphf) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'RHS before MO transformation'
          call util_flush(LuOut)
        endif
        call ga_print(g_rhs_x)
      endif
c
      alo(1) = 1
      alo(2) = 1
      ahi(1) = nbf
      dlo(1) = 1
      dlo(2) = 1
      dhi(2) = nbf
      blo(2) = 1
      blo(3) = 1
      bhi(2) = nbf
      bhi(3) = nbf
      clo(2) = 1
      clo(3) = 1
      ind = 0
      do ndns = 1, ndens
        nmo = nocc(ndns) + nvirt(ndns)
        ahi(2) = nmo  ! g_vecs and g_modx dimensions
        chi(2) = nmo
        chi(3) = nmo 
        dhi(1) = nmo
        do in = 1, n3xyz
          ind = ind + 1
          blo(1) = ind
          bhi(1) = ind
          clo(1) = ind
          chi(1) = ind
          call ga_mask_sync(.true.,.false.)
          call ga_zero(g_modx)
          call ga_mask_sync(.true.,.false.)
          call nga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                           g_vecs(ndns), dlo, dhi,
     $                           g_rhs_x     , blo, bhi,
     $                           g_modx      , dlo, dhi)
          call ga_mask_sync(.true.,.false.)
          call nga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                           g_modx      , dlo, dhi,
     $                           g_vecs(ndns), alo, ahi,
     $                           g_rhs_x     , clo, chi)
c
        enddo
      enddo
      call ga_sync()
c
      if (ocphf) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)'Transformed RHS in MO basis'
          call util_flush(LuOut)
        endif
        call ga_print(g_rhs_x)
      endif
c
      return
      end
c
      subroutine hess_moootoao(g_rhs_x, g_vecs, g_modx,
     $          oprint, ocphf,
     $          nocc, nvirt, ndens, n3xyz, nbf)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical oprint, ocphf
      integer ndens, n3xyz, nbf, nocc(ndens), nvirt(ndens)
      integer g_rhs_x, g_vecs(ndens)
      integer g_modx, ind, ndns, in
      integer alo1(2), ahi1(2), blo1(3), bhi1(3)
      integer alo2(2), ahi2(2), blo2(3), bhi2(3)
c
      if (ocphf) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'RHS before AO transformation'
          call util_flush(LuOut)
        endif
        call ga_print(g_rhs_x)
      endif
c
c  Set up some arrays
c
      alo1(1) = 1
      ahi1(1) = nbf
      alo1(2) = 1
      alo2(1) = 1
      alo2(2) = 1
      ahi2(2) = nbf
      blo1(2) = 1
      blo1(3) = 1
      blo2(2) = 1
      bhi2(2) = nbf
      blo2(3) = 1
      bhi2(3) = nbf
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
        if (nocc(ndns).le.0) then
          call nga_zero_patch(g_rhs_x,blo2,bhi2)
        else ! nocc(ndns).le.0
           call ga_mask_sync(.true.,.false.)
        call ga_zero(g_modx)
        
           call ga_mask_sync(.true.,.false.)
        call nga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                g_vecs(ndns), alo1, ahi1,
     $                g_rhs_x     , blo1, bhi1,
     $                g_modx      , alo1, ahi1)
           call ga_mask_sync(.true.,.false.)
        call nga_matmul_patch ('n','t', 1.D0, 0.D0,
     $                g_modx      , alo1, ahi1,
     $                g_vecs(ndns), alo2, ahi2,
     $                g_rhs_x     , blo2, bhi2)

        endif ! nocc(ndns).le.0
c
       enddo
      enddo
      call ga_sync()
c
      if (ocphf) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)'Transformed RHS in AO basis'
          call util_flush(LuOut)
        endif
        call ga_print(g_rhs_x)
      endif
c
      return
      end
c
      subroutine hess_getaawax(g_rhs, g_aa, g_wax, nocc, nvirt,
     $          ndens, n3xyz, nmo, ocphf)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical ocphf
      integer ndens, n3xyz, nocc(ndens), nvirt(ndens), nmo
      integer g_rhs, g_aa, g_wax
      integer ind, ndns, in
      integer alo1(3), ahi1(3), alo2(3), ahi2(3)
      integer blo(3), bhi(3), clo(3), chi(3)
c
      ind = 0
      do ndns = 1, ndens
       alo1(2) = 1
       ahi1(2) = nocc(ndns)
       alo1(3) = nocc(ndns) + 1
       ahi1(3) = nmo
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
          call ga_mask_sync(.true.,.false.)
        call nga_copy_patch('N',g_rhs,alo1, ahi1,
     $     g_wax,blo, bhi)
          call ga_mask_sync(.true.,.false.)
        call nga_copy_patch('N',g_rhs,alo2, ahi2,
     $     g_aa, clo, chi)
       enddo
      enddo
      call ga_sync()
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
     $          ocphfprint, ndens, n3xyz, nbf)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical ocphfprint
      integer ndens, n3xyz, nocc(ndens), nvirt(ndens), nbf
      integer g_rhs, g_wax
      integer ndns,  i
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer ind0,ind1
      double precision evals(*)
c
      blo(3) = 1
      do ndns = 1, ndens
         alo(3) = nocc(ndns) + 1
         ahi(3) = nocc(ndns) + nvirt(ndns)
         bhi(3) = nvirt(ndns)
         ind0 = (ndns-1)*n3xyz+1
         ind1 = ndns*n3xyz
         alo(1) = ind0
         ahi(1) = ind1
         blo(1) = ind0
         bhi(1) = ind1
         do i = 1, nocc(ndns)
            alo(2) = i
            ahi(2) = i
            blo(2) = i
            bhi(2) = i
            call ga_mask_sync(.true.,.false.)
            call nga_add_patch(1d0,g_wax,blo,bhi,
     $           8d0*evals((ndns-1)*nbf+i),g_rhs,alo,ahi,
     $           g_wax,blo,bhi)
         enddo
      enddo
      call ga_sync()
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)'g_wax after T addition'
          call util_flush(LuOut)
        endif
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
      integer ndns, jnd, ind, jindex, index
      integer alo1(3), ahi1(3), alo2(3), ahi2(3)
      integer blo1(3), bhi1(3), blo2(3), bhi2(3)
      double precision vtrace
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
      jindex = (ndns-1)*n3xyz
      do jnd = 1, n3xyz
        jindex = jindex + 1
        alo2(1) = jindex
        ahi2(1) = jindex
        index = (ndns-1)*n3xyz
c        do ind = 1, n3xyz
        do ind = 1, jnd
          index = index + 1
          alo1(1) = index
          ahi1(1) = index
          call ga_mask_sync(.true.,.false.)
          vtrace=
     $      nga_ddot_patch(g_rhs,'N',alo1,ahi1,
     $         g_aa, 'N',alo2,ahi2) 
          if(ind.ne.jnd) then
          call ga_mask_sync(.true.,.false.)
             vtrace=vtrace+
     $      nga_ddot_patch(g_rhs,'N',alo2,ahi2,
     $            g_aa, 'N',alo1,ahi1)
          else
             vtrace=vtrace*2d0
          endif
          exy(ind, jnd) = 
     $      exy(ind, jnd) + vtrace
          if(ind.ne.jnd) exy(jnd, ind) = 
     $      exy(jnd, ind) + vtrace

        enddo
      enddo
      enddo
      call ga_sync()
c
      return
      end
c
      subroutine hess_modfock(g_rhs_x, g_vecs, g_tmp, ndens, n3xyz,
     $                        tol2e, rtdb, geom, basis, jfac, kfac,
     $                        nbf, ocphfprint)
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer ndens, n3xyz, rtdb, geom, basis, nbf
      integer g_rhs_x, g_vecs(ndens), g_tmp
      integer ind, ndns, in, ndensity, nfock
      integer l_jfac, k_jfac, l_kfac, k_kfac, nfacs
      double precision tol2e, jfac(ndens), kfac(ndens)
      double precision multfac
      logical ocphfprint
c
c
c Set up some parameters for the fock_2e call.  Note that I am not
c sure if jfact and kfact should stay the same for DFT!  This needs to
c be dependent on the jfac and kfac in hess_info.fh.
c
c     call ga_maxelt(g_rhs_x,dmax)
c     if (dmax.gt.tol2e) then  ! the value of tol2e may not be the best here
      nfacs = ndens*n3xyz
      if (ndens.eq.2) nfacs = nfacs*2
      if (.not. ma_push_get(MT_DBL, nfacs,'j factors',
     $     l_jfac, k_jfac))
     $     call errquit
     $     ('hess_modfock: could not allocate l_jfac',ndens*n3xyz,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, nfacs,'k factors',
     $     l_kfac, k_kfac))
     $     call errquit
     $     ('hess_modfock: could not allocate l_kfac',ndens*n3xyz,
     &       MA_ERR)
c
      multfac = 8.0d+00
      if (ndens.gt.1) multfac = 4.0d+00
c
      ind = 0
      do ndns = 1, ndens
        do in = 1, n3xyz
          dbl_mb(k_jfac+ind) =  jfac(ndns)
c         dbl_mb(k_kfac+ind) =  0.0d+00*kfac(ndns)
          dbl_mb(k_kfac+ind) =  kfac(ndns)
          ind = ind + 1
        enddo
      enddo

      call ga_copy(g_rhs_x,g_tmp)  ! copy Ta to g_tmp
      call ga_zero(g_rhs_x)   ! zeroed to get the "new" fock matrices

      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*) 'density like matrix before fock_2e'
          call util_flush(LuOut)
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
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
c
      ndensity = ndens*n3xyz
      nfock = ndensity
      if (ndens.gt.1) call uhf_fock_setup(g_tmp, g_rhs_x,
     &            dbl_mb(k_jfac), dbl_mb(k_kfac), ndensity,
     &            nfock, nbf, ocphfprint)
c
      call shell_fock_build(geom,basis,0,ndensity,dbl_mb(k_jfac),
     &            dbl_mb(k_kfac),tol2e,
     &            g_tmp,g_rhs_x,.true.)
c
      if (ndens.gt.1) call  uhf_fock_finish(g_tmp, g_rhs_x,
     &            dbl_mb(k_jfac), dbl_mb(k_kfac), ndensity,
     &            nfock, nbf, ocphfprint)
c
      call ga_dscal(g_rhs_x, multfac)
c
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',1,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
      call schwarz_tidy()
      call int_terminate()
c
      if (.not.ma_pop_stack(l_kfac))
     $  call errquit('hess_cphf: pop problem with l_jfac',555, MA_ERR)
      if (.not.ma_pop_stack(l_jfac))
     $  call errquit('hess_cphf: pop problem with l_kfac',555, MA_ERR)
c
      if (ocphfprint) then
        if(ga_nodeid().eq.0) then
          write(LuOut,*)
          write(LuOut,*) 'fock matrix after shell_fock_build'
          call util_flush(LuOut)
        endif
        call ga_print(g_rhs_x)
      endif
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
#if 0
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
#else
        alo1(1) = ind+1
        ahi1(1) = ind+n3xyz
        blo1(1) = ind+1    
        bhi1(1) = ind+n3xyz 
        blo2(1) = ind+1    
        bhi2(1) = ind+n3xyz
        alo2(1) = ind+1    
        ahi2(1) = ind+n3xyz
        call nga_add_patch(one,g_wax,alo1,ahi1,
     $                     two,g_rhs,blo1,bhi1,
     $                         g_wax,alo1,ahi1)

        call nga_add_patch(one,g_aa,alo2,ahi2,
     $                     one,g_rhs,blo2,bhi2,
     $                         g_aa,alo2,ahi2)
        ind = ind+n3xyz
#endif
      enddo

c
      return
      end
