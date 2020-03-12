      subroutine uhf_hessv2(acc, g_x, g_ax)
*
* $Id$
*
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "global.fh"

c     
      double precision acc      ! [input] required accuracy of products
      integer g_x               ! [input] handle to input vectors
      integer g_ax              ! [input] handle to output products
c
      integer gtype, vlen, nvec, nvecp, ivec
      double precision dnrm
      integer ilo(2), ihi(2)
c
      logical oprint, olprint
c
      oprint = util_print("hessv",print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)
c
c     Multiply a set of vectors by the level-shifted UHF hessian.
c     
c     Check dimensions
c     
c      if (cuhf_init_flag.ne.1)
      if(.not.cuhf_init_flag)
     $     call errquit('uhf_hessv2: UHF internal block invalid',0,
     &       INPUT_ERR)
      call ga_inquire(g_ax,gtype,vlen,nvec)
      if (vlen.ne.cuhf_vlen)
     $     call errquit('uhf_hessv2: invalid vector length',0, GA_ERR)
      call ga_inquire(g_x,gtype,vlen,nvecp)
      if (vlen.ne.cuhf_vlen)
     $     call errquit('uhf_hessv2: invalid vector length',0, GA_ERR)
      if (nvecp .ne. nvec) 
     $     call errquit('uhf_hessv2: invalid nvecp',nvecp, INPUT_ERR)
c
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,"uhf_hessv2: g_x = ",i4,f24.8)')
     +      ivec,dnrm
          endif
        enddo
      endif
c
      call ga_zero(g_ax)
*      write(6,*) ' lshift ', lshift
      if (lshift .ne. 0.0d0)
     $     call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
c     
      call uhf_hessv_1e(acc, g_x, g_ax, nvec)
c     
 200  continue
c
      if (pflg .gt. 1)then 
         call uhf_hessv_2e2(acc, g_x, g_ax, nvec)
      endif
c
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_ax,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,"uhf_hessv2: g_ax = ",i4,f24.8)')
     +      ivec,dnrm
          endif
        enddo
      endif
c
      end
c
      subroutine uhf_hessv_2e2(acc, g_x, g_ax, nvec)
      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "case.fh"
c     
      double precision acc      ! [input] required accuracy of products
      integer g_x               ! [input] handle to input vectors
      integer g_ax              ! [input] handle to output products
      integer nvec              ! [input] No. of vectors
c
c     SAME spin -> Bai,bj Xbj = [4 (ai|bj) - 2 (bi|aj) - 2 (ij|ab)] Xbj
c                             = 4 (Jai - Kai)
c
c     where J and K constructed from symmetrized AO density CXCT and
c     transformed into MO basis of the SAME spin.
c   
c     DIFF spin -> Bai,bj Xbj = 4 (ai|bj) Xbj
c                             = 4 Jai
c     
c     where J constructed from  symmetrized AO density CXCT and J 
c     transformed into MO basis of the OPPOSITE spin.
c
      logical oskel_local
      integer ivec, nfock, ndensity
      integer g_tmp1, g_tmp2, g_tmp3, g_tmp4
      integer nocc, nvir, g_vecs, nocc_diff, nvir_diff, g_vecs_diff
      integer nnocc, nnocc_diff
      integer iset, isetoff, isetoff2
      integer istart, istart_diff, iend, iend_diff, vbase, vbase_diff
      integer g_d, g_j, g_k
      integer g_dens2, g_fock2, g_work2
      integer l_jfac, k_jfac, l_kfac, k_kfac
      integer alo(3), ahi(3), blo(3), bhi(3), dims(3), chunk(3)
      double precision tol2e_local
      double precision itol_floor, itol_ceil
      double precision xc_xfac1
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
c      logical ga_create_atom_blocked
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
c     DIM/QM
      integer g_dens, g_dim
      logical ldimqm
c
c     First allocated arrays and build the AO densities.
c     This will be cleaner when we restructure the fock build to take
c     a single density array argument rather than an array of arrays
c
c     write(6,*)"in uhf_hessv2, g_x"
c     call ga_print(g_x)
      xc_xfac1=bgj_kfac() 
c
c     allocate MAs for jfacs and kfacs
c
      if (.not.MA_Push_Get(mt_dbl,nvec*4,'hessv jfacs',
     &     l_jfac,k_jfac))
     &     call errquit('hessv2: cannot allocate jfac',nvec*4, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,nvec*4,'hessv kfacs',
     &     l_kfac,k_kfac))
     &     call errquit('hessv2: cannot allocate kfac',nvec*4, MA_ERR)
c
c     initialize jfacs and kfacs (will go into the uhf_fock_setup)
c
      do iset = 0, 2*nvec-1
          if (noskew_uhf) then
             dbl_mb(k_jfac+iset) =  1.0d0
          else
             dbl_mb(k_jfac+iset) =  0.0d0
          endif
          dbl_mb(k_kfac+iset) =  xc_xfac1
      enddo
c
c     Allocate arrays and generate the AO densities
c
      dims(1) = 2*nvec
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens',chunk,
     &     g_dens2)) call errquit
     &    ('uhf_hessv_2e2: could not allocate g_dens2',555, GA_ERR)
      call ga_zero(g_dens2)
      g_d = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
      g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
c
      alo(1) = 0
      ahi(1) = 0
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
c
      do iset = 1, 2
         if (iset .eq. 1) then
            isetoff = 0
            istart = 1
            nocc = nalpha
         else
            isetoff = 2
            istart = nalpha*(nmo-nalpha) + 1
            nocc = nbeta
         endif
         g_vecs = g_movecs(iset)
         nvir  = nmo - nocc
         iend = istart + nocc*nvir - 1
         vbase = nocc+1         ! First virtual
*ga:1:0
         nnocc = max(1,nocc) ! nocc can be 0 in UHF 1-electron systems
         if (.not. ga_create(MT_DBL, nvir, nnocc, 'uhf_hv2e: tmp1',
     $        0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0, GA_ERR)
*ga:1:0
         if (.not. ga_create(MT_DBL, nvir, nbf, 'uhf_hv2e: tmp2',
     $        0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
         do ivec = 1, nvec
            call ga_zero(g_d)
            alo(1) = alo(1) + 1
            ahi(1) = ahi(1) + 1
            if (nocc*nvir .gt. 0) then
               call ga_copy_patch('n', ! Reshape vector into matrix Xbj
     $              g_x, istart, iend, ivec, ivec,
     $              g_tmp1, 1, nvir, 1, nocc)
               call ga_matmul_patch('n', 't', one, zero,
     $              g_tmp1, 1, nvir, 1, nocc,
     $              g_vecs, 1, nocc, 1, nbf,
     $              g_tmp2, 1, nvir, 1, nbf)
               call ga_matmul_patch('n', 'n', one, zero,
     $              g_vecs, 1, nbf, vbase, nmo,
     $              g_tmp2, 1, nvir, 1, nbf,
     $              g_d, 1, nbf, 1, nbf)
               if (noskew_uhf) then
                  call ga_symmetrize(g_d)
               else
                  call ga_antisymmetrize(g_d)
               endif
               call nga_copy_patch('N',g_d,blo,bhi,g_dens2,alo,ahi)
            endif
         enddo
         if (.not.ga_destroy(g_tmp1)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
         if (.not.ga_destroy(g_tmp2)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
      enddo
      if (.not.ga_destroy(g_d)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
c
c     DIM/QM JEM
c     Calculate DIM dipoles and DIM potential from TOTAL density
      if (.not.rtdb_get(bgj_get_rtdb_handle(), 'dimqm:lrsp', mt_log,
     $                  1, ldimqm)) ldimqm = .false.
      if (ldimqm) then
        dims(1) = nvec
        chunk(1) = dims(1)
        if (.not. nga_create(MT_DBL, 3, dims, 'tdens', chunk, g_dens))
     $    call errquit('uhf_hessv2: DIM dens create', 0, GA_ERR)
c        if (.not. nga_create(MT_DBL, 3, dims, 'g_dim', chunk, g_dim))
c     $    call errquit('uhf_hessv2: DIM pot create', 0, GA_ERR)
        alo(1) = 1
        ahi(1) = nvec
c       Copy Alpha
        call nga_copy_patch('n', g_dens2, alo, ahi, g_dens, alo, ahi)
        alo(1) = nvec+1
        ahi(1) = 2*nvec
        blo(1) = 1
        bhi(1) = nvec
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = 1
        bhi(3) = nbf
c       Add Beta
c       Seems like we need to divide both alpha and beta by 4
c       to match closed shell routine?
        call nga_add_patch(0.25d0, g_dens2, alo, ahi,
     $                     0.25d0, g_dens,  blo, bhi,
     $                             g_dens,  blo, bhi)
c jbecca START: added .false. variable at end
        call dimqm_rohf_calcDIM(bgj_get_rtdb_handle(), geom, basis,
     $                          nbf, nvec, g_dens, g_dim, .false.)
c jbecca END
c        call dimqm_indDipoles(bgj_get_rtdb_handle(), nbf, geom,
c     $                          basis, g_dens)
c        call ga_zero(g_dim)
c        call fock_dim(geom, nbf, basis, nvec, g_dim, 0, 1)
c        call ga_symmetrize(g_dim)
        if(.not. ga_destroy(g_dens))
     $    call errquit('uhf_hessv2: DIM destroy dens', 0, GA_ERR)
c       Reset dimensions for the rest of the subroutine
        dims(1) = 2*nvec
        chunk(1) = dims(1)
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        bhi(3) = nbf
      end if
c
c     Now make all the J & K matrices with one pass thru the integrals
c
c     NOTE that symmetrization is not yet implemented except for
c     totally symmetric products.  Assume for the time being that
c     if only 1 RHS is being requested then use symmtery, but disable
c     it if there is more than 1 RHS.
c
      oskel_local = oskel .and. (nvec.eq.1)
c
      tol2e_local = min(max(acc,itol_floor),itol_ceil)
c
c     Get Fock matrices
c
      if (.not. nga_create (MT_DBL, 3, dims, 'Fockv',chunk,
     &     g_fock2)) call errquit
     &    ('uhf_hessv_2e2: could not allocate g_fock2',555, GA_ERR)
      call ga_zero(g_fock2)
c
      if (.not. nga_create (MT_DBL, 3, dims, 'work2',chunk,
     &     g_work2)) call errquit
     &    ('uhf_hessv_2e2: could not allocate g_work2',555, GA_ERR)
      call ga_zero(g_work2)
c
c     Do some setup
c
      ndensity = 2*nvec
      nfock = ndensity
      call uhf_fock_setup(g_dens2, g_fock2,
     &            dbl_mb(k_jfac), dbl_mb(k_kfac), ndensity,
     &            nfock, nbf, .false.)
c
c     Actually do the Fock build
c
      if (.not.cam_exch) then
        call shell_fock_build(geom, basis, 0, nfock,
     $     dbl_mb(k_jfac), dbl_mb(k_kfac),
     $     tol2e, g_dens2, g_fock2, noskew_uhf)
      else
        call shell_fock_build_cam(geom, basis, 0, nfock,
     $     dbl_mb(k_jfac), dbl_mb(k_kfac),
     $     tol2e, g_dens2, g_fock2, noskew_uhf)
      end if ! cam_exch
      !call ga_print(g_dens2)

c
      if (.not.ma_pop_stack(l_kfac))
     $  call errquit('uhf_hessv_2e2: pop problem with l_jfac',555,
     &       MA_ERR)
      if (.not.ma_pop_stack(l_jfac))
     $  call errquit('uhf_hessv_2e2: pop problem with l_kfac',555,
     &       MA_ERR)
c
      if (.not. ga_destroy(g_dens2))
     $     call errquit('uhf_hessv_2e2: GA corrupt?',0, GA_ERR)
      if (.not. ga_destroy(g_work2)) 
     &     call errquit('uhf_hessv_2e2: GA corrupt?',0, GA_ERR)
c
c     Now transform to the MO basis and assemble the results
c     I would like this to eventually be part of uhf_fock_finish
c     -TLW
c
      do iset = 1, 2
         if (iset .eq. 1) then
            isetoff = 0
            isetoff2 = nvec
            istart = 1
            nocc = nalpha
            g_vecs = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff = 2*nvec
            isetoff2 = 3*nvec
            istart = nalpha*(nmo-nalpha) + 1
            nocc = nbeta
            g_vecs = g_movecs(2)
            istart_diff = 1
            nocc_diff = nalpha
            g_vecs_diff = g_movecs(1)
        endif
        nvir  = nmo - nocc
        nvir_diff = nmo - nocc_diff
        iend = istart + nocc*nvir - 1
        iend_diff = istart_diff + nocc_diff*nvir_diff - 1
        vbase = nocc+1          ! First virtual
        vbase_diff = nocc_diff+1          ! First virtual
        nnocc = max(1,nocc)
        nnocc_diff = max(1,nocc_diff)
*ga:1:0
        if (.not. ga_create(MT_DBL, nbf, nnocc, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0,
     &       GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL, nvir, nnocc, 'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL, nbf, nnocc_diff, 'uhf_hv2e: tmp3',
     $       0, 0, g_tmp3)) call errquit('uhf_hv2e: tmp3', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL,nvir_diff,nnocc_diff,
     $       'uhf_hv2e: tmp4',
     $       0, 0, g_tmp4)) call errquit('uhf_hv2e: tmp4', 0, GA_ERR)
        do ivec = 1, nvec
           alo(1) = isetoff + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock2,alo,ahi,g_j,blo,bhi)
           alo(1) = isetoff2 + ivec
           ahi(1) = alo(1)
           call nga_copy_patch('N',g_fock2,alo,ahi,g_k,blo,bhi)
c     
c     Same spin 2-e contributions
c     
           if (nocc*nvir .gt. 0) then
              call ga_dadd(one, g_j, mone, g_k, g_k) ! K <-- J-K
              call ga_matmul_patch('n', 'n', one, zero,
     $             g_k, 1, nbf, 1, nbf,
     $             g_vecs, 1, nbf, 1, nocc,
     $             g_tmp1, 1, nbf, 1, nocc)
              call ga_matmul_patch('t', 'n', one, zero,
     $             g_vecs, vbase, nmo, 1, nbf,
     $             g_tmp1, 1, nbf, 1, nocc,
     $             g_tmp2, 1, nvir, 1, nocc)
              call ga_dadd_patch(
     $             four, g_tmp2, 1, nvir, 1, nocc,
     $             one, g_ax, istart, iend, ivec, ivec,
     $             g_ax, istart, iend, ivec, ivec)
           endif
c     
c     Different spin
c
           if (nocc_diff*nvir_diff .gt. 0) then
              call ga_matmul_patch('n', 'n', one, zero,
     $             g_j, 1, nbf, 1, nbf,
     $             g_vecs_diff, 1, nbf, 1, nocc_diff,
     $             g_tmp3, 1, nbf, 1, nocc_diff)
              call ga_matmul_patch('t', 'n', one, zero,
     $             g_vecs_diff, vbase_diff, nmo, 1, nbf,
     $             g_tmp3, 1, nbf, 1, nocc_diff,
     $             g_tmp4, 1, nvir_diff, 1, nocc_diff)
              call ga_dadd_patch(
     $             four, g_tmp4, 1, nvir_diff, 1, nocc_diff,
     $             one, g_ax, istart_diff, iend_diff, ivec, ivec,
     $                  g_ax, istart_diff, iend_diff, ivec, ivec)
           endif
c
c     DIM/QM JEM
c
           if (ldimqm) then
             alo(1) = ivec
             ahi(1) = alo(1)
             call nga_copy_patch('N', g_dim, alo, ahi, g_k, blo, bhi)
c            Mulitply DIM potential by 4 due to earlier divide.
             call ga_matmul_patch('n', 'n', 4.0d0, zero,
     $             g_k, 1, nbf, 1, nbf,
     $             g_vecs, 1, nbf, 1, nocc,
     $             g_tmp1, 1, nbf, 1, nocc)
              call ga_matmul_patch('t', 'n', one, zero,
     $             g_vecs, vbase, nmo, 1, nbf,
     $             g_tmp1, 1, nbf, 1, nocc,
     $             g_tmp2, 1, nvir, 1, nocc)
              call ga_dadd_patch(
     $             four, g_tmp2, 1, nvir, 1, nocc,
     $             one, g_ax, istart, iend, ivec, ivec,
     $             g_ax, istart, iend, ivec, ivec)

           end if
 200       continue
        enddo
        if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp3)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp4)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      enddo
c
      if (.not. ga_destroy(g_j)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_k)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock2)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
      if (ldimqm) then
        if (.not. ga_destroy(g_dim))
     $    call errquit('uhf_hessv2 destroy g_dim failed', 1, GA_ERR)
      end if
c
      end
