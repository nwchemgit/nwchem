      subroutine uhf_hessv(acc, g_x, g_ax)
*
* $Id$
*
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"

c     
      double precision acc      ! [input] required accuracy of products
      integer g_x               ! [input] handle to input vectors
      integer g_ax              ! [input] handle to output products
c
      integer gtype, vlen, nvec, nvecp, g_tmp, ivec
      double precision dnrm
      integer ilo(2),ihi(2)
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
     $     call errquit('uhf_hessv: UHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_ax,gtype,vlen,nvec)
      if (vlen.ne.cuhf_vlen)
     $     call errquit('uhf_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,vlen,nvecp)
      if (vlen.ne.cuhf_vlen)
     $     call errquit('uhf_hessv: invalid vector length',0,
     &       UNKNOWN_ERR)
      if (nvecp .ne. nvec) 
     $     call errquit('uhf_hessv: invalid nvecp',nvecp,
     &       UNKNOWN_ERR)
c
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,*)"uhf_hessv: g_x = ",ivec,dnrm
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
         call uhf_hessv_2e(acc, g_x, g_ax, nvec)
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
            write(LuOut,*)"uhf_hessv: g_ax = ",ivec,dnrm
          endif
        enddo
      endif
c
      end
      subroutine uhf_hessv_1e(acc, g_x, g_ax, nvec)
      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "global.fh"
#include "mafdecls.fh"
c     
      double precision acc      ! [input] required accuracy of products
      integer g_x               ! [input] handle to input vectors
      integer g_ax              ! [input] handle to output products
      integer nvec              ! [input] No. of vectors
c
c     For both alpha and beta blocks
c
c     sum(bj) Bai,bj Xbj = 2 sum(b) Fab Xbi - 2 sum(j) Xaj Fji
c
      integer nocc, nvir, g_fock, ivec, g_tmp1, g_tmp2
      integer iend, istart, vbase, iset
      double precision zero, one, two, mtwo
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0, mtwo=-2.0d0)
c     
      do iset = 1, 2
         if (iset .eq. 1) then
            istart = 1
            nocc = nalpha
            g_fock = cuhf_g_falpha
         else
            istart = nalpha*(nmo-nalpha) + 1
            nocc = nbeta
            g_fock = cuhf_g_fbeta
         endif
         nvir  = nmo - nocc
*         write(6,*) ' iset, nocc, nvir', iset, nocc, nvir
*         call ga_print(g_fock)
         iend = istart + nocc*nvir - 1
         vbase = nocc+1         ! First virtual
c
         if (nocc*nvir .gt. 0) then
*ga:1:0
            if (.not. ga_create(MT_DBL, nvir, nocc, 'uhf_hv1e: tmp1',
     $           0, 0, g_tmp1)) call errquit('uhf_hv1e: tmp1', 0,
     &       GA_ERR)
*ga:1:0
            if (.not. ga_create(MT_DBL, nvir, nocc, 'uhf_hv1e: tmp2',
     $           0, 0, g_tmp2)) call errquit('uhf_hv1e: tmp2', 0,
     &       GA_ERR)
c
            do ivec = 1, nvec
               call ga_copy_patch('n',
     $              g_x, istart, iend, ivec, ivec,
     $              g_tmp1, 1, nvir, 1, nocc )
               call ga_zero(g_tmp2)
               call ga_matmul_patch('n', 'n', two, zero, ! 2 Fab Xbi
     $              g_fock, vbase, nmo, vbase, nmo,
     $              g_tmp1, 1, nvir, 1, nocc,
     $              g_tmp2, 1, nvir, 1, nocc )
               call ga_matmul_patch('n', 'n', mtwo, one, ! -2 Xaj Fji
     $              g_tmp1, 1, nvir, 1, nocc, 
     $              g_fock, 1, nocc, 1, nocc,
     $              g_tmp2, 1, nvir, 1, nocc )
               call ga_dadd_patch(
     $              one, g_tmp2, 1, nvir, 1, nocc,
     $              one, g_ax, istart, iend, ivec, ivec,
     $                   g_ax, istart, iend, ivec, ivec)
            enddo
            if (.not. ga_destroy(g_tmp1)) call errquit
     $           ('uhf_hessv_1e: GA corrupt?', 0, GA_ERR)
            if (.not. ga_destroy(g_tmp2)) call errquit
     $           ('uhf_hessv_1e: GA corrupt?', 0, GA_ERR)
         endif
      enddo
c
      end
      subroutine uhf_hessv_2e(acc, g_x, g_ax, nvec)
      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
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
      integer ivec, nfock
      integer g_tmp1, g_tmp2, g_tmp3, g_tmp4
      integer nocc, nvir, g_vecs, nocc_diff, nvir_diff, g_vecs_diff
      integer iset, isetoff, i
      integer istart, istart_diff, iend, iend_diff, vbase, vbase_diff
      integer max_nvec
      parameter (max_nvec=300)
      integer g_dens(4,max_nvec), g_fock(4,max_nvec), g_d, g_j, g_k
      double precision jfac(4,max_nvec), kfac(4,max_nvec)
      double precision tol2e_local
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
c      logical ga_create_atom_blocked
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision zero, one, mone, four
      parameter (zero=0.0d0, one=1.0d0, mone=-1.0d0, four=4.0d0)
c
c     First allocated arrays and build the AO densities.
c     This will be cleaner when we restructure the fock build to take
c     a single density array argument rather than an array of arrays
c
c      write(*,*)"in uhf_hessv, g_x"
c      call ga_print(g_x)
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
         if (.not. ga_create(MT_DBL, nvir, nocc, 'uhf_hv2e: tmp1',
     $        0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0, GA_ERR)
*ga:1:0
         if (.not. ga_create(MT_DBL, nvir, nbf, 'uhf_hv2e: tmp2',
     $        0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
         do ivec = 1, nvec
            g_d = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
            g_j = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
            g_k = ga_create_atom_blocked(geom, basis, 'uhf_h2e: dens')
            call ga_zero(g_d)
            call ga_zero(g_j)
            call ga_zero(g_k)
            g_dens(isetoff+1,ivec) = g_d
            g_dens(isetoff+2,ivec) = g_d
            g_fock(isetoff+1,ivec) = g_j
            g_fock(isetoff+2,ivec) = g_k
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
               call ga_symmetrize(g_d)
            endif
         enddo
         if (.not. ga_destroy(g_tmp1)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_tmp2)) call errquit('uhf_hv2e: GA?',0,
     &       GA_ERR)
      enddo
c
c     Now make all the J & K matrices with one pass thru the integrals
c
c     NOTE that symmetrization is not yet implemented except for
c     totally symmetric products.  Assume for the time being that
c     if only 1 RHS is being requested then use symmtery, but disable
c     it if there is more than 1 RHS.
c
      do ivec = 1, nvec
         jfac(1,ivec) = one
         jfac(2,ivec) = zero
         jfac(3,ivec) = one
         jfac(4,ivec) = zero
         kfac(1,ivec) = zero
         kfac(2,ivec) = one
         kfac(3,ivec) = zero
         kfac(4,ivec) = one
      enddo
      oskel_local = oskel .and. (nvec.eq.1)
c
      tol2e_local = min(max(acc,itol_floor),itol_ceil)
c
      nfock = 4*nvec
c
      call fock_2e(geom, basis, nfock, jfac, kfac,
     $     tol2e_local, oskel_local, g_dens, g_fock, .false.)
c
c     Again, note that this is NOT what you want except for
c     a totally symmetric input vector
c
      if (oskel_local) then
         if (oscfps) call pstat_on(ps_sym_sym)
         do ivec = 1, nvec
            do iset = 1, 4
               call sym_symmetrize(geom, basis, .false.,
     $              g_fock(iset,ivec))
            enddo
         enddo
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c
      do ivec = 1, nvec         ! Don't need the AO densities any more
         if (.not. ga_destroy(g_dens(1,ivec)))
     $        call errquit('uhf_hessv: GA corrupt?',0, GA_ERR)
         if (.not. ga_destroy(g_dens(3,ivec)))
     $        call errquit('uhf_hessv: GA corrupt?',0, GA_ERR)
      enddo
c
c     Now transform to the MO basis and assemble the results
c
      do iset = 1, 2
         if (iset .eq. 1) then
            isetoff = 0
            istart = 1
            nocc = nalpha
            g_vecs = g_movecs(1)
            istart_diff = nalpha*(nmo-nalpha) + 1
            nocc_diff = nbeta
            g_vecs_diff = g_movecs(2)
         else
            isetoff = 2
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
*ga:1:0
        if (.not. ga_create(MT_DBL, nbf, nocc, 'uhf_hv2e: tmp1',
     $       0, 0, g_tmp1)) call errquit('uhf_hv2e: tmp1', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL, nvir, nocc, 'uhf_hv2e: tmp2',
     $       0, 0, g_tmp2)) call errquit('uhf_hv2e: tmp2', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL, nbf, nocc_diff, 'uhf_hv2e: tmp3',
     $       0, 0, g_tmp3)) call errquit('uhf_hv2e: tmp3', 0, GA_ERR)
*ga:1:0
        if (.not. ga_create(MT_DBL,nvir_diff,nocc_diff,'uhf_hv2e: tmp4',
     $       0, 0, g_tmp4)) call errquit('uhf_hv2e: tmp4', 0, GA_ERR)
        do ivec = 1, nvec
           g_j = g_fock(isetoff+1,ivec)
           g_k = g_fock(isetoff+2,ivec)
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
     $                four, g_tmp2, 1, nvir, 1, nocc,
     $                one, g_ax, istart, iend, ivec, ivec,
     $                g_ax, istart, iend, ivec, ivec)
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
c     write(*,*)"g_ax: differet spin, iset = ", iset
c              call ga_print(g_ax)
           endif
 200       continue
           if (.not. ga_destroy(g_j)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
           if (.not. ga_destroy(g_k)) call errquit('uhf_hessv: GA?',0,
     &       GA_ERR)
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
      end

           
