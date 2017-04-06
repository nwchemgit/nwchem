      subroutine dft_gradients_so(rtdb)
c     $Id$
c
c calculate energy gradients with respect to nuclear coordinates
c modified from scf version for use in DFT gradients
c
c------------------------------------------------------------------------------
c         ___                 ___                         ___
c dE      \         dh(i,j)   \             d(mn|op)      \          dS(i,j)
c -- = 2   > D(i,j) ------- +  > P(m,n,o,p) -------- - 2   > Dw(i,j) -------
c dA      /           dA      /                dA         /            dA
c         ---                 ---                         ---
c         i,j                ijkl                         i,j
c
c
c        dV(nuc-nuc)
c      + ----------     + exchange-correlation terms
c            dA
c
c------------------------------------------------------------------------------
c                                1
c P(i,j,k,l) = [2 D(i,j)D(k,l) - - (D(i,k)D(j,l) + D(i,l)D(j,k))
c                                2
c------------------------------------------------------------------------------
c
c     This version computes the pieces specific to DFT (XC on grid
c     and CD-fit) and call the standard grad_force() to do the rest
c     including writing out the results.
c
      implicit none
#include "errquit.fh"
c      
      integer rtdb
c
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "util.fh"
#include "sym.fh"
#include "cdft.fh"
c
c!!! BGJ test !!!
#include "bgj.fh"
      integer l_hess, k_hess, g_rhs(3,100), j
c!!! BGJ test !!!
      double precision  zero, one, two
      Parameter (zero=0.d0, one=1.d0, two=2.d0)
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_header_so, movecs_read_so
      external movecs_read_header_so, movecs_read_so
      integer nmo(2)
      integer iga_dens(2), g_vecs(2), g_force, g_densso(2)
      Integer k_dens, l_dens,
     $     k_evals(2), l_evals(2)
      double precision grad_norm, grad_max
      external grad_norm, grad_max
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      logical status
      integer me, nproc, max_sh_bf, max_at_bf, nat, max_sh_bfcd,
     $     nocc, lforce, nactive, i, nbf_vecs,
     $     max1e, max2e, mscratch_1e, mscratch_2e, 
     $     max2e3c, mscratch_2e3c, lbuf, lscratch, lsqa
      integer l_force, k_force, l_occ, k_occ, l_act, k_act, 
     $     l_buf, k_buf, l_scr, k_scr, l_wdens, k_wdens, 
     $     l_cdcoef, i_cdcoef, ippp, isvec, lsvec,
     $     ilo, ihi, 
     $     k_frc_2el, k_frc_xc, 
     $     l_frc_2el, l_frc_xc
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     $     lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     $     nscr, lscr, iscr, nTotEl 
      double precision charge, charge_nuc, rhffact, tol2e, onem
c     
      nproc = ga_nnodes()
      me=ga_nodeid()
c     
c     Print options
c     
      if (.not. geom_ncent(geom, nat)) 
     $     call errquit('dft_gradient: could not get natoms',0,
     &       GEOM_ERR)
c     
      if (.not. bas_nbf_cn_max(ao_bas_han, max_sh_bf))
     $     call errquit('dft_gradient: could not get max_sh_bf',0,
     &       BASIS_ERR)
      max_at_bf = 0
      do i = 1, nat
         if (.not. bas_ce2bfr(ao_bas_han, i, ilo, ihi))
     $        call errquit('dft_gradient: bas_ce2bfr failed', i,
     &       BASIS_ERR)
         max_at_bf = max(max_at_bf, ihi-ilo+1)
      enddo
c     
c     use of scratch array in cdfit ... needs (3,max_at_bf)
c     
      max_at_bf = max(max_at_bf,3)
c     
      charge = rcharge
      status = geom_nuc_charge(geom, charge_nuc)
      if (.not.status)then
         call errquit('dft_gradient: no nuclear charge',0, GEOM_ERR)
      endif 
c     
      if (ipol .eq. 1)then
         nocc = max(noc(1),noc(2))
      endif 
c     
c     allocate and initialize global and local memory
c     
c     mo-vectors
c     
      if (ipol .eq. 1)then
         g_vecs(1) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'mo vectors')
      else 
         status = ga_create(mt_dbl, 2*nbf, 2*nbf,'Movecs Re', 0, 0, 
     &        g_vecs(1))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' dft_gradients_so: ', 'Movecs Re'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('dft_gradients_so: ga_create ', 0, GA_ERR)
         endif
         
         status = ga_create(mt_dbl, 2*nbf, 2*nbf,'Movecs Im', 0, 0, 
     &        g_vecs(2))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' dft_gradients_so: ', 'Movecs Im'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('dft_gradients_so: ga_create ', 0, GA_ERR)
         endif
         
         call ga_zero(g_vecs(1))
         call ga_zero(g_vecs(2))
      endif 
      call ga_sync() 
c     
c     global density
c     
      if (ipol .eq. 1)then
         rhffact = two
      else 
         rhffact = one
      endif 
      if (ipol .eq. 1)then
         iga_dens(1) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'density matrix')
      else 
         status = ga_create
     $        (mt_dbl, nbf, nbf,'alpha density matrix', 0, 0, 
     &        iga_dens(1))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' dft_gradients_so: ', 
     $              'alpha density matrix'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('dft_gradients_so: ga_create ', 0, GA_ERR)
         endif
         
         
         status = ga_create
     $        (mt_dbl, nbf, nbf,'beta density matrix', 0, 0, 
     &        iga_dens(2))
         if (.not. status) then
            if (ga_nodeid() .eq. 0) then
               write(6,*) ' dft_gradients_so: ', 
     $              'beta density matrix'
               call util_flush(6)
            endif
            call ga_sync()
            call errquit('dft_gradients_so: ga_create ', 0, GA_ERR)
         endif
      endif 

      status = ga_create(mt_dbl, 2*nbf, 2*nbf,'DenMx Re', 0, 0, 
     &                   g_densso(1))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' dft_gradients_so: ', 'DenMx Re'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('dft_gradients_so: ga_create ', 0, GA_ERR)
      endif


      status = ga_create(mt_dbl, 2*nbf, 2*nbf,'DenMx Im', 0, 0, 
     &                   g_densso(2))
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' dft_gradients_so: ', 'DenMx Im'
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('dft_gradients_so: ga_create ', 0, GA_ERR)
      endif

c     
c     forces on atoms (3xnat)
c     
      status = ga_create (mt_dbl, 3, nat, 'forces', 3, 0, g_force)
      call ga_zero (g_force)
c     
c     local replication (separate for the different pieces)
c     
      lforce = nat * 3
      if (.not.ma_alloc_get(mt_dbl, lforce, 'forces',l_force, k_force))
     $     call errquit('could not allocate l_force',1, MA_ERR)
      call dfill(lforce, 0.0d0, dbl_mb(k_force), 1)
c     
      if (.not.ma_alloc_get(mt_dbl,lforce,'forces',l_frc_2el,k_frc_2el))
     $     call errquit('could not allocate l_frc_2el',1, MA_ERR)
      call dfill(lforce, 0.0d0, dbl_mb(k_frc_2el), 1)
c     
      if (.not.ma_alloc_get(mt_dbl,lforce,'forces',l_frc_xc,k_frc_xc))
     $     call errquit('could not allocate l_frc_xc',1, MA_ERR)
      call dfill(lforce, 0.0d0, dbl_mb(k_frc_xc), 1)
c     
c     eigenvalues
c     
      status = ma_alloc_get(mt_dbl, 2*nbf_ao, 'alpha MO evals',
     $     l_evals(1), k_evals(1))
c     
c     occupation numbers (not used, but necessary for movecs_read_so)
c     
c     should do k_occ for both spins, in case used at some point...
c     
      if (.not. ma_alloc_get(mt_dbl, 2*nbf_ao, 'occ. numbers',
     $     l_occ, k_occ))
     $     call errquit('dft_gradient: could not allocate l_occ',1,
     &       MA_ERR)
c     
c     lookup table and list of active atoms
c     
      if (.not. ma_alloc_get(MT_LOG, nat, 'active atoms',
     $     l_act, k_act))
     $     call errquit('grad: could not allocate l_act',1,
     &       MA_ERR)
      call grad_active_atoms(rtdb, nat, log_mb(k_act), nactive)
c     
c     get MO vectors from file
c     
c      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
c     $     call errquit('dft_gradient: DFT MO vectors not defined',0)
      status = movecs_read_header_so
     $     (movecs_in, title_vecs, basis_vecs,
     $     scftype_vecs, nbf_vecs)
c     
c     Should check much more info than just nbf for consistency
c     
c     
c     get mo eigevectors
c     
      if (2*nbf_ao .ne. nbf_vecs)then
         write(6,*)'dft_gradient movecs output = ',movecs_in
         call errquit('dft_gradient: could not read mo vectors',911,
     &       DISK_ERR)
      else 
         status = .true.
         status = status .and.
     $        movecs_read_so(movecs_in, dbl_mb(k_occ),
     $        dbl_mb(k_evals(1)), g_vecs)
      endif
c     
      if (.not.status)then
         write(6,*)'dft_gradient movecs output = ',movecs_in
         call errquit('dft_gradient: could not read mo vectors',917,
     &       DISK_ERR)
      endif
 
      call ga_zero(g_densso(1))
      call ga_zero(g_densso(2))
cso
cso   Form Re and Im of density matrix
cso
      nTotEl = noc(1)+noc(2) 
      call dft_densm_so(g_densso, g_vecs, nbf_ao, nTotEl)
      call ga_sync()

      call ga_zero(iga_dens(1))
      call ga_zero(iga_dens(2))
      call ga_dens_sf(iga_dens, g_densso, nbf_ao)
      call ga_sync

      call ga_symmetrize(iga_dens(1))
      call ga_symmetrize(iga_dens(2))

c     
c     Pre-compute mapping vectors
c     
      if (.not.ma_push_get
     $     (mt_int,nat*2,'cntoce map',lcetobfr,icetobfr))
     $     call errquit('dft_scf:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     $     (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     $     call errquit('dft_scf:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     $     (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     $     call errquit('dft_scf:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     $     int_mb(icetobfr), nat, nshells_ao)
      if (.not.ma_pop_stack(lcntobfr))
     $     call errquit('dft_gradient: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcntoce))
     $     call errquit('dft_gradient: cannot pop stack',0, MA_ERR)
c     
c     Pre-compute reduced density matrices over atoms
c     
      if (.not.ma_push_get(mt_dbl,ipol*nat*nat,'rdens_atom',
     $     lrdens_atom,irdens_atom))
     $     call errquit('dft_scf: cannot allocate rdens_atom',0, MA_ERR)
      call dfill(ipol*nat*nat, 0.0d0, dbl_mb(irdens_atom), 1)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.ma_push_get(mt_dbl,nscr,'scr',lscr,iscr))
     $     call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
      call util_ga_mat_reduce(nbf_ao, nat, int_mb(icetobfr), 
     $     iga_dens, ipol, dbl_mb(irdens_atom), 
     $     'absmax', dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
      if (.not.ma_pop_stack(lscr))
     $     call errquit('dft_scf: cannot pop stack',0, MA_ERR)
c     
      status = ga_destroy(g_vecs(2))
      status = ga_destroy(g_vecs(1))
c     
      if (.not.status)then
         call errquit('dft_gradient: could not destroy g_eigen_diag',1,
     &       GA_ERR)
      endif 
c     
      status = ma_free_heap(l_occ)
      status = ma_free_heap (l_evals(1))
c     
c     determine memory requirements for integral gradients
c     
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call int_mem_2e3c(max2e3c, mscratch_2e3c)
      lbuf = max(max1e, max2e)
      if (CDFIT) lbuf = max(lbuf, max2e3c) + 500
      lscratch = max(mscratch_1e, mscratch_2e)
      if (CDFIT) lscratch = max(lscratch, mscratch_2e3c)
c     
c     one-electron contribution
c     buffers for one electron integral derivatives
c     
      status = ma_push_get(mt_dbl, lbuf, 'deriv buffer', l_buf, k_buf)
      if (.not.status) 
     $     call errquit('dft_gradient: could not allocate buffer',1,
     &       MA_ERR)
c     
      status = ma_push_get(mt_dbl, lscratch, 'deriv scratch', 
     $     l_scr, k_scr)
      if (.not.status) 
     $     call errquit('dft_gradient: could not allocate scratch',1,
     &       MA_ERR)
c     
c     allocate local density matrix block
c     
      lsqa = max_at_bf * max_at_bf
c     
      status = ma_push_get(mt_dbl, lsqa, 'local_w_density', 
     $     l_wdens, k_wdens)
      status = status .and. ma_push_get(mt_dbl, lsqa, 'local_density',
     $     l_dens, k_dens)
      if (.not.status)call errquit('could not allocate l_dens',1,
     &       MA_ERR)
c     
c     store total DM in ga_dens(1)
c     
      call ga_dadd (one,iga_dens(1),one,iga_dens(2),iga_dens(1))
c     
c     define threshold for Schwarz screening (same as in DFT)
c     
      tol2e=10.d0**(-itol2e)
c     
c     zero force component vector
c     
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_2el), 1)
c     
      IF (CDFIT)then
c     
c     compute 3 center coulomb derivative
c     
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c     
         if (.not. ma_push_get(mt_dbl, nbf_cd, 'CD coef',
     $        l_cdcoef, i_cdcoef))
     $        call errquit('dft_gradient: could not alloc CD coef',0,
     &       MA_ERR)
c     
         if(.not.bas_nbf_cn_max(cd_bas_han, max_sh_bfcd))
     $        call errquit('dftg_force: basnbfcdmax broken?',0,
     &       BASIS_ERR)
         if (.not. MA_Push_get(MT_DBL, 3*max_sh_bfcd*2, 'svec',
     $        lsvec, isvec))
     $        call errquit('dftg_force: could not alloc svec',0, MA_ERR)
         ippp=k_wdens
c     
         call dftg_cdfit(geom,ao_bas_han, cd_bas_han,
     $        nbf_cd, nat, tol2e, dbl_mb(k_scr),
     $        lscratch, dbl_mb(k_buf), lbuf,
     $        dbl_mb(isvec), dbl_mb(ippp), max_sh_bf,
     $        iga_dens, dbl_mb(k_frc_2el),
     $        DBL_MB(i_cdcoef), oskel)
c     
         call ga_dgop(msg_grad_2el, dbl_mb(k_frc_2el), lforce, '+')
      endif
c     
      if (.not.ma_chop_stack(l_buf))
     $     call errquit('dft_gradient: cannot chop stack',0, MA_ERR)
c     
c     restore alpha DM in g_dens(1)
c     
      onem = -1.d0
      call ga_dadd(one, iga_dens(1), onem, iga_dens(2), iga_dens(1))
c     
c     get exchange-correlation contribution to the gradient
c     
c     zero force component vector
c     
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_xc),  1)
c     
      call dftg_getxc(rtdb, nat, iga_dens, dbl_mb(k_frc_xc), 
     $     log_mb(k_act), nactive, 
     $     dbl_mb(irdens_atom), int_mb(icetobfr))
      call ga_dgop(msg_grad_xc,  dbl_mb(k_frc_xc),  lforce, '+')
c     
c     
      if (ga_nodeid() .eq. 0)then
         status = rtdb_parallel (.false.)
         do i = 0, lforce-1
            dbl_mb(k_force+i) = dbl_mb(k_frc_2el+i) + 
     $           dbl_mb(k_frc_xc+i) 
         enddo
         if (.not. rtdb_put(rtdb, 'dft:cd+xc gradient', mt_dbl,
     $        lforce, dbl_mb(k_force))) call errquit
     $        ('dft_gradient: failed storing cd+xc gradient',0,
     &       RTDB_ERR)
c     
         status = rtdb_parallel (.true.)
      endif
c     
      call ga_sync()
c     
      if (.not.ma_pop_stack(lrdens_atom))
     $     call errquit('dft_gradient: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcetobfr))
     $     call errquit('dft_gradient: cannot pop stack',0, MA_ERR)
c     
c!!! BGJ test !!!
c
c     store total DM in ga_dens(1)
c
      call ga_dadd (1d0,iga_dens(1),1d0,iga_dens(2),iga_dens(1))
c
c     J hesssian test calculation done by setting bgj:j_hessian
c     to true
c
      if (.not. rtdb_get(rtdb, 'bgj:j_hessian', mt_log,
     &     1, status)) status = .false.
      if (status) then
         call schwarz_tidy()
         call intd_terminate()
         call int_init(rtdb, 1, ao_bas_han)
         call schwarz_init (geom, ao_bas_han)
         call int_terminate()
         if (CDFIT) then
           nmo(1) = ao_bas_han
           nmo(2) = cd_bas_han
           call intdd_init(rtdb,2,nmo)
         else
           call intdd_init(rtdb,1,ao_bas_han)
         endif

         status = MA_push_get(MT_DBL, 3*nat*3*nat,
     &        'j hessian', l_hess, k_hess)
         if (.not.status) 
     &        call errquit('dft_gradients_so: could not alloc j 
     &        hessian', 1, MA_ERR)
         if (bgj_print() .gt. 0)
     &        write(6,*)'*** In dft_gradients_so: calling j_hessian'
         call j_hessian(iga_dens, log_mb(k_act), nactive,
     &        dbl_mb(k_hess))
         status = MA_pop_stack(l_hess)
         if (.not.status) call
     &        errquit('dft_gradients_so: could not pop j hessian',
     &        1, MA_ERR)

         call schwarz_tidy()
         call intdd_terminate()

      endif
c
c     J CPKS RHS test calculation done by setting bgj:j_cpks_rhs
c     to true
c
      if (.not. rtdb_get(rtdb, 'bgj:j_cpks_rhs', mt_log,
     &     1, status)) status = .false.
      if (status) then

         call schwarz_tidy()
         call intd_terminate()

         call int_init(rtdb, 1, ao_bas_han)
         call schwarz_init (geom, ao_bas_han)
         call int_terminate()
         if (CDFIT) then
           nmo(1) = ao_bas_han
           nmo(2) = cd_bas_han
           call intd_init(rtdb,2,nmo)
         else
           call intd_init(rtdb,1,ao_bas_han)
         endif
c        !!! Do this to be consistent with DFT gradient
c!!!         call int_app_set_no_texas(rtdb)
c
c        Allocate temp GA's for RHS
c
         if (bgj_print() .gt. 0)
     &        write(*,*)'*** j cpks rhs test: nactive =',nactive
         if (nat.gt.100)
     &        call errquit('dft_gradients_so: 
     &                      dimension error in test',0, UNKNOWN_ERR)
         do i = 1, nat
            if (log_mb(k_act+i-1)) then
               g_rhs(1,i) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test x')
               g_rhs(2,i) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test y')
               g_rhs(3,i) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test z')
            endif
         enddo

         if (bgj_print() .gt. 0)
     &        write(6,*)'*** In dft_gradients_so: calling j_cpks_rhs'
         call j_cpks_rhs(iga_dens, log_mb(k_act), nactive, g_rhs)

         do i = 1, nat
            if (log_mb(k_act+i-1)) then
               do j = 1, 3
                  if (.not.ga_destroy(g_rhs(j,i))) then
                     call errquit('j_cpks_rhs: problem destroying ga',1,
     &       GA_ERR)
                  endif
               enddo
            endif
         enddo
      endif
c!!! BGJ test !!!
      status = ma_free_heap (l_act)
      status = ma_free_heap (l_frc_xc)
      status = ma_free_heap (l_frc_2el)
      status = ma_free_heap (l_force)
      status = ga_destroy (iga_dens(2))
      status = ga_destroy (iga_dens(1))
      status = ga_destroy (g_densso(2))
      status = ga_destroy (g_densso(1))
      status = ga_destroy (g_force)
c     
      return
      end
