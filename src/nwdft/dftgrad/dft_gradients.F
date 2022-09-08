      subroutine dft_gradients(rtdb)
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
*
* $Id$
*
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
#include "sym.fh"
#include "cdft.fh"
#include "util.fh"
#include "dftps.fh"
c
c!!! BGJ test !!!
#include "bgj.fh"
      integer l_hess, k_hess, g_rhs(3,nw_max_atom), j
c!!! BGJ test !!!
      double precision  zero, one, two
      Parameter (zero=0.d0, one=1.d0, two=2.d0)
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_header, movecs_read,xc_rep_close
      external movecs_read_header, movecs_read,xc_rep_close
c     integer noc(2)
      integer nmo(2)
      integer iga_dens(2), g_vecs(2), g_force
      integer idum(4), ndum
      double precision edum
      Integer k_evals(2), l_evals(2)
      double precision grad_norm, grad_max
      external grad_norm, grad_max
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      character*80 scftype
      integer ifocc
      logical status,frac_occ
      integer me, nproc, max_sh_bf, max_at_bf, nat, max_sh_bfcd,
     $     lforce, nactive, i, nbf_vecs, nsets, ispin, 
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
     $     nscr, lscr, iscr
      integer ipoint, itmpm,ltmpm,g_tmp(2),lenvec,maavail
      double precision charge, charge_nuc, rhffact, tol2e, onem,
     , toll
c
c     vdw
      double precision dum
      logical cgmin
      logical disp
      logical xc_chkdispauto
      external xc_chkdispauto
c
c     xdm
      integer xdmdisp
      integer nxdm
      integer ixdm_v, ixdm_ml, lxdm_v, lxdm_ml
      common /xdmd/ nxdm, ixdm_v, lxdm_v, ixdm_ml, lxdm_ml
c
      double precision fant_a,fant_d
      parameter(toll=1.d-9)
      logical  oprint_force_comps
c     
      logical has_frac_occ
      external has_frac_occ
c
      nproc = ga_nnodes()
      me=ga_nodeid()
      oprint_force_comps = util_print('force components', print_debug)
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
c     check for cgmin since it breaks cdfit
c
      if (.not.rtdb_get(rtdb,'dft:cgmin', mt_log, 1, cgmin))
     &   cgmin=.false.

c     if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
c    $     call errquit('dft_gradient: rtdb_get of noc failed', 0,
c    &       RTDB_ERR)
c     
c     check if fractional occupation is on
c
      frac_occ = .false.
      if (has_frac_occ(rtdb)) frac_occ = .true.
c     
c     allocate and initialize global and local memory
c     
c     mo-vectors
c     
      if (ipol .eq. 1)then
         g_vecs(1) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'mo vectors')
      else 
         g_vecs(1) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'alpha mo vectors')
         g_vecs(2) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'beta mo vectors')
      endif 
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
         iga_dens(1) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'alpha density matrix')
         iga_dens(2) = ga_create_atom_blocked (geom, ao_bas_han,
     $        'beta density matrix')
      endif 
c     
c     forces on atoms (3xnat)
c     
*ga:1:0
      status = ga_create(mt_dbl, 3, nat, 'forces', 3, 0, g_force)
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
      if (ipol .eq. 1)then
         if (.not. ma_alloc_get(mt_dbl, nbf_ao, 'MO evals', l_evals(1),
     $        k_evals(1)))
     $        call errquit('dft_gradient: could not allocate l_evals',1,
     &       MA_ERR)
      else 
         status = ma_alloc_get(mt_dbl, nbf_ao, 'alpha MO evals',
     $        l_evals(1), k_evals(1))
         status = status .and.
     $        ma_alloc_get(mt_dbl, nbf_ao, 'beta MO evals',
     $        l_evals(2), k_evals(2))
         if (.not. status)then
            call errquit('dft_gradient: could not allocate l_evals',1,
     &       MA_ERR)
         endif 
      endif
c     
c     occupation numbers (not used, but necessary for movecs_read)
c     
c     should do k_occ for both spins, in case used at some point...
c     
      if (.not. ma_alloc_get(mt_dbl, nbf_ao*ipol, 'occ. numbers',
     $     l_occ, k_occ))
     $     call errquit('dft_gradient: could not allocate l_occ',1,
     &       MA_ERR)
c     
c     lookup table and list of active atoms
c     
      if (.not. ma_alloc_get(MT_LOG, nat, 'active atoms',
     $     l_act, k_act))
     $     call errquit('grad: could not allocate l_act',1, MA_ERR)
      call grad_active_atoms(rtdb, nat, log_mb(k_act), nactive)
c     
c     get MO vectors from file
c     
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $     call errquit('dft_gradient: DFT MO vectors not defined',0,
     &       RTDB_ERR)
      status = movecs_read_header(movecs_in, title_vecs, basis_vecs,
     $         scftype_vecs, nbf_vecs, nsets, nmo, 2)
c
c     ipol  - number of spin channels   (RKS=1, ROKS=2, UKS=2)
c     nsets - number of sets of vectors (RKS=1, ROKS=1, UKS=2)
c
      if (.not. rtdb_cget(rtdb, 'dft:scftype', 1, scftype))
     $     call errquit('dft_gradient: DFT scftype not defined',0,
     &       RTDB_ERR)
      if (scftype.eq.'RHF') then
        if (ipol .ne. nsets .or. ipol.ne.1)then  
          write (LuOut,*) 'dft_gradient:  ERROR ipol, nsets:',ipol,nsets
          call errquit('dft_gradient:  ERROR ipol, nsets disagree',2,
     &                 INPUT_ERR)
        endif 
      elseif (scftype.eq.'ROHF') then
        if (nsets.ne.1.or.ipol.ne.2) then
          write (LuOut,*) 'dft_gradient:  ERROR ipol, nsets:',ipol,nsets
          call errquit('dft_gradient:  ERROR ipol, nsets disagree',2,
     &                 INPUT_ERR)
        endif 
      elseif (scftype.eq.'UHF') then
        if (nsets.ne.2.or.ipol.ne.2) then
          write (LuOut,*) 'dft_gradient:  ERROR ipol, nsets:',ipol,nsets
          call errquit('dft_gradient:  ERROR ipol, nsets disagree',2,
     &                 INPUT_ERR)
        endif
      else
        call errquit('dft_gradient: illegal scftype',0,UERR)
      endif
c     
c     Should check much more info than just nbf for consistency
c     
c     
c     get mo eigenvectors
c     
      if (nbf_ao .ne. nbf_vecs)then
         write(LuOut,*)'dft_gradient movecs output = ',movecs_in
         call errquit('dft_gradient: could not read mo vectors',911,
     &       DISK_ERR)
      else 
         status = .true.
         do ispin = 1, ipol
            status = status .and.
     $           movecs_read(movecs_in, min(ispin,nsets), 
     &                       dbl_mb(k_occ+(ispin-1)*nbf_ao),
     $                       dbl_mb(k_evals(ispin)), g_vecs(ispin))
         enddo 
      endif
c     
      if (.not.status)then
         write(LuOut,*)'dft_gradient movecs output = ',movecs_in
         call errquit('dft_gradient: could not read mo vectors',917,
     &       DISK_ERR)
      endif
c
      if(frac_occ) then
c
c       fractional occupation, therefore check new nocs
c
        if (.not. MA_Push_Get(MT_Dbl, nbf_ao, 'tmpm', ltmpm, itmpm))
     &     call errquit('dftgforce: failed to alloc tmpm',0, MA_ERR)
           rhffact = one
c
        do ispin=1,ipol
          g_tmp(ispin) = ga_create_atom_blocked(geom, ao_bas_han,
     &         'frac vecs')
          call ga_zero(g_tmp(ispin))
          ipoint=k_occ+(ispin-1)*nbf_ao-1
c
          do i = ga_nodeid()+1, nbf_ao, ga_nnodes()
            call get_col(g_vecs(ispin), nbf_ao, i, DBL_MB(itmpm))
            call dscal(nbf_ao, dbl_mb(i+ipoint), DBL_MB(itmpm), 1)
            call put_col(g_tmp(ispin), nbf_ao, i, DBL_MB(itmpm))
          enddo
          do i=1,nbf_ao
            if(dbl_mb(ipoint+i).ge.toll) noc(ispin)=i
          enddo
        enddo
        if (.not.ma_pop_stack(ltmpm))
     &     call errquit('dftg_force: cannot pop stack',0, MA_ERR)
      
      else
        do ispin=1,ipol
          g_tmp(ispin)=g_vecs(ispin)
        enddo
      endif
c
      do ispin = 1, ipol
c     
c        dens = vecs*vecs
c     
         if (odftps) call pstat_on(ps_dgemm)
         call ga_dgemm('n', 't', nbf_ao, nbf_ao, noc(ispin), rhffact, 
     $        g_tmp(ispin), g_vecs(ispin), 0.0d0, iga_dens(ispin))
         if (odftps) call pstat_off(ps_dgemm)
         call ga_symmetrize(iga_dens(ispin))
c     
c     free temporary arrays
c     
         if(frac_occ) then
           if(.not.ga_destroy (g_tmp(ispin))) call
     *          errquit('dftg_force: could not gadestr gtmp',ispin,
     &       GA_ERR)
         endif
      enddo   !ispin
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
      if (.not.ma_chop_stack(lcntoce))
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
      if (ipol .eq. 2)status = ga_destroy(g_vecs(2))
      status = ga_destroy(g_vecs(1))
c     
      if (.not.status)then
         call errquit('dft_gradient: could not destroy g_eigen_diag',1,
     &       GA_ERR)
      endif 
c     
      status = ma_free_heap(l_occ)
      if (ipol .eq. 2)then
         status = ma_free_heap (l_evals(2))
      endif 
      status = ma_free_heap (l_evals(1))
c
      if (CDFIT.and.(.not.cgmin)) then
c     
c     determine memory requirements for integral gradients
c     
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call int_mem_2e3c(max2e3c, mscratch_2e3c)
      lbuf = max(max1e, max2e)
      lbuf = max(lbuf, max2e3c) + 500
      lscratch = max(mscratch_1e, mscratch_2e)
      lscratch = max(lscratch, mscratch_2e3c)
c
c     fitting coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf_cd*3*ipol, 'CD coef',
     $     l_cdcoef, i_cdcoef))
     $     call errquit('dft_gradient: could not alloc CD coef',0,
     &     MA_ERR)
c     
c     one-electron contribution
c     buffers for one electron integral derivatives
c     
      status = ma_push_get(mt_dbl, lbuf, 'deriv buffer', l_buf, k_buf)
      if(.not.status) then
         maavail=MA_inquire_avail(mt_dbl)
         call errquit('dft_gradient: could not allocate buffer',
     &     8*(lbuf-maavail),  MA_ERR)
      endif
c     
      status = ma_push_get(mt_dbl, lscratch, 'deriv scratch', 
     $     l_scr, k_scr)
      if (.not.status) 
     $     call errquit('dft_gradient: could not allocate scratch',1,
     &       MA_ERR)
c     
c     allocate local density matrix block
c     we can eventually need up and down spins at the same time
c     
      lsqa = 2 * max_at_bf * max_at_bf
c     
      status = ma_push_get(mt_dbl, lsqa, 'local_w_density', 
     $     l_wdens, k_wdens)
      if (.not.status)call errquit('could not allocate l_dens',1,
     &       MA_ERR)
c     
c     store total DM in ga_dens(1) for standard CDFIT cases,
c     otherwise have the spin-resolved case
c     
      if ((.not.adft) .and. ipol .eq. 2)then
         call ga_dadd (one,iga_dens(1),one,iga_dens(2),iga_dens(1))
      endif 
c     
c     define threshold for Schwarz screening (same as in DFT)
c     
      tol2e=10.d0**(-itol2e)
c     
c     
c     compute 3 center coulomb derivative
c     
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c     
c     
         if(.not.bas_nbf_cn_max(cd_bas_han, max_sh_bfcd))
     $        call errquit('dftg_force: basnbfcdmax broken?',0,
     &       BASIS_ERR)
         
         lenvec=max(6*max_sh_bfcd,3*max_sh_bf*max_sh_bf)
         if (.not. MA_Push_get(MT_DBL, lenvec, 'svec',
     $        lsvec, isvec))
     $        call errquit('dftg_force: could not alloc svec',0, MA_ERR)
         ippp=k_wdens
c     
         if (odftps) call pstat_on(ps_vcoul)
         call xc_rep_init(rtdb, geom, ao_bas_han,iga_dens,iga_dens,
     &        nbf_ao,ipol,.true.,.true.)
         call dftg_cdfit(geom,ao_bas_han, cd_bas_han,
     $        nbf_cd, nat, tol2e, dbl_mb(k_scr), 
     $        lscratch, dbl_mb(k_buf), lbuf,
     $        dbl_mb(isvec), dbl_mb(ippp), max_sh_bf,
     $        iga_dens, dbl_mb(k_frc_2el), 
     $        DBL_MB(i_cdcoef), oskel, ipol, adft)
         if(.not.xc_rep_close(rtdb, nbf_ao,ipol,ipol,
     D        iga_dens,iga_dens,.true.)) call
     .        errquit(' dftggrad: xcrepclose failed ',0, 0)
c     
         call ga_dgop(msg_grad_2el, dbl_mb(k_frc_2el), lforce, '+')
         if (odftps) call pstat_off(ps_vcoul)
c     
         if (.not.ma_chop_stack(l_buf))
     $     call errquit('dft_gradient: cannot chop stack',0, MA_ERR)
c     
c        restore alpha DM in g_dens(1)
c     
         if ((.not.adft).and.ipol .eq. 2)then
            onem = -1.d0
            call ga_dadd(one, iga_dens(1), onem, iga_dens(2), 
     &                   iga_dens(1))
         endif 
      else
        i_cdcoef = 1
      endif  ! cdfit
c     
c     get exchange-correlation contribution to the gradient
c     
c     
c$$$      write(LuOut,*) ' BEFORE CALL TO GETXC'
c$$$      call ga_print(iga_dens(1))
c$$$      call ga_print(iga_dens(2))
c$$$      write(LuOut,*) ' nactive '
c$$$      call output(dbl_mb(irdens_atom), 1, nat, 1, nat, nat, nat, 1)
c$$$      write(LuOut,*) (int_mb(icntoce+i),i=0,nshells_ao-1)
c$$$      write(LuOut,*) (int_mb(icntobfr+i),i=0,2*nshells_ao-1)
c$$$      write(LuOut,*) (int_mb(icetobfr+i),i=0,2*nat-1)
      if (odftps) call pstat_on(ps_xc)
      call dftg_getxc(rtdb, nat, iga_dens, dbl_mb(k_frc_xc), 
     $     log_mb(k_act), nactive, 
     $     dbl_mb(irdens_atom),  int_mb(icetobfr), dbl_mb(i_cdcoef))

      if (cdfit.and..not.cgmin) then
         if (.not.ma_chop_stack(l_cdcoef))
     $     call errquit('dft_gradient: cannot chop stack',0, MA_ERR)
      endif

      call ga_dgop(msg_grad_xc,  dbl_mb(k_frc_xc),  lforce, '+')
      if(oprint_force_comps.and.me.eq.0)then
            write(luout,2200) 
     $           'XC gradient',
     $           ((dbl_mb(k_frc_xc+i-1+3*(j-1)),i=1,3),j=1,nat)
            write(luout,2200) 
     $            'CD gradient',
     $           ((dbl_mb(k_frc_2el+i-1+3*(j-1)),i=1,3),j=1,nat)
      endif
 2200       format(A/,1000(3(1x,F12.6),/))
      
      if (odftps) call pstat_off(ps_xc)
c     
c     vdW bit
c
      if (.not.rtdb_get(rtdb, 'dft:disp', mt_log, 1, disp))
     &   disp=.false.
      if(disp.or.xc_chkdispauto())
     &   call xc_vdw(rtdb,geom,dum,dbl_mb(k_frc_xc), 'forces')
c     
c     xdm bit
c
      if (.not.rtdb_get(rtdb, 'dft:xdm', mt_int, 1, xdmdisp))
     &   xdmdisp=0
      if (xdmdisp.ne.0) then
         call xc_xdm(rtdb,iga_dens,idum,nat,ndum,edum,
     &        dbl_mb(k_frc_xc),dbl_mb(ixdm_v),
     &        dbl_mb(ixdm_ml),'forces')
       if (.not. rtdb_put(rtdb,'dft:xdmsave', mt_log, 1, .false.)) 
     $     call errquit('dft_gradient: cannot rtdb_put',0, RTDB_ERR)
          call xc_xdm_cleanup(rtdb)
      endif
c
c     add Bonacic-Fantucci repulsive term
c
      if (.not.rtdb_get(rtdb, 'dft:fant_d', mt_dbl, 1,
     &   fant_d)) fant_d=-1d0
      if (.not.rtdb_get(rtdb, 'dft:fant_a', mt_dbl, 1,
     &   fant_a)) fant_a=-1d0
      if(fant_a.ne.-1d0.and.fant_d.ne.-1d0)
     A     call dftg_fant(geom,nat,fant_a,fant_d,
     A     dbl_mb(k_frc_xc))
c     
      if (ga_nodeid() .eq. 0)then
         status = rtdb_parallel (.false.)
         do i = 0, lforce-1
          dbl_mb(k_force+i) = dbl_mb(k_frc_2el+i) + dbl_mb(k_frc_xc+i) 
         enddo
         if (.not. rtdb_put(rtdb, 'dft:cd+xc gradient', mt_dbl,
     $        lforce, dbl_mb(k_force))) call errquit
     $        ('dft_gradient: failed storing cd+xc gradient',0,
     &       UNKNOWN_ERR)
c     
         status = rtdb_parallel (.true.)
      endif
      if(oprint_force_comps.and.me.eq.0)then
         do 31 i = 1, nat
            write (luout,2000) i, 
     &                         (dbl_mb(k_frc_xc+3*(i-1)+j),j=0,2)
 2000       format (1X,I3,3(1X,F10.6))
   31    continue
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
      if (ipol .eq. 2)then
         call ga_dadd (1d0,iga_dens(1),1d0,iga_dens(2),iga_dens(1))
      endif 
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
     &        call errquit('dft_gradients: could not alloc j hessian',
     &        1, MA_ERR)
         call dfill(9*nat*nat, 0.0d0, dbl_mb(k_hess), 1)
         if (bgj_print() .gt. 0)
     &        write(LuOut,*)'*** In dft_gradients: calling j_hessian'
         call j_hessian(iga_dens, log_mb(k_act), nactive,
     &        dbl_mb(k_hess))
         status = MA_pop_stack(l_hess)
         if (.not.status) call
     &        errquit('dft_gradients: could not pop j hessian',
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
c        Allocate and initialize temp GA's for RHS
c
         if (bgj_print() .gt. 0)
     &        write(*,*)'*** j cpks rhs test: nactive =',nactive
         if (nat.gt.100)
     &        call errquit('dft_gradients: dimension error in test',0,
     &       UNKNOWN_ERR)
         do i = 1, nat
            if (log_mb(k_act+i-1)) then
               g_rhs(1,i) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test x')
               g_rhs(2,i) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test y')
               g_rhs(3,i) = ga_create_atom_blocked
     &              (geom, ao_bas_han, 'CPKS RHS test z')
               call ga_zero(g_rhs(1,i))
               call ga_zero(g_rhs(2,i))
               call ga_zero(g_rhs(3,i))
            endif
         enddo

         if (bgj_print() .gt. 0)
     &        write(LuOut,*)'*** In dft_gradients: calling j_cpks_rhs'
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
c
      status = ma_free_heap (l_act)
      status = ma_free_heap (l_frc_xc)
      status = ma_free_heap (l_frc_2el)
      status = ma_free_heap (l_force)
      if (ipol .eq. 2) status = ga_destroy (iga_dens(2))
      status = ga_destroy (iga_dens(1))
      status = ga_destroy (g_force)
c     
      return
      end
c
      subroutine dftg_fant(geom,natoms,a,d,forces)
C     Bonacic-Kouteck, V.; Cespiva, L.; Fantucci, P.; Pittner, J.; Kouteck, J. J. 
C     Chem. Phys. 1994, 98, 490.
C     Mitric, R.; Hartmann, M.; Stanca, B.; Bonacic-Koutecky, V.; Fantucci, P.;
C     J. Phys. Chem. A.; (Article); 2001; 105(39); 8892-8905
C     The core-core repulsion has been corrected according to
C     (CC(rij) ) l/rij + D exp(-arij)). 
C     Constants D and a obtained for 1e-RECP's from the fitting procedure
C     described in Appendix A necessery for corection of core-core potential: 
C     DAg-Ag ) 1619.887392, aAg-Ag ) 2.49630105, 
C     DAu-Au ) 1911.131090, aAu-Au ) 2.46590129, 
C     DAu-Ag ) 1765.509532, aAu-Ag ) 2.48110117. 
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
      integer geom,natoms
      double precision a,d,forces(3,*)
c
      integer i,j,B
      double precision qi,qj
      double precision xB(3),xj(3),rBj,ffant,drbjdxb
      character*16 tag
c     
      do B=1,natoms
         if (.not. geom_cent_get(geom, B, tag, 
     &        xB, qi))call errquit
     &        ('grid_acc_def: geom_cent_get failed', 0, GEOM_ERR)
         do j=1,natoms
            if(B.ne.j) then
               if (.not. geom_cent_get(geom, j, tag, 
     &              xj, qj))call errquit
     &              ('grid_acc_def: geom_cent_get failed', 0, GEOM_ERR)
               rBj=sqrt((xB(1)-xj(1))**2+(xB(2)-xj(2))**2+
     +              (xB(3)-xj(3))**2)
                  ffant=d*exp(-a*rBj)*(-a)
               do i=1,3
                  drBjdxb=1d0/rBj*(xB(i)-xj(i))
                  forces(i,B)=forces(i,B) + ffant*drBjdxB
               enddo
            endif

         enddo
      enddo
!      if(ga_nodeid().eq.0)
!     W     write(luout,*) 
!     W     ' Bonacic-Koutecky Fantucci-Repulsive Term ',dft_fant
      return
      end
