      logical function int_ecp_init(ecpidin,soidin,indx_grad)
* $Id$
      implicit none
*::cr::7
*--------------------------------------------------*
* COPYRIGHT (C) 1994, 1995, 1996, 1997, 1998, 1999 *
*         Pacific Northwest National Laboratory,   * 
*         Battelle Memorial Institute.             *
*--------------------------------------------------*
*------------> All Rights Reserved <---------------*
*--------------------------------------------------*
*
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "geom.fh"
#include "geomP.fh"
#include "apiP.fh"
#include "ecp_nwc.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "sym.fh"
*::passed
      integer ecpidin   ! [input] ecp basis handle
      integer soidin    ! [input] so potential handle
      integer indx_grad ! [input] gradient level 0=energy 1=grad 2=hess
*::local
c
      integer ecpid                 ! lexical index for ecp basis handle
      integer soid                  ! lexical index for ecp basis handle
      integer geom                  ! geometry handle
      integer geomecp               ! geometry handle from ecp
      integer geomso                ! geometry handle from so pot.
*      integer num_ecp               ! number of ecp centers (counted)
      integer ncenters              ! number of centers
      integer n_zeta_c_e            ! length of ecp exponent array
      integer n_zeta_c_c            ! length of ecp coefficient array
      integer nz_add                ! increment for exp/coef counters
*      integer icent                 ! counter for center loop
*      integer iucent                ! unique center of icent
      integer l_ecp_sz              ! size of ang info for ecp pointer arrays
      integer l_so                  ! maximum angular momentum of so pot.
      integer lmax_both
      integer idum
      integer basis_handle
      integer l_bas
      integer bas_id
      integer basis_ncont
      integer i_cont, ic, iuc, ie, icf, inp, ing, il, ig
      integer j_cont, jc, juc, je, jcf, jnp, jng, jl, jg
      integer lecp_mem, mem_ecp, lecp_mem1
*      integer h_tlist, k_tlist
      integer nV
      logical o_both, o_ecp, o_so
      double precision dpdum(2),qfact
      logical oskel
      integer ictr,ic1,ic2
      integer jctr,jc1,jc2
*
      external ecp_init_bd
*
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
*<<<------------------------------------------------------------>
*<<<-First pass is to explode all information into arrays
*<<<-The memory used can be reduced by keeping track of the
*<<< tag to basis-unique tag information so that the pointers
*<<< for exponent and coefficient arrays is similar to the 
*<<< storage in the basis set object.  It may be possible to 
*<<< pass the entire exndcf array for the ecp basis to the 
*<<< integral routines with the proper pointer mechanism in 
*<<< place.
*<<<
*<<< RAK Apr 1996
*<<<------------------------------------------------------------>
c
c
c.. check initialization
      if (init_ecp_init) then
        write(6,*)' already called int_ecp_init' 
        call errquit('int_ecp_init error',911, INT_ERR)
      endif
      init_ecp_init = .true.
c.. determine if ecp/so or both are needed
      o_ecp = ecpidin.ne.0
      o_so  = soidin.ne.0
      o_both = o_ecp.and.o_so
c.. lexical basis index
      ecpid = -1
      soid = -1
      if (o_ecp) ecpid = ecpidin + BASIS_HANDLE_OFFSET
      if (o_so)  soid  = soidin  + BASIS_HANDLE_OFFSET

c.. determine geom
      if (o_ecp) then
        geomecp = ibs_geom(ecpid)
        geom = geomecp
      endif
      if (o_so) then
        geomso  = ibs_geom(soid)
        geom = geomso
      endif
      if (o_both) then
        if (geomso.ne.geomecp) then
          write(luout,*)' geometry used to load ecp ',geomecp
          write(luout,*)' geometry used to load so  ',geomso
          call errquit('int_ecp_init: ecp/so geometry mismatch',911,
     &       GEOM_ERR)
        endif
        geom = geomecp
      endif
*
*      build up ecp/so information for heap
c....
c.. get total number of centers
      if (.not.geom_ncent(geom,ncenters)) call errquit
     &    ('int_ecp_init: geom_ncent failed?',911, GEOM_ERR)
c.. get number of ecp/so/both centers
      n_ecp = 0
      if (o_both) then
        if (.not.ecpso_ncent
     &      (geom,soidin,ecpidin,n_ecp)) call errquit
     &      ('int_ecp_init: ecpso_ncent failed',911, INT_ERR)
      else if (o_ecp) then
        if (.not.ecp_ncent
     &      (geom,ecpidin,n_ecp)) call errquit
     &      ('int_ecp_init: ecp_ncent failed?',911, INT_ERR)
      else if (o_so) then
        if (.not.so_ncent
     &      (geom,soidin,n_ecp)) call errquit
     &      ('int_ecp_init: so_ncent failed',911, INT_ERR)
      else
        call errquit('int_ecp_init: no so or ecp basis',911,
     &       BASIS_ERR)
      endif
      if (n_ecp.eq.0) call errquit
     &    ('int_ecp_init: fatal error no ecp/so centers',911,
     &       BASIS_ERR)

c.. allocate space for coordintates of ecp/so centers
      if (.not.ma_alloc_get(mt_dbl,
     &    (3*n_ecp),
     &    'ecp center coords',
     &    h_xyzecp, k_xyzecp)) call errquit
     &    (' int_ecp_init: ecp center coords ma failed ',911,
     &       BASIS_ERR)
      call dfill((3*n_ecp),0.0d00,dbl_mb(k_xyzecp),1)

c.. get coordinates for ecp centers.
      if (o_both) then
        if (.not.ecpso_get_coords(geom,soidin,ecpidin,
     &      n_ecp,dbl_mb(k_xyzecp)))
     &    call errquit('int_ecp_init: ecpso_get_coords failed',911,
     &       INT_ERR)
      else if (o_ecp) then
        if (.not.ecp_get_coords(geom,ecpidin,
     &      n_ecp,dbl_mb(k_xyzecp)))
     &    call errquit('int_ecp_init: ecp_get_coords failed',911,
     &       INT_ERR)
      else if (o_so) then
        if (.not.so_get_coords(geom,soidin,
     &      n_ecp,dbl_mb(k_xyzecp)))
     &    call errquit('int_ecp_init: so_get_coords failed',911,
     &       INT_ERR)
      else
        call errquit('int_ecp_init: no so or ecp basis',911,
     &       INT_ERR)
      endif
*      write(6,*)' coordinates after read '
*      call output(dbl_mb(k_xyzecp),1,3,1,n_ecp,3,n_ecp,1)
c
c...       now comes the tricky part
c
c allocate and fill an exponent and coeff array for the ecp basis
c
*      write(6,*)'inside ecp init'
*      if (.not.ecp_print(ecpidin)) stop ' ecp_print error'
      if (o_both) then
        if (.not.bas_num_prim(ecpidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_prim failed?',911, BASIS_ERR)
        n_zeta_c_e = nz_add
        if (.not.bas_num_prim(soidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_prim failed?',911,
     &       BASIS_ERR)
        n_zeta_c_e = n_zeta_c_e + nz_add
        if (.not.bas_num_coef(ecpidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_coef failed?',911, BASIS_ERR)
        n_zeta_c_c = nz_add
        if (.not.bas_num_coef(soidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_prim failed?',911, BASIS_ERR)
        n_zeta_c_c = n_zeta_c_c + nz_add
      else if (o_ecp) then
        if (.not.bas_num_prim(ecpidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_prim failed?',911, BASIS_ERR)
        n_zeta_c_e = nz_add
        if (.not.bas_num_coef(ecpidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_coef failed?',911, BASIS_ERR)
        n_zeta_c_c = nz_add
      else if (o_so) then
        if (.not.bas_num_prim(soidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_prim failed?',911, BASIS_ERR)
        n_zeta_c_e = nz_add
        if (.not.bas_num_coef(soidin,nz_add)) call errquit
     &      ('int_ecp_init:bas_num_coef failed?',911, BASIS_ERR)
        n_zeta_c_c = nz_add
      else
        call errquit('int_ecp_init: no so or ecp basis',911, BASIS_ERR)
      endif
*.. error check for general contraction
      if (n_zeta_c_e.ne.n_zeta_c_c) then
        write(6,*)' possible general contraction on ecp/so basis'
        call errquit ('int_ecp_init: n_zeta_c_e .ne. n_zeta_c_c',
     &      (n_zeta_c_e- n_zeta_c_c), BASIS_ERR)
      endif
      n_zeta_c = n_zeta_c_e ! use e one
c
cc AJL/Begin/SPIN ECPs
c
c... determine the number of channels needed (Spin-polarised ECPs)
      ecp_channels = 1 
      if (o_ecp) then
        if (.not.ecp_get_high_chan(ecpidin,ecp_channels))
     &    call errquit('int_ecp_init error',911, INT_ERR)
      end if
cc AJL/Debug
c      write(6,*) 'Number of ECP Channels: ', ecp_channels
c      write(6,*) '(1 = Spin Independent, Default; 2 = Spin Polarised)'
cc AJL/End
c
c.. allocate space for ecp/so exponents and ecp/so coefficients
      if (.not.ma_alloc_get(mt_dbl,
     &    n_zeta_c,
     &    'ecp exponents',
     &    h_ecp_e, k_ecp_e)) call errquit
     &    (' int_ecp_init: ecp exponent ma failed ',911, BASIS_ERR)
      call dfill(n_zeta_c,0.0d00,dbl_mb(k_ecp_e),1)
      if (.not.ma_alloc_get(mt_dbl,
     &    n_zeta_c,
     &    'ecp coefficients',
     &    h_ecp_c, k_ecp_c)) call errquit
     &    (' int_ecp_init: ecp coefficients ma failed ',911, BASIS_ERR)
      call dfill(n_zeta_c,0.0d00,dbl_mb(k_ecp_c),1)
c
cc AJL/Begin/SPIN ECPs
      if (ecp_channels.gt.1) then
        if (.not.ma_alloc_get(mt_dbl,
     &      n_zeta_c,
     &      'beta ecp exponents',
     &      h_ecp_e_beta, k_ecp_e_beta)) call errquit
     &      (' int_ecp_init: beta ecp exponent ma failed ',
     &       911, BASIS_ERR)
        call dfill(n_zeta_c,0.0d00,dbl_mb(k_ecp_e_beta),1)
        if (.not.ma_alloc_get(mt_dbl,
     &      n_zeta_c,
     &      'beta ecp coefficients',
     &      h_ecp_c_beta, k_ecp_c_beta)) call errquit
     &      (' int_ecp_init: beta ecp coefficients ma failed ',
     &      911, BASIS_ERR)
        call dfill(n_zeta_c,0.0d00,dbl_mb(k_ecp_c_beta),1)      
      end if
cc AJL/End
c
c... determine maximum angular momentum of ecp basis
      l_ecp = 0
      if (o_ecp) then
        if (.not.bas_high_angular(ecpidin,l_ecp)) call errquit
     &      ('int_ecp_init:ecp: bas_high_angular failed',911, BASIS_ERR)
      endif
      if (o_so) then
        if (.not.bas_high_angular(soidin,l_so)) call errquit
     &      ('int_ecp_init:so: bas_high_angular failed',911, BASIS_ERR)
        l_ecp = max(l_ecp,l_so)
      endif
c
      l_ecp_sz = l_ecp + 2       !    (-1->Lval)
c
c... allocate space for n_prim_C(0:4,-1:l_ecp_max,n_C*2),
      if (.not.ma_alloc_get(mt_int,
     &    (5*l_ecp_sz*n_ecp*2),
     &    'ecp n_prim_C',
     &    h_ecp_nprim_c, k_ecp_nprim_c)) call errquit
     &    (' int_ecp: ecp nprim_c ma failed ',911, BASIS_ERR)
      call ifill((5*l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_nprim_c),1)
c... allocate space for n_coef_C(-1:l_ecp_max,n_C*2)
      if (.not.ma_alloc_get(mt_int,
     &    (l_ecp_sz*n_ecp*2),
     &    'ecp n_coef_C',
     &    h_ecp_ncoef_c, k_ecp_ncoef_c)) call errquit
     &    (' int_ecp: ecp ncoef_c ma failed ',911, BASIS_ERR)
      call ifill((l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_ncoef_c),1)
c... allocate space for ind_C      
      if (.not.ma_alloc_get(mt_int,
     &    (l_ecp_sz*n_ecp*2),
     &    'ecp ind_C',
     &    h_ecp_ind_c, k_ecp_ind_c)) call errquit
     &      (' int_ecp: ecp ind_c ma failed ',911, BASIS_ERR)
      call ifill((l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_ind_c),1)
c... allocate space for ind_z      
      if (.not.ma_alloc_get(mt_int,
     &    (l_ecp_sz*n_ecp*2),
     &    'ecp ind_z',
     &    h_ecp_ind_z, k_ecp_ind_z)) call errquit
     &      (' int_ecp: ecp ind_z ma failed ',911, BASIS_ERR)
      call ifill((l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_ind_z),1)
c... allocate space for l_C      
      if (.not.ma_alloc_get(mt_int,
     &    n_ecp,
     &    'ecp l_C',
     &    h_ecp_l_c, k_ecp_l_c)) call errquit
     &    (' int_ecp: ecp l_c ma failed ',911, BASIS_ERR)
      call ifill(n_ecp,0,int_mb(k_ecp_l_c),1)
c... allocate space for ecp center pointer list
      if (.not.ma_alloc_get(mt_int,
     &    n_ecp,
     &    'ecp lexical indeces for ecp centers',
     &    h_ecp_lip, k_ecp_lip)) call errquit
     &    (' int_ecp: ecp_lip ma failed',911, BASIS_ERR)
      call ifill(n_ecp,0,int_mb(k_ecp_lip),1)
c
cc AJL/Begin/SPIN ECPs
c... allocate space for n_prim_C(0:4,-1:l_ecp_max,n_C*2),
      if (ecp_channels.gt.1) then
        if (.not.ma_alloc_get(mt_int,
     &      (5*l_ecp_sz*n_ecp*2),
     &      'beta ecp n_prim_C',
     &      h_ecp_nprim_c_beta, k_ecp_nprim_c_beta)) call errquit
     &      (' int_ecp: beta ecp nprim_c ma failed ',911, BASIS_ERR)
        call ifill((5*l_ecp_sz*n_ecp*2),
     &    0,int_mb(k_ecp_nprim_c_beta),1)
c... allocate space for n_coef_C(-1:l_ecp_max,n_C*2)
        if (.not.ma_alloc_get(mt_int,
     &      (l_ecp_sz*n_ecp*2),
     &      'beta ecp n_coef_C',
     &      h_ecp_ncoef_c_beta, k_ecp_ncoef_c_beta)) call errquit
     &      (' int_ecp: beta ecp ncoef_c ma failed ',911, BASIS_ERR)
        call ifill((l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_ncoef_c_beta),1)
c... allocate space for ind_C      
        if (.not.ma_alloc_get(mt_int,
     &      (l_ecp_sz*n_ecp*2),
     &      'beta ecp ind_C',
     &      h_ecp_ind_c_beta, k_ecp_ind_c_beta)) call errquit
     &      (' int_ecp: beta ecp ind_c ma failed ',911, BASIS_ERR)
        call ifill((l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_ind_c_beta),1)
c... allocate space for ind_z      
        if (.not.ma_alloc_get(mt_int,
     &      (l_ecp_sz*n_ecp*2),
     &      'beta ecp ind_z',
     &      h_ecp_ind_z_beta, k_ecp_ind_z_beta)) call errquit
     &      (' int_ecp: beta ecp ind_z ma failed ',911, BASIS_ERR)
        call ifill((l_ecp_sz*n_ecp*2),0,int_mb(k_ecp_ind_z_beta),1)
c... allocate space for l_C      
        if (.not.ma_alloc_get(mt_int,
     &      n_ecp,
     &      'beta ecp l_C',
     &      h_ecp_l_c_beta, k_ecp_l_c_beta)) call errquit
     &      (' int_ecp: beta ecp l_c ma failed ',911, BASIS_ERR)
        call ifill(n_ecp,0,int_mb(k_ecp_l_c_beta),1)
c... allocate space for ecp center pointer list
        if (.not.ma_alloc_get(mt_int,
     &      n_ecp,
     &      'beta ecp lexical indices for ecp centers',
     &      h_ecp_lip_beta, k_ecp_lip_beta)) call errquit
     &      (' int_ecp: beta ecp_lip ma failed',911, BASIS_ERR)
        call ifill(n_ecp,0,int_mb(k_ecp_lip_beta),1)
      endif
cc AJL/End
c
c
      call int_ecp_build_ecp_ptrs(ecpidin,soidin,geom,
     &    o_both,o_ecp,o_so,
     &    ncenters,
     &    n_ecp,
     &    l_ecp,
     &    n_zeta_c,
cc AJL/Begin/SPIN ECPs
     &    1, ! Channel to be filled. Would be 2 for beta channel
     &    ecp_channels, ! Total no. of channels, for print_ptrs
cc AJL/End
     &    int_mb(k_ecp_nprim_c),
     &    int_mb(k_ecp_ncoef_c),
     &    int_mb(k_ecp_ind_c),
     &    int_mb(k_ecp_ind_z),
     &    int_mb(k_ecp_l_c),
     &    int_mb(k_ecp_lip),
     &    dbl_mb(k_ecp_e),
     &    dbl_mb(k_ecp_c) )
c
cc AJL/Begin/SPIN ECPs
      if (ecp_channels.gt.1) ! Setup beta options
     &    call int_ecp_build_ecp_ptrs(ecpidin,soidin,geom,
     &    o_both,o_ecp,o_so,
     &    ncenters,
     &    n_ecp,
     &    l_ecp,
     &    n_zeta_c,
     &    2, ! Beta channel
     &    ecp_channels, ! Total no. of channels, for print_ptrs
     &    int_mb(k_ecp_nprim_c_beta),
     &    int_mb(k_ecp_ncoef_c_beta),
     &    int_mb(k_ecp_ind_c_beta),
     &    int_mb(k_ecp_ind_z_beta),
     &    int_mb(k_ecp_l_c_beta),
     &    int_mb(k_ecp_lip_beta),
     &    dbl_mb(k_ecp_e_beta),
     &    dbl_mb(k_ecp_c_beta))
cc AJL/End
c
*...   allocate space for c2s and s2c internal ecp transformation routines
c
c determine lmax among ao basis and ecp basis
c l_ecp currently has Lmax for ecp basis
      if (.not.ecp_get_parent_handle(ecpidin,basis_handle))
     &      call errquit
     &      ('int_ecp_init: ecp_get_parent_handle failed',911,
     &       BASIS_ERR)
      if (.not.bas_high_angular(basis_handle,l_bas)) call errquit
     &      ('int_ecp_init: bas_high_angular failed for ao handle',
     &      911, BASIS_ERR)
      lmax_both = max(l_ecp,l_bas) + l_bas + 2 
      call ecp_init_c2s(lmax_both,dpdum,dpdum,idum,1,1,.true.,mem_c2s)
      if (.not.ma_alloc_get(mt_dbl,
     &      mem_c2s,
     &      'ecp c2s routines',
     &      h_ecp_c2s, k_ecp_c2s)) call errquit
     &      ('int_ecp_init: ma failed for c2s',911, MA_ERR)
      call dfill(mem_c2s,0.0d00,dbl_mb(k_ecp_c2s),1)
      call ecp_init_c2s(lmax_both,
     &      dbl_mb(k_ecp_c2s),dbl_mb(k_ecp_c2s),mem_c2s,1,1,.false.,
     &      idum)
c
c initialize constants for ecp integral code
c
      call ecp_init_con()
c
c determine maximum memory for ecp integrals
c
      bas_id = basis_handle + basis_handle_offset
      basis_ncont = ncont_tot_gb(bas_id)
      ig = ibs_geom(bas_id)
      jg = ig
      mem_ecp = 0
*      if (.not.bas_print(basis_handle)) stop ' error'
*      if (.not.gbs_map_print(basis_handle)) stop ' error '
*      if (.not.bas_print(ecpidin)) stop 'error'
*      if (.not.gbs_map_print(ecpidin)) stop ' error'
      oskel=.true.
      do ictr=1,ncenters
         if (sym_atom(ig, ictr, qfact))  then
            if(.not.(bas_ce2cnr(basis_handle,ictr,ic1,ic2)))
     .           call errquit('intecp:basce2cnr failed',0, BASIS_ERR)
            do i_cont = ic1,ic2
               ic  = sf_ibs_cn2ce(i_cont,bas_id)
               iuc = sf_ibs_cn2ucn(i_cont,bas_id)
               ie  = infbs_cont(cont_iexp,iuc,bas_id)
               icf = infbs_cont(cont_icfp,iuc,bas_id)
               inp = infbs_cont(cont_nprim,iuc,bas_id)
               ing = infbs_cont(cont_ngen,iuc,bas_id)
               il  = infbs_cont(cont_type,iuc,bas_id) 
               do jctr=ictr,ncenters
                  if (sym_atom(ig, jctr, qfact))  then
                     if(.not.(bas_ce2cnr(basis_handle,jctr,jc1,jc2)))
     .                    call errquit('intecp:basce2cnr failed',0,
     &       BASIS_ERR)
                     do j_cont = jc1,jc2
                        jc = sf_ibs_cn2ce(j_cont,bas_id)
                        juc = sf_ibs_cn2ucn(j_cont,bas_id)
                        je  = infbs_cont(cont_iexp,juc,bas_id)
                        jcf = infbs_cont(cont_icfp,juc,bas_id)
                        jnp = infbs_cont(cont_nprim,juc,bas_id)
                        jng = infbs_cont(cont_ngen,juc,bas_id)
                        jl  = infbs_cont(cont_type,juc,bas_id)
                        nV  = int_nbf_x(il)*int_nbf_x(jl)*ing*jng
                        lecp_mem = 90000
                        if (indx_grad.eq.0) then
                           if (o_ecp) then
                             call int_ecp_hf1(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,dpdum,lecp_mem1,.true.)
                             lecp_mem = max(lecp_mem,lecp_mem1)
                           endif
                           if (o_so) then
                             call intso_hf1(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,dpdum,lecp_mem1,.true.)
                             lecp_mem = max(lecp_mem,lecp_mem1)
                           endif
                        elseif (indx_grad.eq.1) then
                           if (o_ecp) then
                             call intd_ecp_hf1(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,ncenters,dpdum,lecp_mem1,
     &                          .true.)
                             lecp_mem = max(lecp_mem,lecp_mem1)
                           endif
                           if (o_so) then
                             call intd_so_hf1(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,ncenters,dpdum,lecp_mem1,
     &                          .true.)
                             lecp_mem = max(lecp_mem,lecp_mem1)
                           endif
                        elseif (indx_grad.eq.2) then
                           if (o_ecp) then
                             call intdd_ecp_hf1(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,ncenters,dpdum,lecp_mem1,
     &                          .true.)
                             lecp_mem = max(lecp_mem,lecp_mem1)
                           endif
                           if (o_so) then
                             call intdd_so_hf1(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,ncenters,dpdum,lecp_mem1,
     &                          .true.)
                             lecp_mem = max(lecp_mem,lecp_mem1)
                           endif
                        else
                           call errquit(
     .                          'int_ecp_init:unknown initializ',911,
     &       INT_ERR)
                        endif
c     ! nint block used in api routine int_hf1sp                  
                        lecp_mem = lecp_mem + nV 
c
cc AJL/Begin/SPIN ECPs
cc Need to account for scenario where we have two channels so
cc double the number of ECPs could be generated
                        if (o_ecp.and.ecp_channels.gt.1) then
                          if (indx_grad.eq.0) then
                            call int_ecp_hf1_beta(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,dpdum,lecp_mem,.true.)
                          elseif (indx_grad.eq.1) then
                            call intd_ecp_hf1_beta(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,ncenters,dpdum,lecp_mem,.true.)
                          elseif (indx_grad.eq.2) then
                            call intdd_ecp_hf1_beta(
     &                          coords(1,ic,ig),
     &                          dbl_mb(mb_exndcf(ie,bas_id)),
     &                          dbl_mb(mb_exndcf(icf,bas_id)),
     &                          inp,ing,il,ic,
     &                          coords(1,jc,jg),
     &                          dbl_mb(mb_exndcf(je,bas_id)),
     &                          dbl_mb(mb_exndcf(jcf,bas_id)),
     &                          jnp,jng,jl,jc,
     &                          dpdum,nV,ncenters,dpdum,lecp_mem,.true.)
                          endif
c Add data storage space to lecp_mem
                          lecp_mem = lecp_mem + nV
                        endif
cc AJL/End
c
*                        write(6,*)' i_cont ',i_cont
*                        write(6,*)' j_cont ',j_cont
                        mem_ecp = max(mem_ecp,lecp_mem)
*                        write(6,*)' lecp_mem from int_ecp_init is ',
*     &        lecp_mem,mem_ecp
                     enddo
                  endif
               enddo
            enddo
         endif
      enddo
*     write(6,*)' scr 1e memory without ecp',mem_1e
      mem_1e = max(mem_1e,mem_ecp)
*     write(6,*)' scr 1e memory with    ecp',mem_1e
*     
      call int_nbf_max(basis_handle,mem_ecp)
c
      if (indx_grad.eq.0) then
         mem_ecp=mem_ecp*mem_ecp
      elseif (indx_grad.eq.1) then
         mem_ecp=mem_ecp*mem_ecp*3*ncenters
      else
         mem_ecp=mem_ecp*mem_ecp*3*3*(ncenters*(ncenters-1)/2+ncenters)
         mem_1e = max(mem_1e,mem_ecp)
      endif
      isz_1e = max (isz_1e,mem_ecp)
c
      int_ecp_init = .true.
*     
      end
      subroutine int_ecp_terminate()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "ecp_nwc.fh"

      logical status

      status = .true.
      status = status .and. MA_free_heap(h_xyzecp)
      status = status .and. MA_free_heap(h_ecp_e)
      status = status .and. MA_free_heap(h_ecp_c)
      status = status .and. MA_free_heap(h_ecp_nprim_c)
      status = status .and. MA_free_heap(h_ecp_ncoef_c)
      status = status .and. MA_free_heap(h_ecp_ind_c)
      status = status .and. MA_free_heap(h_ecp_ind_z)
      status = status .and. MA_free_heap(h_ecp_l_c)
      status = status .and. MA_free_heap(h_ecp_c2s)
      status = status .and. MA_free_heap(h_ecp_lip)
      h_xyzecp      = 0
      k_xyzecp      = 0
      h_ecp_e       = 0
      k_ecp_e       = 0
      h_ecp_c       = 0
      k_ecp_c       = 0
      h_ecp_nprim_c = 0
      k_ecp_nprim_c = 0
      h_ecp_ncoef_c = 0
      k_ecp_ncoef_c = 0
      h_ecp_ind_c   = 0
      k_ecp_ind_c   = 0
      h_ecp_ind_z   = 0
      k_ecp_ind_z   = 0
      h_ecp_l_c     = 0
      k_ecp_l_c     = 0
      h_ecp_c2s     = 0
      k_ecp_c2s     = 0
      h_ecp_lip     = 0
      k_ecp_lip     = 0
      n_zeta_c      = 0
      l_ecp         = 0
      n_ecp         = 0
      init_ecp_init = .false.
c
cc AJL/Begin
      if (ecp_channels.gt.1) then
        status = status .and. MA_free_heap(h_ecp_e_beta)
        status = status .and. MA_free_heap(h_ecp_c_beta)
        status = status .and. MA_free_heap(h_ecp_nprim_c_beta)
        status = status .and. MA_free_heap(h_ecp_ncoef_c_beta)
        status = status .and. MA_free_heap(h_ecp_ind_c_beta)
        status = status .and. MA_free_heap(h_ecp_ind_z_beta)
        status = status .and. MA_free_heap(h_ecp_l_c_beta)
        status = status .and. MA_free_heap(h_ecp_lip_beta)

        h_ecp_e_beta       = 0
        k_ecp_e_beta       = 0
        h_ecp_c_beta       = 0
        k_ecp_c_beta       = 0
        h_ecp_nprim_c_beta = 0
        k_ecp_nprim_c_beta = 0
        h_ecp_ncoef_c_beta = 0
        k_ecp_ncoef_c_beta = 0
        h_ecp_ind_c_beta   = 0
        k_ecp_ind_c_beta   = 0
        h_ecp_ind_z_beta   = 0
        k_ecp_ind_z_beta   = 0
        h_ecp_l_c_beta     = 0
        k_ecp_l_c_beta     = 0
        h_ecp_lip_beta     = 0
        k_ecp_lip_beta     = 0
        ecp_channels       = 0 ! Clear out ECP channels and we are done
      end if
cc AJL/End
c
      if (status) return
      call errquit
     &    ('int_ecp_terminate: error freeing heap',911, MEM_ERR)
      end
      subroutine int_ecp_build_ecp_ptrs(ecpidin,soidin,geom,
     &    o_both, o_ecp, o_so,
     &    ncenters,
     &    n_ecp,
     &    l_ecp,
     &    nz_ecp,
cc AJL/Begin/SPIN ECPs
     &    channel, ! spin channel
     &    channels, ! total number of channels
cc AJL/End
     &    n_prim_C,
     &    n_coef_C,
     &    ind_C,
     &    ind_Z,
     &    l_C,
     &    i_cent_C,
     &    c_exp,
     &    c_coef)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      integer ecpidin  ! [input] ecp basis set handle
      integer soidin   ! [input] so pot. handle
      logical o_both, o_ecp, o_so ![input] which ones to set up
      integer ncenters ! [input] number of centers
      integer n_ecp    ! [input] number of ecp centers
      integer l_ecp    ! [input] maximum angular momentum in ecp basis
      integer nz_ecp   ! [input] number of prims/coeffs in stored data structure. 
cc AJL/Begin/SPIN ECPs
      integer channel  ! [input] Current spin channel used for ECPs
      integer channels ! [input] Total channels, used for print
cc AJL/End
      integer n_prim_C(0:4,-1:l_ecp,n_ecp,2) ! [output] 
      integer n_coef_C(-1:l_ecp,n_ecp,2)     ! [output] 
      integer ind_C(-1:l_ecp,n_ecp,2)        ! [output] 
      integer ind_z(-1:l_ecp,n_ecp,2)        ! [output] 
      integer l_C(n_ecp)                      ! [output] 
      integer i_cent_C(n_ecp)                 ! [output]
      double precision c_exp(nz_ecp)          ! [output] 
      double precision c_coef(nz_ecp)         ! [output] 
c
      logical on_ecp, on_so
      integer geom        ! geometry handle
      integer ecpid       ! lexical basis set handle
      integer soid        ! lexical so pot handle
      integer icent       ! counter for centers
      integer iucent      ! unique map of icent
      integer p_nprim     ! counter/pointer into exp/coeff array
      integer num_ecp     ! ecp center as counted
      integer F_cont      ! first contraction on center iucent
      integer L_cont      ! last  contraction on center iucent
      integer iucont      ! contraction counter
      integer type        ! function type (-1 = local, 0-lval is non-local)
      integer nprim       ! number of prims in a contraction
      integer ncoef       ! number of coefficients in a contraction
      integer iexp        ! pointer into exndcf for exponents
      integer icfp        ! pointer into exndcf for coefficients
      integer irexp       ! pointer into exndcf for r-exponents
      integer n0,n1,n2    ! r-exponent count 
      integer n3,n4       ! r-exponent count 
      integer ie          ! ecp/so counter
cc AJL/Begin/SPIN ECPs
      integer mychan      ! channel of contracted ECP
cc AJL/End
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
*      write(6,*)' inside build'
      ecpid = ecpidin + BASIS_HANDLE_OFFSET
      soid  = soidin  + BASIS_HANDLE_OFFSET
c... zero arrays
      call ifill((5*(l_ecp+2)*n_ecp*2),0,n_prim_C,1)
      call ifill(((l_ecp+2)*n_ecp*2),0,n_coef_C,1)
      call ifill(((l_ecp+2)*n_ecp*2),0,ind_C,1)
      call ifill(((l_ecp+2)*n_ecp*2),0,ind_z,1)
      call ifill(n_ecp,0,l_C,1)
      call ifill(n_ecp,0,i_cent_C,1)
      call dfill(nz_ecp,0.0d00,c_exp,1)
      call dfill(nz_ecp,0.0d00,c_coef,1)
c
      if (o_both) then
        if (.not.
     &      ecpso_list_ncent(geom,soidin,ecpidin,num_ecp,i_cent_C))
     &      call errquit('ecp_ptrs: ecpso_list failed',911, BASIS_ERR)
      else if (o_ecp) then
        if (.not.ecp_list_ncent(geom,ecpidin,num_ecp,i_cent_C))
     &      call errquit('ecp_ptrs: ecp_list failed',911, BASIS_ERR)
      else if (o_so) then
        if (.not.so_list_ncent(geom,soidin,num_ecp,i_cent_C))
     &      call errquit('ecp_ptrs: ecpso_list failed',911, BASIS_ERR)
      else 
        call errquit('ecp_ptrs: should not get here',911, BASIS_ERR)
      endif
      if (num_ecp.ne.n_ecp) call errquit
     &    ('ecp_ptrs: mismatch in num_ecp and n_ecp',
     &    (num_ecp-n_ecp), BASIS_ERR)
c
      p_nprim = 0
      do ie = 1,n_ecp
        icent = i_cent_C(ie)
        if (o_ecp) then
          on_ecp = bas_isce(geom,ecpidin,icent)
        else
          on_ecp = .false.
        endif
        if (o_so) then
          on_so  = bas_isce(geom,soidin,icent)
        else
          on_so  = .false.
        endif
        if (on_ecp) then
          iucent = sf_ibs_ce2uce(icent,ecpid)
          F_cont = infbs_tags(TAG_FCONT,iucent,ecpid)
          L_cont = infbs_tags(TAG_LCONT,iucent,ecpid)
          do iucont = F_cont, L_cont
            type  = infbs_cont(CONT_TYPE,iucont,ecpid)
            nprim = infbs_cont(CONT_NPRIM,iucont,ecpid)
            ncoef = nprim*infbs_cont(CONT_NGEN,iucont,ecpid)
            if (nprim.ne.ncoef) then
              write(6,*)'general contraction ecp basis are invalid now'
              call errquit('int_ecp_build_ecp_ptrs: error',911,
     &       BASIS_ERR)
            endif
            iexp  = infbs_cont(CONT_IEXP,iucont,ecpid)
            icfp  = infbs_cont(CONT_ICFP,iucont,ecpid)
            irexp = infbs_cont(CONT_IREXP,iucont,ecpid)
cc AJL/Begin/SPIN ECPs
            mychan = infbs_cont(CONT_CHANNEL,iucont,ecpid)
            if (mychan.eq.channel.or.mychan.eq.0) then
cc AJL: If mychan = 0, applies the ECP to all channels
              if ((p_nprim+nprim).gt.nz_ecp) call errquit
     &            ('int_ecp_build_ecp_ptrs:too many coefficients',911,
     &             BASIS_ERR)
              call ecp_get_n3(
     &            c_exp((p_nprim+1)),
     &            dbl_mb(mb_exndcf(iexp,ecpid)),
     &            c_coef((p_nprim+1)),
     &            dbl_mb(mb_exndcf(icfp,ecpid)),
     &            dbl_mb(mb_exndcf(irexp,ecpid)),
     &            nprim,n0,n1,n2,n3,n4)
              n_prim_C(0,type,ie,1) = n0
              n_prim_C(1,type,ie,1) = n1
              n_prim_C(2,type,ie,1) = n2
              n_prim_C(3,type,ie,1) = n3
              n_prim_C(4,type,ie,1) = n4
              ind_C(type,ie,1)      = p_nprim+1
              ind_Z(type,ie,1)      = p_nprim+1
              n_coef_c(type,ie,1)   = nprim
              l_C(ie)               = max(type,l_C(ie))
              p_nprim = p_nprim+nprim
            endif
cc AJL/End
          enddo
        endif
        if (on_so) then
          iucent = sf_ibs_ce2uce(icent,soid)
          F_cont = infbs_tags(TAG_FCONT,iucent,soid)
          L_cont = infbs_tags(TAG_LCONT,iucent,soid)
          do iucont = F_cont, L_cont
            type  = infbs_cont(CONT_TYPE,iucont,soid)
            nprim = infbs_cont(CONT_NPRIM,iucont,soid)
            ncoef = nprim*infbs_cont(CONT_NGEN,iucont,soid)
            if (nprim.ne.ncoef) then
              write(6,*)'general contraction so basis are invalid now'
              call errquit('int_ecp_build_ecp_ptrs: error',911,
     &       BASIS_ERR)
            endif
            iexp  = infbs_cont(CONT_IEXP,iucont,soid)
            icfp  = infbs_cont(CONT_ICFP,iucont,soid)
            irexp = infbs_cont(CONT_IREXP,iucont,soid)
            if ((p_nprim+nprim).gt.nz_ecp) call errquit
     &          ('int_ecp_build_ecp_ptrs:too many coefficients',911,
     &          BASIS_ERR)
            call ecp_get_n3(
     &          c_exp((p_nprim+1)),
     &          dbl_mb(mb_exndcf(iexp,soid)),
     &          c_coef((p_nprim+1)),
     &          dbl_mb(mb_exndcf(icfp,soid)),
     &          dbl_mb(mb_exndcf(irexp,soid)),
     &          nprim,n0,n1,n2,n3,n4)
            n_prim_C(0,type,ie,2) = n0
            n_prim_C(1,type,ie,2) = n1
            n_prim_C(2,type,ie,2) = n2
            n_prim_C(3,type,ie,2) = n3
            n_prim_C(4,type,ie,2) = n4
            ind_C(type,ie,2)      = p_nprim+1
            ind_Z(type,ie,2)      = p_nprim+1
            n_coef_c(type,ie,2)   = nprim
            l_C(ie)               = max(type,l_C(ie))
            p_nprim = p_nprim+nprim
          enddo
        endif
      enddo
c
cc AJL/Begin/SPIN ECPS (Debug)
c      call print_ecp_ptrs(n_ecp,
c     &    l_ecp,nz_ecp,
c     &    channel, ! current channel [alpha/beta]
c     &    channels, ! total number of channels
c     &    n_prim_C,
c     &    n_coef_C,
c     &    ind_C,
c     &    ind_Z,
c     &    l_C,
c     &    c_exp,
c     &    c_coef)
cc AJL/End
c
      end
      subroutine print_ecp_ptrs(n_ecp,
     &    l_ecp,nz_ecp,
cc AJL/Begin/SPIN ECPs
     &    channel, ! current channel [alpha/beta]
     &    channels, ! total number of channels
cc AJL/End
     &    n_prim_C,
     &    n_coef_C,
     &    ind_C,
     &    ind_Z,
     &    l_C,
     &    c_exp,
     &    c_coef)
      implicit none
c     
      integer n_ecp    ! [input] number of ecp centers
      integer l_ecp    ! [input] maximum angular momentum in ecp basis
      integer nz_ecp   ! [input] number of prims/coeffs in stored data structure. 
cc AJL/Begin/SPIN ECPs
      integer channel  ! [input] Current spin channel used for ECPs 
      integer channels ! [input] Total number of spin channels used
cc AJL/End
      integer n_prim_C(0:4,-1:l_ecp,n_ecp,2) ! [output] 
      integer n_coef_C(-1:l_ecp,n_ecp,2)     ! [output] 
      integer ind_C(-1:l_ecp,n_ecp,2)        ! [output] 
      integer ind_Z(-1:l_ecp,n_ecp,2)        ! [output] 
      integer l_C(n_ecp)                     ! [output] 
      double precision c_exp(nz_ecp)         ! [output] 
      double precision c_coef(nz_ecp)        ! [output] 
*
      integer i, j, k, l
*      integer low0, high0, low1, high1, low2, high2
*      integer ir
      integer pe, pek
*
      write(6,*)' print_ecp_ptrs: start'
c
cc AJL/Begin/SPIN ECPs
      if (channels.eq.1) then ! ECP variable for number of channels
        write(6,*) 'Both Channels'
      else if (channel.eq.1) then ! Spin-polarised check?
        write(6,*) 'Alpha Channel'
      else if (channel.eq.2) then
        write(6,*) 'Beta Channel'
      endif
cc AJL/End
c
      write(6,*)' exponents and coefficients array'
      do i = 1,nz_ecp
        write(6,10000)i,c_exp(i),i,c_coef(i)
      enddo
10000 format(1x,'exp(',i5,') =',f12.6,2x,'coeff(',i5,') =',f12.6)
      do i = 1,n_ecp
        write(6,*)' l_c(',i,')',l_c(i)
      enddo
      do l = 1,2
        do i=1,n_ecp
          do j=-1,l_ecp
            write(6,*)' n_coef_C(',j,',',i,',',l,') = ',
     &          n_coef_C(j,i,l)
            write(6,*)' ind_C(',j,',',i,',',l,') = ',
     &          ind_C(j,i,l)
            write(6,*)' ind_z(',j,',',i,',',l,') = ',
     &          ind_z(j,i,l)
            do k = 0,4
              write(6,*)' n_prim_C(',
     &            k,',',j,',',i,',',l,') =',
     &            n_prim_C(k,j,i,l)
            enddo
          enddo
        enddo
      enddo
      write(6,*)' --- '
      write(6,*)' i = 1:',n_ecp
      do l = 1,2
        if (l.eq.1)
     &      write(6,*)'                                            ecp'
        if (l.eq.2)
     &      write(6,*)'                                            sop'
        do i=1,n_ecp
          write(6,*)' angular momentum max on center ',
     &        i,'is',l_c(i)
          write(6,*)' j = -1:',l_ecp
          do j=-1,l_ecp
            write(6,*)' lval = ',j
            pe = ind_c(j,i,l)-1
            do k=1,n_prim_C(0,j,i,l)
              pek = pe + k
              write(6,*)' 0 ',c_exp(pek),c_coef(pek)
            enddo
            pe = pe + n_prim_C(0,j,i,l)
            do k=1,n_prim_C(1,j,i,l)
              pek = pe + k
              write(6,*)' 1 ',c_exp(pek),c_coef(pek)
            enddo
            pe = pe + n_prim_C(1,j,i,l)
            do k=1,n_prim_C(2,j,i,l)
              pek = pe + k
              write(6,*)' 2 ',c_exp(pek),c_coef(pek)
            enddo
            pe = pe + n_prim_C(2,j,i,l)
            do k=1,n_prim_C(3,j,i,l)
              pek = pe + k
              write(6,*)' 3 ',c_exp(pek),c_coef(pek)
            enddo
            pe = pe + n_prim_C(3,j,i,l)
            do k=1,n_prim_C(4,j,i,l)
              pek = pe + k
              write(6,*)' 4 ',c_exp(pek),c_coef(pek)
            enddo
            write(6,*)' '
          enddo
          write(6,*)' '
        enddo
      enddo
      write(6,*)' print_ecp_ptrs: finish'
      end
      subroutine ecp_get_n3(
     &    c_exp,
     &    ecp_exp,
     &    c_coef,
     &    ecp_coef,
     &    grexp,
     &    nprim,new0,new1,new2,new3,new4)
      implicit none
#include "errquit.fh"
c
c data structures in the ecp code assume all 0 exponents, all 1 exponents
c and all 2 exponents are contiguous.  Not enforced in the input.
c 
c
      integer nprim, new0, new1, new2, new3, new4
      double precision c_exp(nprim)
      double precision ecp_exp(nprim)
      double precision c_coef(nprim)
      double precision ecp_coef(nprim)
      double precision grexp(nprim)
c
      integer i, j, ival, iptr
      integer new(0:4)
c
c assumes no general contraction in the ecp basis
c
*:debugn3:      write(6,*)' '
*:debugn3:      write(6,*)'-----------------------------------------------------'
*:debugn3:      write(6,*)' ecp exp '
*:debugn3:      call output (ecp_exp,1,nprim,1,1,nprim,1,1)
*:debugn3:      write(6,*)' ecp coef '
*:debugn3:      call output (ecp_coef,1,nprim,1,1,nprim,1,1)
*:debugn3:      write(6,*)' grexp '
*:debugn3:      call output (grexp,1,nprim,1,1,nprim,1,1)
      call ifill(5,0,new,1)
      iptr = 0
      do i = 0,4
        do j = 1,nprim
          ival = int(grexp(j) + 0.00002d00)
*:debugn3:          write(6,*)' i, j, ival ',i, j, ival
          if (ival.lt.0.or.ival.gt.4) then
            write(6,*)' ival    = ',ival
            write(6,*)' grexp(j) = ',grexp(j)
            call errquit
     &          ('ecp_get_n3: r-exponent not equal to 0,1,2,3,4',911,
     &       BASIS_ERR)
          endif
          if (ival.eq.i) then
            new(i) = new(i) + 1
            iptr = iptr + 1
*:debugn3:            write(6,*)' i, j, ival,iptr,new ',i, j, ival,iptr,new
            c_exp(iptr) = ecp_exp(j)
            c_coef(iptr) = ecp_coef(j)
          endif
        enddo
      enddo
*:debugn3:      write(6,*) ' iprt after loop,nprim ',iptr,nprim
      new0 = new(0)
      new1 = new(1)
      new2 = new(2)
      new3 = new(3)
      new4 = new(4)
*:debugn3:      write(6,*)' c_exp  '
*:debugn3:      call output (c_exp,1,nprim,1,1,nprim,1,1)
*:debugn3:      write(6,*)' c_coef  '
*:debugn3:      call output (c_coef,1,nprim,1,1,nprim,1,1)
*:debugn3:      write(6,*)'-----------------------------------------------------'
*:debugn3:      write(6,*)' '
c
      end
      block data ecp_init_bd
c
c Block data structure to initialize common block for ecp
c
#include "ecp_nwc.fh"
      data h_xyzecp      /0/       ! MA handle for ecp center coordinates
      data k_xyzecp      /0/       ! MA index  for ecp center coordinates
      data h_ecp_e       /0/       ! MA handle for ecp exponents
      data k_ecp_e       /0/       ! MA index  for ecp exponents
      data h_ecp_c       /0/       ! MA handle for ecp coefficients
      data k_ecp_c       /0/       ! MA index  for ecp coefficients
      data h_ecp_nprim_c /0/       ! MA handle for n_prim_C (see ecp_integral)
      data k_ecp_nprim_c /0/       ! MA index  for n_prim_C (see ecp_integral)
      data h_ecp_ncoef_c /0/       ! MA handle for n_coef_C (see ecp_integral)
      data k_ecp_ncoef_c /0/       ! MA index  for n_coef_C (see ecp_integral)
      data h_ecp_ind_c   /0/       ! MA handle for int_C (see ecp_integral)
      data k_ecp_ind_c   /0/       ! MA index  for int_C (see ecp_integral)
      data h_ecp_l_c     /0/       ! MA handle for l_C (see ecp_integral)
      data k_ecp_l_c     /0/       ! MA index  for l_C (see ecp_integral)
      data h_ecp_c2s     /0/       ! MA handle for c2s routines
      data k_ecp_c2s     /0/       ! MA index  for c2s routines
      data h_ecp_lip     /0/       ! MA handle for ecp center index
      data k_ecp_lip     /0/       ! MA index for ecp center index
cc AJL/Begin/SPIN ECPs
      data h_ecp_e_beta  /0/        ! MA handle for ecp exponents
      data k_ecp_e_beta  /0/        ! MA index  for ecp exponents
      data h_ecp_c_beta  /0/        ! MA handle for ecp coefficients
      data k_ecp_c_beta  /0/        ! MA index  for ecp coefficients
      data h_ecp_nprim_c_beta /0/   ! MA handle for n_prim_C (see ecp_integral)
      data k_ecp_nprim_c_beta /0/   ! MA index  for n_prim_C (see ecp_integral)
      data h_ecp_ncoef_c_beta /0/   ! MA handle for n_coef_C (see ecp_integral)
      data k_ecp_ncoef_c_beta /0/   ! MA index  for n_coef_C (see ecp_integral)
      data h_ecp_ind_c_beta   /0/   ! MA handle for ind_C (see ecp_integral)
      data k_ecp_ind_c_beta   /0/   ! MA index  for ind_C (see ecp_integral)
      data h_ecp_ind_z_beta   /0/   ! MA handle for ind_Z (see ecp_integral)
      data k_ecp_ind_z_beta   /0/   ! MA index  for ind_Z (see ecp_integral)
      data h_ecp_l_c_beta     /0/   ! MA handle for l_C (see ecp_integral)
      data k_ecp_l_c_beta     /0/   ! MA index  for l_C (see ecp_integral)
      data h_ecp_lip_beta     /0/   ! MA handle for ecp center index
      data k_ecp_lip_beta     /0/   ! MA index for ecp center index
cc AJL/End
      data n_zeta_c      /0/       ! length of ecp exp/coef array
      data l_ecp         /0/       ! high ang for ecp basis
      data n_ecp         /0/       ! number of ecp centers (from API)
cc AJL/Begin/SPIN ECPs
      data ecp_channels  /0/       ! number of ecp channels [Alpha/Beta]
cc AJL/End
      data init_ecp_init /.false./ ! logical saying if ecp is init-ed
      end
      subroutine int_ecp_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_ints,sz_ints,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_ints  ! [input] buffer size for ecp_ints
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_ints(sz_ints) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' int_ecp_hf1: coords a ',xyza
*        write(6,*)' int_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_integral(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),  ! new name is n_colc_C
     &      int_mb(k_ecp_ind_z),
     &      int_mb(k_ecp_ind_c),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_ints,sz_ints,1,   ! nblk 1 for ecp integrals only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intdd_ecp_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_grad,sz_grad,nat,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      integer nat
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_grad  ! [input] buffer size for ecp_grad
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_grad(sz_grad,3,3,(nat*(nat-1)/2+nat)) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN d_ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' intd_ecp_hf1: coords a ',xyza
*        write(6,*)' intd_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_hessian(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),
     &      int_mb(k_ecp_ind_z),
     &      int_mb(k_ecp_ind_c),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,1,nat,          ! nblock = 1 for ECP only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intd_ecp_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_grad,sz_grad,nat,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      integer nat
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_grad  ! [input] buffer size for ecp_grad
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_grad(sz_grad,3,nat) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN d_ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' intd_ecp_hf1: coords a ',xyza
*        write(6,*)' intd_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_gradient(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),
     &      int_mb(k_ecp_ind_z),
     &      int_mb(k_ecp_ind_c),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,1,nat,          ! nblock = 1 for ECP only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intso_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_ints,sz_ints,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_ints  ! [input] buffer size for ecp_ints
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_ints(sz_ints) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN so_hf1:',lscr
c
      call ecp_integral(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),  ! new name is n_colc_C
     &      int_mb(k_ecp_ind_z),
     &      int_mb(k_ecp_ind_c),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_ints,sz_ints,3,   ! nblk 1 for ecp integrals only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intd_so_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_grad,sz_grad,nat,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      integer nat
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_grad  ! [input] buffer size for ecp_grad
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_grad(sz_grad,3,nat) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN d_ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' intd_ecp_hf1: coords a ',xyza
*        write(6,*)' intd_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_gradient(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),
     &      int_mb(k_ecp_ind_z),
     &      int_mb(k_ecp_ind_c),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,3,nat,          ! nblock = 3 for SO only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intdd_so_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_grad,sz_grad,nat,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      integer nat
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_grad  ! [input] buffer size for ecp_grad
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_grad(sz_grad,3,3,(nat*(nat-1)/2+nat)) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN d_ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' intd_ecp_hf1: coords a ',xyza
*        write(6,*)' intd_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_hessian(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),
     &      int_mb(k_ecp_ind_z),
     &      int_mb(k_ecp_ind_c),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,3,nat,          ! nblock = 1 for ECP only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
c
cc AJL/Begin/SPIN ECPs
cc All the subroutines designed for use with the beta ECPs

      subroutine int_ecp_hf1_beta(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_ints,sz_ints,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_ints  ! [input] buffer size for ecp_ints
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_ints(sz_ints) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' int_ecp_hf1: coords a ',xyza
*        write(6,*)' int_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_integral(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e_beta),
     &      dbl_mb(k_ecp_c_beta),
     &      int_mb(k_ecp_nprim_c_beta),
     &      int_mb(k_ecp_ncoef_c_beta),  ! new name is n_colc_C
     &      int_mb(k_ecp_ind_z_beta),
     &      int_mb(k_ecp_ind_c_beta),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c_beta),
     &      int_mb(k_ecp_lip_beta), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_ints,sz_ints,1,   ! nblk 1 for ecp integrals only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intdd_ecp_hf1_beta(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_grad,sz_grad,nat,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      integer nat
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_grad  ! [input] buffer size for ecp_grad
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_grad(sz_grad,3,3,(nat*(nat-1)/2+nat)) 
                                         ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN d_ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' intd_ecp_hf1: coords a ',xyza
*        write(6,*)' intd_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_hessian(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e_beta),
     &      dbl_mb(k_ecp_c_beta),
     &      int_mb(k_ecp_nprim_c_beta),
     &      int_mb(k_ecp_ncoef_c_beta),
     &      int_mb(k_ecp_ind_z_beta),
     &      int_mb(k_ecp_ind_c_beta),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c_beta),
     &      int_mb(k_ecp_lip_beta), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,1,nat,          ! nblock = 1 for ECP only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
      subroutine intd_ecp_hf1_beta(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,ictra,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,ictrb,
     &    ecp_grad,sz_grad,nat,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La, ictra
      integer b_nprim, b_ngen, Lb, ictrb
      integer nat
      double precision expa(a_nprim), expb(b_nprim)
      double precision coefa(a_nprim,a_ngen), coefb(b_nprim,b_ngen)
      integer sz_grad  ! [input] buffer size for ecp_grad
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_grad(sz_grad,3,nat) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*
*      write(6,*)' lscr IN d_ecp_hf1:',lscr
*      if (.not.dryrun) then
*        write(6,*)' intd_ecp_hf1: coords a ',xyza
*        write(6,*)' intd_ecp_hf1: coords b ',xyzb
*      endif
c
      call ecp_gradient(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,ictra,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,ictrb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e_beta),
     &      dbl_mb(k_ecp_c_beta),
     &      int_mb(k_ecp_nprim_c_beta),
     &      int_mb(k_ecp_ncoef_c_beta),
     &      int_mb(k_ecp_ind_z_beta),
     &      int_mb(k_ecp_ind_c_beta),
     &      n_zeta_c,
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c_beta),
     &      int_mb(k_ecp_lip_beta), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,1,nat,          ! nblock = 1 for ECP only 
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
cc AJL/End
