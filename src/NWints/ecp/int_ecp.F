      logical function int_ecp_init(ecpidin,indx_grad)
* $Id: int_ecp.F,v 1.5 1996-10-20 06:10:07 d3e129 Exp $
      implicit none
*------------------------------------------------*
* COPYRIGHT (C) 1994, 1995, 1996,                *
*         Pacific Northwest National Laboratory, *
*         Battelle Memorial Institute.           *
*------------------------------------------------*
*------------> All Rights Reserved <-------------*
*------------------------------------------------*
#include "mafdecls.fh"
#include "bas.fh"
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
*::passed
      integer ecpidin   ! [input] ecp basis handle
      integer indx_grad ! [input] gradient level 0=energy 1=grad 2=hess
*::local
c
      integer ecpid                 ! lexical index for ecp basis handle
      integer geom                  ! geometry handle
      integer num_ecp               ! number of ecp centers (counted)
      integer ncenters              ! number of centers
      integer n_zeta_c_e            ! length of ecp exponent array
      integer n_zeta_c_c            ! length of ecp coefficient array
      integer nz_add                ! increment for exp/coef counters
      integer icent                 ! counter for center loop
      integer iucent                ! unique center of icent
      integer l_ecp_sz              ! size of ang info for ecp pointer arrays
      integer lmax_both
      integer idum
      integer basis_handle
      integer l_bas
      integer bas_id
      integer basis_ncont
      integer i_cont, ic, iuc, ie, icf, inp, ing, il, ig
      integer j_cont, jc, juc, je, jcf, jnp, jng, jl, jg
      integer lecp_mem, mem_ecp
      integer nV
      double precision dpdum(2)
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
        call errquit('int_ecp_init error',911)
      endif
      init_ecp_init = .true.
c.. lexical basis index
      ecpid = ecpidin + BASIS_HANDLE_OFFSET

c.. determine geom
      geom = ibs_geom(ecpid)
      
*      build up ecp information for heap
c....

c.. get number of ecp centers
      if (.not.geom_ncent_ecp(geom,n_ecp)) call errquit
     &    ('int_ecp: geom_ncent_ecp failed?',911)

c.. get number of centers
      if (.not.geom_ncent(geom,ncenters)) call errquit
     &    ('int_ecp: geom_ncent failed?',911)

c.. allocate space for coordintates of ecp centers
      if (.not.ma_alloc_get(mt_dbl,
     &    (3*n_ecp),
     &    'ecp center coords',
     &    h_xyzecp, k_xyzecp)) call errquit
     &    (' int_ecp: ecp center coords ma failed ',911)

c.. get coordinates for ecp centers.
      if (.not.geom_coords_ecp(geom,dbl_mb(k_xyzecp),n_ecp))
     &    call errquit('int_ecp: geom_coords_ecp failed',911)
*      write(6,*)' coordinates after read '
*      call output(dbl_mb(k_xyzecp),1,3,1,n_ecp,3,n_ecp,1)
c
c...       now comes the tricky part
c
c allocate and fill an exponent and coeff array for the ecp basis
c.. this MUST be changed when general contractions for ecp bases are used!!!!!
c
*      write(6,*)'inside ecp init'
      if (.not.ecp_print(ecpidin)) stop ' ecp_print error'
      num_ecp = 0
      n_zeta_c_e = 0
      n_zeta_c_c = 0
      do icent=1,ncenters
*        write(6,*)' icent = ',icent,oecpcent(icent,geom)
        if (oecpcent(icent,geom)) then
          num_ecp = num_ecp + 1
          iucent = sf_ibs_ce2uce(icent,ecpid)
*          write(6,*)' iucent',iucent
          nz_add = infbs_tags(TAG_NPRIM,iucent,ecpid)
*          write(6,*)' nz_add 1 exp ',nz_add
          n_zeta_c_e = n_zeta_c_e + nz_add
          nz_add = infbs_tags(TAG_NCOEF,iucent,ecpid)
*          write(6,*)' nz_add 2 coef ',nz_add
          n_zeta_c_c = n_zeta_c_c + nz_add
        endif
      enddo
*.. error check ecp centers
      if (num_ecp.ne.n_ecp) call errquit
     &    ('int_ecp_init: diff number of ecp centers',
     &    (num_ecp-n_ecp))

*.. error check for general contraction
      if (n_zeta_c_e.ne.n_zeta_c_c) then
        write(6,*)' possible general contraction on ecp basis'
        call errquit ('int_ecp_init: n_zeta_c_e .ne. n_zeta_c_c',
     &      (n_zeta_c_e- n_zeta_c_c))
      endif
      n_zeta_c = n_zeta_c_e ! use e one
c
c.. allocate space for ecp exponents and ecp coefficients
      if (.not.ma_alloc_get(mt_dbl,
     &    n_zeta_c,
     &    'ecp exponents',
     &    h_ecp_e, k_ecp_e)) call errquit
     &    (' int_ecp: ecp exponent ma failed ',911)
      if (.not.ma_alloc_get(mt_dbl,
     &    n_zeta_c,
     &    'ecp coefficients',
     &    h_ecp_c, k_ecp_c)) call errquit
     &    (' int_ecp: ecp coefficients ma failed ',911)
c
c... determine maximum angular momentum of ecp basis
      if (.not.bas_high_angular(ecpidin,l_ecp)) call errquit
     &    ('int_ecp: bas_high_angular failed',911)
c
      l_ecp_sz = l_ecp + 2       !    (-1->Lval)
c
c... allocate space for n_prim_C(0:2,-1:l_ecp_max,n_C),
      if (.not.ma_alloc_get(mt_int,
     &    (3*l_ecp_sz*n_ecp),
     &    'ecp n_prim_C',
     &    h_ecp_nprim_c, k_ecp_nprim_c)) call errquit
     &    (' int_ecp: ecp nprim_c ma failed ',911)
c... allocate space for n_coef_C(-1:l_ecp_max,n_C)
      if (.not.ma_alloc_get(mt_int,
     &    (l_ecp_sz*n_ecp),
     &    'ecp n_coef_C',
     &    h_ecp_ncoef_c, k_ecp_ncoef_c)) call errquit
     &    (' int_ecp: ecp ncoef_c ma failed ',911)
c... allocate space for n_tot_C      
      if (.not.ma_alloc_get(mt_int,
     &    n_ecp,
     &    'ecp n_tot_C',
     &    h_ecp_ntot_c, k_ecp_ntot_c)) call errquit
     &    (' int_ecp: ecp ncoef_c ma failed ',911)
c... allocate space for ind_C      
      if (.not.ma_alloc_get(mt_int,
     &    (l_ecp_sz*n_ecp),
     &    'ecp ind_C',
     &    h_ecp_ind_c, k_ecp_ind_c)) call errquit
     &      (' int_ecp: ecp ind_c ma failed ',911)
c... allocate space for l_C      
      if (.not.ma_alloc_get(mt_int,
     &    n_ecp,
     &    'ecp l_C',
     &    h_ecp_l_c, k_ecp_l_c)) call errquit
     &    (' int_ecp: ecp l_c ma failed ',911)
c... allocate space for ecp center pointer list
      if (.not.ma_alloc_get(mt_int,
     &    n_ecp,
     &    'ecp lexical indeces for ecp centers',
     &    h_ecp_lip, k_ecp_lip)) call errquit
     &    (' int_ecp: ecp_lip ma failed',911)
c
      call int_ecp_build_ecp_ptrs(ecpidin,
     &    ncenters,
     &    n_ecp,
     &    l_ecp,
     &    n_zeta_c,
     &    int_mb(k_ecp_nprim_c),
     &    int_mb(k_ecp_ncoef_c),
     &    int_mb(k_ecp_ind_c),
     &    int_mb(k_ecp_ntot_c),
     &    int_mb(k_ecp_l_c),
     &    int_mb(k_ecp_lip),
     &    dbl_mb(k_ecp_e),
     &    dbl_mb(k_ecp_c) )
c
*...   allocate space for c2s and s2c internal ecp transformation routines
c
c determine lmax among ao basis and ecp basis
c l_ecp currently has Lmax for ecp basis
      if (.not.ecp_get_parent_handle(ecpidin,basis_handle))
     &      call errquit
     &      ('int_ecp_init: ecp_get_parent_handle failed',911)
      if (.not.bas_high_angular(basis_handle,l_bas)) call errquit
     &      ('int_ecp_init: bas_high_angular failed for ao handle',
     &      911)
      lmax_both = l_ecp + l_bas + 2 
      call ecp_init_c2s(lmax_both,dpdum,dpdum,idum,1,1,.true.,mem_c2s)
      if (.not.ma_alloc_get(mt_dbl,
     &      mem_c2s,
     &      'ecp c2s routines',
     &      h_ecp_c2s, k_ecp_c2s)) call errquit
     &      ('int_ecp_init: ma failed for c2s',911)
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
      do i_cont = 1,basis_ncont
        ic  = sf_ibs_cn2ce(i_cont,bas_id)
        iuc = sf_ibs_cn2ucn(i_cont,bas_id)
        ie  = infbs_cont(cont_iexp,iuc,bas_id)
        icf = infbs_cont(cont_icfp,iuc,bas_id)
        inp = infbs_cont(cont_nprim,iuc,bas_id)
        ing = infbs_cont(cont_ngen,iuc,bas_id)
        il  = infbs_cont(cont_type,iuc,bas_id) + indx_grad ! for grads
        do j_cont = 1,basis_ncont
*          if (.not.bas_print(basis_handle)) stop ' error'
*          if (.not.gbs_map_print(basis_handle)) stop ' error '
          jc = sf_ibs_cn2ce(j_cont,bas_id)
          juc = sf_ibs_cn2ucn(j_cont,bas_id)
          je  = infbs_cont(cont_iexp,juc,bas_id)
          jcf = infbs_cont(cont_icfp,juc,bas_id)
          jnp = infbs_cont(cont_nprim,juc,bas_id)
          jng = infbs_cont(cont_ngen,juc,bas_id)
          jl  = infbs_cont(cont_type,juc,bas_id) + indx_grad  ! one for grads
          nV  = int_nbf_x(il)*int_nbf_x(jl)
          lecp_mem = 90000
          call int_ecp_hf1(
     &        coords(1,ic,ig),
     &        dbl_mb(mb_exndcf(ie,bas_id)),
     &        dbl_mb(mb_exndcf(icf,bas_id)),
     &        inp,ing,il,
     &        coords(1,jc,jg),
     &        dbl_mb(mb_exndcf(je,bas_id)),
     &        dbl_mb(mb_exndcf(jcf,bas_id)),
     &        jnp,jng,jl,
     &        dpdum,nV,dpdum,lecp_mem,.true.)
          lecp_mem = lecp_mem + nV  ! nint block used in api routine int_hf1sp
          mem_ecp = max(mem_ecp,lecp_mem)
*          write(6,*)' i_cont ',i_cont
*          write(6,*)' j_cont ',j_cont
*          write(6,*)' lecp_mem from int_ecp_init is ',
*     &        lecp_mem,mem_ecp
        enddo
      enddo
      write(6,*)' scr 1e memory without ecp',mem_1e
      mem_1e = max(mem_1e,mem_ecp)
      write(6,*)' scr 1e memory with    ecp',mem_1e
*
      int_ecp_init = .true.
*
      end
      subroutine int_ecp_terminate()
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"

      logical status

      status = .true.
      status = status .and. MA_free_heap(h_xyzecp)
      status = status .and. MA_free_heap(h_ecp_e)
      status = status .and. MA_free_heap(h_ecp_c)
      status = status .and. MA_free_heap(h_ecp_nprim_c)
      status = status .and. MA_free_heap(h_ecp_ncoef_c)
      status = status .and. MA_free_heap(h_ecp_ntot_c)
      status = status .and. MA_free_heap(h_ecp_ind_c)
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
      h_ecp_ntot_c  = 0
      k_ecp_ntot_c  = 0
      h_ecp_ind_c   = 0
      k_ecp_ind_c   = 0
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
      if (status) return
      call errquit
     &    ('int_ecp_terminate: error freeing heap',911)
      end
      subroutine int_ecp(
     &    basisina, basisinb, ecpidin,
     &    sha,xyza,shb,xyzb,
     &    ecp_ints,sz_ints,
     &    scr,lscr,dryrun)
      implicit none
*rak:#include "basP.fh"
*rak:#include "basdeclsP.fh"
*rak:#include "geobasmapP.fh"
*rak:#include "mafdecls.fh"
*rak:#include "bas_exndcf_dec.fh"
*rak:#include "bas_ibs_dec.fh"
*rak:#include "ecp_nwc.fh"
*rak:#include "geom.fh"
*rak:*
      integer basisina ! [input] basis set handle a center
      integer basisinb ! [input] basis set handle b center
      integer ecpidin  ! [input] ecp handle
      integer sha, shb ! [input] lexical index of shells a and b
      integer sz_ints  ! [input] buffer size for ecp_ints
      integer lscr     ! [input] length of scratch array
      double precision xyza(3), xyzb(3)  ! [input] a and b center coords.
      double precision ecp_ints(sz_ints) ! [output] ecp integrals
      double precision scr(lscr)         ! [scratch] array
      logical dryrun   ! [input] compute vs calculate memory requirements.
*rak:*
*rak:      integer basisa ! lexical basis index center a
*rak:      integer basisb ! lexical basis index center b
*rak:      integer ecpid  ! lexical ecp index
*rak:      integer geo_bas, geo_ecp, geoma, geomb, geom ! geometry handles 
*rak:      integer ucont
*rak:      integer la, npa, nga, iea, ica
*rak:      integer lb, npb, ngb, ieb, icb
*rak:*
*rak:#include "bas_exndcf_sfn.fh"
*rak:#include "bas_ibs_sfn.fh"
*rak:*
*rak:      basisa = basisina + BASIS_HANDLE_OFFSET
*rak:      basisb = basisinb + BASIS_HANDLE_OFFSET
*rak:      ecpid  = ecpidin  + BASIS_HANDLE_OFFSET
*rak:*
*rak:      geo_bas = ibs_geom(basisa)
*rak:      geo_ecp = ibs_geom(ecpid)
*rak:      if (geo_bas.eq.geo_ecp) then
*rak:        geoma = geo_bas
*rak:      else
*rak:        write(6,*)' basis and ecp loaded with different geometries'
*rak:        write(6,*)' geo_bas: ',geo_bas
*rak:        write(6,*)' geo_ecp: ',geo_ecp
*rak:        call errquit(' int_ecp: error ',911)
*rak:      endif
*rak:      geo_bas = ibs_geom(basisb)
*rak:      geo_ecp = ibs_geom(ecpid)
*rak:      if (geo_bas.eq.geo_ecp) then
*rak:        geomb = geo_bas
*rak:      else
*rak:        write(6,*)' basis and ecp loaded with different geometries'
*rak:        write(6,*)' geo_bas: ',geo_bas
*rak:        write(6,*)' geo_ecp: ',geo_ecp
*rak:        call errquit(' int_ecp: error ',911)
*rak:      endif
*rak:      if (geoma.ne.geomb) call errquit
*rak:     &    (' different geometries for a and b center??',911)
*rak:*
*rak:      geom = geoma  ! use a
*rak:*
*rak:      ucont = sf_ibs_cn2ucn(sha,basisa)
*rak:      La  = infbs_cont(CONT_TYPE,ucont,basisa)
*rak:      Npa = infbs_cont(CONT_NPRIM,ucont,basisa)
*rak:      Nga = infbs_cont(CONT_NGEN,ucont,basisa)
*rak:      iea = infbs_cont(CONT_IEXP,ucont,basisa)
*rak:      ica = infbs_cont(CONT_ICFP,ucont,basisa)
*rak:
*rak:      ucont = sf_ibs_cn2ucn(shb,basisb)
*rak:      Lb  = infbs_cont(CONT_TYPE,ucont,basisb)
*rak:      Npb = infbs_cont(CONT_NPRIM,ucont,basisb)
*rak:      Ngb = infbs_cont(CONT_NGEN,ucont,basisb)
*rak:      ieb = infbs_cont(CONT_IEXP,ucont,basisb)
*rak:      icb = infbs_cont(CONT_ICFP,ucont,basisb)
*rak:c
*rak:      call ecp_integral(
*rak:     &    xyza,
*rak:     &    dbl_mb(mb_exndcf(iea,basisa)),
*rak:     &    dbl_mb(mb_exndcf(ica,basisa)),
*rak:     &    Npa,Nga,La,
*rak:     &    xyzb,
*rak:     &    dbl_mb(mb_exndcf(ieb,basisb)),
*rak:     &    dbl_mb(mb_exndcf(icb,basisa)),
*rak:     &    Npb,Ngb,Lb,
*rak:     &    dbl_mb(k_xyzecp),
*rak:     &    dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
*rak:     &    int_mb(k_ecp_nprim_c),
*rak:     &    int_mb(k_ecp_ncoef_c),
*rak:     &    int_mb(k_ecp_ind_c),
*rak:     &    int_mb(k_ecp_ntot_c),
*rak:     &    n_zeta_c,
*rak:     &    int_mb(k_ecp_l_c),
*rak:     &    n_ecp,l_ecp,
*rak:     &    0,ecp_ints,sz_ints,
*rak:     &    dryrun,scr,lscr)
*rak:*
      call errquit('int_ecp called why?',911)
      end
      subroutine int_ecp_build_ecp_ptrs(ecpidin,
     &    ncenters,
     &    n_ecp,
     &    l_ecp,
     &    nz_ecp,
     &    n_prim_C,
     &    n_coef_C,
     &    ind_C,
     &    n_tot_C,
     &    l_C,
     &    i_cent_C,
     &    c_exp,
     &    c_coef)
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
      integer ecpidin  ! [input] ecp basis set handle
      integer ncenters ! [input] number of centers
      integer n_ecp    ! [input] number of ecp centers
      integer l_ecp    ! [input] maximum angular momentum in ecp basis
      integer nz_ecp   ! [input] number of prims/coeffs in stored data structure. 
      integer n_prim_C(0:2,-1:l_ecp,n_ecp) ! [output] 
      integer n_coef_C(-1:l_ecp,n_ecp)     ! [output] 
      integer ind_C(-1:l_ecp,n_ecp)        ! [output] 
      integer n_tot_C(n_ecp)                  ! [output] 
      integer l_C(n_ecp)                      ! [output] 
      integer i_cent_C(n_ecp)                 ! [output]
      double precision c_exp(nz_ecp)          ! [output] 
      double precision c_coef(nz_ecp)         ! [output] 
c
      integer geom        ! geometry handle
      integer ecpid       ! lexical basis set handle
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
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
*      write(6,*)' inside build'
      ecpid = ecpidin + BASIS_HANDLE_OFFSET
      geom = ibs_geom(ecpid)
c... zero arrays
      call ifill((3*(l_ecp+2)*n_ecp),0,n_prim_C,1)
      call ifill(((l_ecp+2)*n_ecp),0,n_coef_C,1)
      call ifill(((l_ecp+2)*n_ecp),0,ind_C,1)
      call ifill(n_ecp,0,n_tot_C,1)
      call ifill(n_ecp,0,l_C,1)
c
      if (.not.ma_verify_allocator_stuff()) stop ' error'
      num_ecp = 0
      p_nprim = 0
*      write(6,*)' ncenters',ncenters
      do icent=1,ncenters
        if (oecpcent(icent,geom)) then
          num_ecp = num_ecp + 1
          i_cent_C(num_ecp) = icent
          iucent = sf_ibs_ce2uce(icent,ecpid)
          F_cont = infbs_tags(TAG_FCONT,iucent,ecpid)
          L_cont = infbs_tags(TAG_LCONT,iucent,ecpid)
*          write(6,*)' num_ecp :',num_ecp
*          write(6,*)' iucent  :',iucent
*          write(6,*)' f_cont  :',F_cont
*          write(6,*)' l_cont  :',L_cont
          do iucont = F_cont, L_cont
*            write(6,*)' iucont :',iucont
            type  = infbs_cont(CONT_TYPE,iucont,ecpid)
            nprim = infbs_cont(CONT_NPRIM,iucont,ecpid)
            ncoef = nprim*infbs_cont(CONT_NGEN,iucont,ecpid)
*            write(6,*)' type  :',type
*            write(6,*)' nprim :',nprim
*            write(6,*)' ncoef :',ncoef
            if (nprim.ne.ncoef) then
              write(6,*)'general contraction ecp basis are invalid now'
              call errquit('int_ecp_build_ecp_ptrs: error',911)
            endif
            iexp  = infbs_cont(CONT_IEXP,iucont,ecpid)
            icfp  = infbs_cont(CONT_ICFP,iucont,ecpid)
            irexp = infbs_cont(CONT_IREXP,iucont,ecpid)
            if ((p_nprim+nprim).gt.nz_ecp) call errquit
     &          ('int_ecp_build_ecp_ptrs:too many coefficients',911)
            call dcopy(nprim,
     &          dbl_mb(mb_exndcf(iexp,ecpid)),1,
     &          c_exp(p_nprim+1),1)
            call dcopy(nprim,
     &          dbl_mb(mb_exndcf(icfp,ecpid)),1,
     &          c_coef(p_nprim+1),1)
            call ecp_get_n3(dbl_mb(mb_exndcf(irexp,ecpid)),
     &          nprim,n0,n1,n2)
            n_prim_C(0,type,num_ecp) = n0
            n_prim_C(1,type,num_ecp) = n1
            n_prim_C(2,type,num_ecp) = n2
            ind_C(type,num_ecp)      = p_nprim+1
            n_coef_C(type,num_ecp)   = nprim
            n_tot_C(num_ecp)          = nprim + n_tot_C(num_ecp)
            l_C(num_ecp)              = max(type, l_C(num_ecp))
c.. update pointer to c_exp and c_coef
            p_nprim = p_nprim + nprim
          enddo
        endif
      enddo
c
*debug:      call print_ecp_ptrs(n_ecp,
*debug:     &    l_ecp,nz_ecp,
*debug:     &    n_prim_C,
*debug:     &    n_coef_C,
*debug:     &    ind_C,
*debug:     &    n_tot_C,
*debug:     &    l_C,
*debug:     &    c_exp,
*debug:     &    c_coef)
c
      end
      subroutine print_ecp_ptrs(n_ecp,
     &    l_ecp,nz_ecp,
     &    n_prim_C,
     &    n_coef_C,
     &    ind_C,
     &    n_tot_C,
     &    l_C,
     &    c_exp,
     &    c_coef)
      implicit none
c     
      integer n_ecp    ! [input] number of ecp centers
      integer l_ecp    ! [input] maximum angular momentum in ecp basis
      integer nz_ecp   ! [input] number of prims/coeffs in stored data structure. 
      integer n_prim_C(0:2,-1:l_ecp,n_ecp) ! [output] 
      integer n_coef_C(-1:l_ecp,n_ecp)     ! [output] 
      integer ind_C(-1:l_ecp,n_ecp)        ! [output] 
      integer n_tot_C(n_ecp)                  ! [output] 
      integer l_C(n_ecp)                      ! [output] 
      double precision c_exp(nz_ecp)          ! [output] 
      double precision c_coef(nz_ecp)         ! [output] 
*
      integer i, j, k
*      integer low0, high0, low1, high1, low2, high2
*      integer ir
      integer pe, pek
*
      write(6,*)' print_ecp_ptrs: start'
      write(6,*)' exponents and coefficients array'
      do i = 1,nz_ecp
        write(6,10000)i,c_exp(i),i,c_coef(i)
      enddo
10000 format(1x,'exp(',i5,') =',f12.6,2x,'coeff(',i5,') =',f12.6)
      do i = 1,n_ecp
        write(6,*)' n_tot_c(',i,')',n_tot_c(i)
      enddo
      do i = 1,n_ecp
        write(6,*)' l_c(',i,')',l_c(i)
      enddo
      do i=1,n_ecp
        do j=-1,l_ecp
          write(6,*)' n_coef_C(',j,',',i,') = ',n_coef_C(j,i)
          write(6,*)' ind_C(',j,',',i,') = ',ind_C(j,i)
          do k = 0,2
            write(6,*)' n_prim_C(',k,',',j,',',i,') =',n_prim_C(k,j,i)
          enddo
        enddo
      enddo
      write(6,*)' --- '
      write(6,*)' i = 1:',n_ecp
      do i=1,n_ecp
        write(6,*)' total number of coefficients on center',
     &      i,'is',n_tot_c(i)
        write(6,*)' angular momentum max on center ',
     &      i,'is',l_c(i)
        write(6,*)' j = -1:',l_ecp
        do j=-1,l_ecp
          write(6,*)' lval = ',j
          pe = ind_c(j,i)-1
          do k=1,n_prim_C(0,j,i)
            pek = pe + k
            write(6,*)' 0 ',c_exp(pek),c_coef(pek)
          enddo
          pe = pe + n_prim_C(0,j,i)
          do k=1,n_prim_C(1,j,i)
            pek = pe + k
            write(6,*)' 1 ',c_exp(pek),c_coef(pek)
          enddo
          pe = pe + n_prim_C(1,j,i)
          do k=1,n_prim_C(2,j,i)
            pek = pe + k
            write(6,*)' 2 ',c_exp(pek),c_coef(pek)
          enddo
          write(6,*)' '
        enddo
        write(6,*)' '
      enddo
      write(6,*)' print_ecp_ptrs: finish'
      end
      subroutine ecp_get_n3(rexp,nprim,new0,new1,new2)
      implicit none
c
      integer nprim, new0, new1, new2
      double precision rexp(nprim)
c
      integer i, ival
c
      new0 = 0
      new1 = 0
      new2 = 0
      do i = 1,nprim
        ival = int(rexp(i) + 0.00002d00)
        if      (ival.eq.0) then
          new0 = new0 + 1
        else if (ival.eq.1) then
          new1 = new1 + 1
        else if (ival.eq.2) then
          new2 = new2 + 1
        else
          write(6,*)' ival    = ',ival
          write(6,*)' rexp(i) = ',rexp(i)
          call errquit('ecp_get_n3: r-exponent not equal to 0,1,2',911)
        endif
      enddo
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
      data h_ecp_ntot_c  /0/       ! MA handle for n_tot_C (see ecp_integral)
      data k_ecp_ntot_c  /0/       ! MA index  for n_tot_C (see ecp_integral)
      data h_ecp_ind_c   /0/       ! MA handle for int_C (see ecp_integral)
      data k_ecp_ind_c   /0/       ! MA index  for int_C (see ecp_integral)
      data h_ecp_l_c     /0/       ! MA handle for l_C (see ecp_integral)
      data k_ecp_l_c     /0/       ! MA index  for l_C (see ecp_integral)
      data h_ecp_c2s     /0/       ! MA handle for c2s routines
      data k_ecp_c2s     /0/       ! MA index  for c2s routines
      data h_ecp_lip     /0/       ! MA handle for ecp center index
      data k_ecp_lip     /0/       ! MA index for ecp center index
      data n_zeta_c      /0/       ! length of ecp exp/coef array
      data l_ecp         /0/       ! high ang for ecp basis
      data n_ecp         /0/       ! number of ecp centers (from API)
      data init_ecp_init /.false./ ! logical saying if ecp is init-ed
      end
      subroutine int_ecp_hf1(
     &    xyza,expa,coefa,a_nprim,a_ngen,La,
     &    xyzb,expb,coefb,b_nprim,b_ngen,Lb,
     &    ecp_ints,sz_ints,scr,lscr,dryrun)
      implicit none
#include "mafdecls.fh"
#include "ecp_nwc.fh"
*
      integer a_nprim, a_ngen, La
      integer b_nprim, b_ngen, Lb
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
c
      call ecp_integral(
     &      xyza,
     &      expa,
     &      coefa,
     &      a_nprim,a_ngen,La,
     &      xyzb,
     &      expb,
     &      coefb,
     &      b_nprim,b_ngen,Lb,
     &      dbl_mb(k_xyzecp),
     &      dbl_mb(k_ecp_e),dbl_mb(k_ecp_c),
     &      int_mb(k_ecp_nprim_c),
     &      int_mb(k_ecp_ncoef_c),
     &      int_mb(k_ecp_ind_c),
     &      int_mb(k_ecp_ntot_c),
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_ints,sz_ints,
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
*      write(6,*)' lscr IN ecp_hf1:',lscr
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
     &      int_mb(k_ecp_ind_c),
     &      int_mb(k_ecp_ntot_c),
     &      n_zeta_c,
     &      int_mb(k_ecp_l_c),
     &      int_mb(k_ecp_lip), 
     &      n_ecp,l_ecp,
     &      0,
     &      dbl_mb(k_ecp_c2s),mem_c2s,
     &      ecp_grad,sz_grad,nat,
     &      dryrun,scr,lscr,
     &      0)  ! ibug
*
      end
