* $Id$
      subroutine atscf_ecp(geom,basis,tag,hatom,nbas,
     &    zeta,cont,ucont,nz)
*
*.... computes the ecp contribution to the atomic scf potential
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
c
c::-passed
      integer geom  ! [input] geometry handle
      integer basis ! [input] basis set handle
      character*16  tag ! [input] tag for atom
      double precision hatom(*) ! [output] ecp contribution to the V
      integer nbas(4) ! [input] atomic scf internal structure of 
*                               number of basis functions per sym
*                               or shell type s, p, d, f
      integer nz      ! [input] number of exponents and coefs
      double precision zeta(nz) ! [input] exponents
      double precision cont(nz) ! [input] scaled atomic scf coefs
      double precision ucont(nz)! [input] unscaled coefs (from bas)
c::-local
      integer bgeom ! geometry handle that basis was loaded with
      integer nat   ! number of atoms
      integer icent ! center index
      character*16 gtag ! dummy geometry tag
      character*16 btag ! dummy basis tag
      integer clo   ! first non-uniuqe contraction on atom
      integer chi   ! last non-uniuqe contraction on atom
      integer i,j   ! dummy index
      integer sz_hatom ! size of hatom array
      integer ecpis    ! ecp basis handle
c
cMG_start   pseudopotentials do not work for periodic systems yet
c
      integer isys
c
      if (.not. geom_systype_get(geom, isys)) call errquit
     $     ('atscf_ecp: systype?', 0, UNKNOWN_ERR)
      if (isys.gt.0) return
cMG_end  
c
      gtag = ' '
      btag = ' '
c
c compare geometry handles input vs. what bas used to read basis info
c
      if (.not.bas_geom(basis,bgeom)) call errquit
     &    ('atscf_ecp: bas_geom failed ',911, BASIS_ERR)
      if (geom.ne.bgeom) call errquit
     &    ('atscf_ecp: geom .ne. bgeom Fatal mismatch error',911,
     &       GEOM_ERR)
c
c determine the center index in the global picture
c
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('atscf_ecp: geom_ncent failed ',911, GEOM_ERR)
c
      do icent = 1,nat
        if (.not.geom_cent_tag(geom, icent, gtag)) call errquit
     &      ('atscf_ecp: geom_cent_tag failed ',911, GEOM_ERR)
        if (gtag.eq.tag) goto 00001
      enddo
00001 continue
*debug:      write(6,*)'tags 1', tag, gtag
c
c determine size of hatom (primitive triangles in s,p,d,f)
c
      sz_hatom = 0
      do i = 1,4
        j = nbas(i)
        sz_hatom = sz_hatom + j*(j+1)/2
      enddo
*debug:  write(6,*)' atscf_ecp: size of hatom is ',sz_hatom
      call dfill(sz_hatom,0.0d00,hatom,1)
c
c check if tag is an ecp center ?  If not do not compute anything
c
      if (.not.geom_ecp_get(geom,icent)) return
c
c get ecpis .. ecp basis set handle
c
      if (.not.bas_get_ecp_handle(basis,ecpis)) call errquit
     &    ('atscf_ecp: bas_get_ecp_handle failed ',911, BASIS_ERR)
c
c gather info for call to ecp_integral
c
      if (.not.bas_ce2cnr(ecpis,icent,clo,chi))
     &    call errquit('atscf_ecp:bas_ce2cnr failed',911, BASIS_ERR)
c Check if center is in the given ecp basis set
      if (.not.(bas_isce(geom,ecpis,icent))) then
        write(luout,*)' the lexical center ',icent, ' with tag ',
     &      gtag,' is not in the ecp basis '
        call errquit('atscf_ecp: tag mismatch fatal error',911,
     &       BASIS_ERR)
      endif
      call atscf_ecp_setup_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent)
      end
      subroutine atscf_ecp_setup_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent)
      implicit none
#include "errquit.fh"
*
* routine that sets up pointers for atomic ecp computation
*
#include "mafdecls.fh"
#include "bas.fh"
c::-passed
      integer geom  ! [input] geometry handle
      integer basis ! [input] basis set handle
      integer ecpis ! [input] ecp basis set handle
      integer nbas(4) ! [input] atomic scf internal structure of 
*                               number of basis functions per sym
*                               or shell type s, p, d, f
      integer nz      ! [input] number of exponents and coefs
      double precision zeta(nz) ! [input] exponents
      double precision cont(nz) ! [input] scaled atomic scf coefs
      double precision ucont(nz)! [input] unscaled coefs (from bas)
      integer sz_hatom ! [input] size of hatom array
      double precision hatom(sz_hatom) ! [output] ecp part of V
      integer icent ! [input] lexical index of center
c::-local
      integer sz_zeta ! size of exponent and coef arrays for ecp info
      integer lmaxbs  ! high angular momentum in basis set
      integer lmax_all_ecp ! high angular momentum in ecp basis
      integer lscr, lbuf   ! size of scratch and buffer arrays
      integer h_zetac,  k_zetac  ! MA handle and pointer ecp exponents
      integer h_coefc,  k_coefc  ! MA handle and pointer ecp coefs
      integer h_nprimc, k_nprimc ! MA handle and pointer ecp n_prim_c
      integer h_ncoefc, k_ncoefc ! MA handle and pointer ecp n_coef_c
      integer h_indc,   k_indc   ! MA handle and pointer ecp ind_c
      integer h_scr,    k_scr    ! MA handle and pointer scratch
      integer h_buf,    k_buf    ! MA handle and pointer buffer
c
c...  allocate space for info on ecp functions
c.. n_zeta_c number of exponents, need a max
c.1. zeta_c = exponents of ecp      (flat for all ecp info)
c.2. coef_c = contraction coeefs    (flat for all ecp info)
c.3. n_prim_c = info about primitives rval, type and center 
*               (0:4,-1:lmax,necp=1)
c.4. n_coef_c = info about prims summed over rval just type 
*               and center (-1:lmax,necp=1)
c.5. ind_c    = first coeff/exponent for each type (-1:lmax,necp=1)
c.. n_zeta_C = total number ecp exponents summed over centers
c.. l_C      = max projector on center C (necp=1)
c.. n_c      = number of ecp centers necp=1=n_c
c.. l_ecp_max = max l val of projector on any center (lmax)
c
c.. only need 5 arrays use arbitrary length for zeta_c and coef_c 
*   with max for test
c.. get lmax in ecp basis and use this.  
c.. 
c get ecp primitive exponents in ma
c
      sz_zeta = 500 ! should be okay for most atoms
      
c..   allocate memory for exponents
      if (.not.
     &    ma_alloc_get(mt_dbl,sz_zeta,
     &    'atscf ecp exponents',
     &    h_zetac,k_zetac)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911, MA_ERR)
      call dfill(sz_zeta,0.0d00,dbl_mb(k_zetac),1)

c..   allocate memory for coeffs
      if (.not.
     &    ma_alloc_get(mt_dbl,sz_zeta,
     &    'atscf ecp coeffs',
     &    h_coefc,k_coefc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911, MA_ERR)
      call dfill(sz_zeta,0.0d00,dbl_mb(k_coefc),1)

c.. get lmax for all ecps
      if (.not.bas_high_angular(ecpis,lmax_all_ecp)) call errquit
     &    ('atscf_ecp_setup_hatom: bas_high_angular failed',911,
     &       BASIS_ERR)

c.. get lmax for all basis set
      if (.not.bas_high_angular(basis,lmaxbs)) call errquit
     &    ('atscf_ecp_setup_hatom: bas_high_angular failed',911,
     &       BASIS_ERR)

c.. allocate memory for n_prim_c array (see ecp_integral for def)
      if (.not.
     &    ma_alloc_get(mt_int,(5*(lmax_all_ecp+2)),
     &    'atscf ecp nprimc',
     &    h_nprimc,k_nprimc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911,
     &       MA_ERR)
      call ifill((5*(lmax_all_ecp+2)),0,int_mb(k_nprimc),1)

c.. allocate memory for n_coef_c array (see ecp_integral for def)
      if (.not.
     &    ma_alloc_get(mt_int,(lmax_all_ecp+2),
     &    'atscf ecp ncoefc',
     &    h_ncoefc,k_ncoefc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911, MA_ERR)
      call ifill((lmax_all_ecp+2),0,int_mb(k_ncoefc),1)

c.. allocate memory for ind_c array (see ecp_integral for def)
      if (.not.
     &    ma_alloc_get(mt_int,(lmax_all_ecp+2),
     &    'atscf ecp nindc',
     &    h_indc,k_indc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911, MA_ERR)
      call ifill((lmax_all_ecp+2),0,int_mb(k_indc),1)

c..   allocate memory for ecp integral buffer array
      lbuf = (lmaxbs+1)*(lmaxbs+2)/2
      lbuf = lbuf*lbuf
      if (.not.
     &    ma_alloc_get(mt_dbl,lbuf,
     &    'atscf ecp integrals',
     &    h_buf,k_buf)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911, MA_ERR)
* does not need to be zeroed      

c..   allocate memory for scratch array used by ecp code
      lscr = max(lbuf*10,20000)
      if (.not.
     &    ma_alloc_get(mt_dbl,lscr,
     &    'atscf ecp scr',
     &    h_scr,k_scr)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911, MA_ERR)
* does not need to be zeroed      

c.. call routine that does computation
      call atscf_ecp_build_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent,
     &    sz_zeta,dbl_mb(k_zetac),dbl_mb(k_coefc),
     &    lmax_all_ecp,
     &    int_mb(k_nprimc),int_mb(k_ncoefc),int_mb(k_indc),
     &    dbl_mb(k_scr),lscr,dbl_mb(k_buf),lbuf)

c.. free allocated memory
      if (.not.ma_free_heap(h_zetac)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_zetac',911, MEM_ERR)
      if (.not.ma_free_heap(h_coefc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_coefc',911, MEM_ERR)
      if (.not.ma_free_heap(h_nprimc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_nprimc',911, MEM_ERR)
      if (.not.ma_free_heap(h_ncoefc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_ncoefc',911, MEM_ERR)
      if (.not.ma_free_heap(h_indc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_indc',911, MEM_ERR)
      if (.not.ma_free_heap(h_buf)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_buf',911, MEM_ERR)
      if (.not.ma_free_heap(h_scr)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_scr',911, MEM_ERR)
      end
      subroutine atscf_ecp_build_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent,
     &    sz_zetac,zeta_c,coef_c,lmax_all,n_prim_c,
     &    n_coef_c, ind_c, scr,lscr, ecp_ints, l_ecp_ints)
      implicit none
#include "errquit.fh"
*
* this routine fills pointer and coeff arrays with appropriate
* ecp information and then generates the ecp integrals
*
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "ecp_nwc.fh"
c
c::-passed
      integer geom  ! [input] geometry handle
      integer basis ! [input] basis set handle
      integer ecpis ! [input] ecp basis set handle
      integer nbas(4) ! [input] atomic scf internal structure of 
*                               number of basis functions per sym
*                               or shell type s, p, d, f
      integer nz      ! [input] number of exponents and coefs
      double precision zeta(nz) ! [input] exponents
      double precision cont(nz) ! [input] scaled atomic scf coefs
      double precision ucont(nz)! [input] unscaled coefs (from bas)
      integer sz_hatom ! [input] size of hatom array
      double precision hatom(sz_hatom) ! [output] ecp part of V
      integer icent ! [input] lexical index of center
      integer sz_zetac ! [input] size of ecp exp/coef arrays
* [s] == [scratch] used/filled locally and down tree but not 
*                  important above
* [i] == [input]
      double precision zeta_c(sz_zetac) ! [s] ecp exponent array
      double precision coef_c(sz_zetac) ! [s] ecp coefs array
      integer lmax_all                  ! [i] ecp max angular momentum
      integer n_prim_c(0:4,-1:lmax_all) ! [s] ecp code pointer array
      integer n_coef_c(-1:lmax_all)     ! [s] ecp code pointer array
      integer ind_c(-1:lmax_all)        ! [s] ecp code pointer array
      integer lscr     ! [input] length of scratch buffer
      double precision scr(lscr) ! [s] scratch buffer for ecp code
      integer l_ecp_ints         ! [input] size of ecp int buffer
      double precision ecp_ints(l_ecp_ints) ! [s] ecp int buffer
c::-local
      double precision C(3)  ! coordinates
      integer type           ! function type
      integer ncoef          ! number of coefs in contraction
      integer nprim          ! number of prims in contraction
      integer ucent          ! unique basis center
      integer nn             ! dummy contraction index
      integer nn_off         ! pointer offset for building ecp arrays
      integer f_cont         ! first contraction on center
      integer l_cont         ! last contraction on center 
      integer ecp            ! lexical basis set index for ecpis
      integer atscf_n_zetac  ! size of n_zeta_c in reality
      integer n0, n1, n2     ! dummy r exponent indecies
      integer n3, n4         ! dummy r exponent indecies
      integer iexp,icfp,irexp ! private basis set pointer info
      integer l_c            ! max ang on ecp center
* for atomic code the next three are the same = 1
      integer i_cent_C       ! lexical atom index of ecp center
      integer i_c_A          ! lexical atom index of basis function center A
      integer i_c_B          ! lexical atom index of basis function center B
      integer i, j, k,       ! loop indices 
     &    lj, lk, ljkoff     ! pointer offsets
      double precision
     &    expa, coefa, expb, coefb  ! prim exp and coefs
      integer sz_ints        ! size of ints block
      integer nbf_s          ! nbf of shell in spherical d=5
      integer nbf_x          ! nbf of shell in cart.     d=6
      integer lscr_guess     ! dummy arg to check memory in ecp code
      integer cnt_hatom      ! counter in hatom array
c::-statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c..  get lexical basis set index for ecp basis set handle
      ecp = ecpis + Basis_handle_offset

c..  get unique center index.
      ucent = sf_ibs_ce2uce(icent,ecp)
*debug:      write(6,*)' icent/ucent ',icent,'/',ucent
c      
c.. get number of contractions/primitives on ecp tag
      atscf_n_zetac = infbs_tags(Tag_Nprim,ucent,ecp)

c.. check to make sure ecp arrays are big enough
      if (atscf_n_zetac.gt.sz_zetac) then
        write(luout,*)' nprimc/ncoefc array too small '
        write(luout,*)'atscf_n_zetac  = ',atscf_n_zetac
        write(luout,*)'sz_zetac = ',sz_zetac
        write(luout,*)'contact nwchem-support@emsl.pnl.gov'
        call errquit
     &      ('atscf_ecp_build_hatom: fatal error: ',911, BASIS_ERR)
      endif

c.. get first/last contraction on ecp tag
      f_cont = infbs_tags(Tag_Fcont,ucent,ecp)
      l_cont = infbs_tags(Tag_Lcont,ucent,ecp)
*debug:      write(6,*)' f/l cont',f_cont,'/',l_cont
c
c.. intialize some variables
      l_c = -1000
      nn_off = 1

c.. loop over contractions
      do nn = f_cont,l_cont
c
        type = infbs_cont(Cont_Type,nn,ecp)
        nprim = infbs_cont(Cont_Nprim,nn,ecp)
        ncoef = nprim*infbs_cont(Cont_Ngen,nn,ecp)
        if (nprim.ne.ncoef) then
          write(luout,*)
     &        'general contraction ecp basis are invalid now'
          call errquit('atscf_ecp_build_hatom: error',911, BASIS_ERR)
        endif
        iexp  = infbs_cont(Cont_Iexp,nn,ecp)
        icfp  = infbs_cont(Cont_Icfp,nn,ecp)
        irexp = infbs_cont(Cont_Irexp,nn,ecp)
        if ((nn_off+nprim-1).gt.sz_zetac) call errquit
     &  ('atscf_ecp_build_hatom: too many exponents/coefficents',911,
     &       BASIS_ERR)
*rak:        call dcopy(nprim,dbl_mb(mb_exndcf(iexp,ecp)),1,
*rak:     &      zeta_c(nn_off),1)
*rak:        call dcopy(nprim,dbl_mb(mb_exndcf(icfp,ecp)),1,
*rak:     &      coef_c(nn_off),1)
        call ecp_get_n3(
     &      zeta_c(nn_off),
     &      dbl_mb(mb_exndcf(iexp,ecp)),
     &      coef_c(nn_off),
     &      dbl_mb(mb_exndcf(icfp,ecp)),
     &      dbl_mb(mb_exndcf(irexp,ecp)),nprim,n0,n1,n2,n3,n4)
        n_prim_c(0,type) = n0
        n_prim_c(1,type) = n1
        n_prim_c(2,type) = n2
        n_prim_c(3,type) = n3
        n_prim_c(4,type) = n4
        ind_c(type)      = nn_off
        n_coef_c(type)   = nprim
        l_c              = max(type,l_c)
        nn_off = nn_off + nprim
      enddo
      i_cent_c = 1
      i_c_A    = 1
      i_c_B    = 1
c
*debug:      write(6,*)' exponent list for ecp '
*debug:      call output(zeta_c,1,atscf_n_zetac,1,1,atscf_n_zetac,1,1)
*debug:      write(6,*)' contraction list for ecp '
*debug:      call output(coef_c,1,atscf_n_zetac,1,1,atscf_n_zetac,1,1)
c      
c atom centered at origin
      call dfill(3,0.0d00,C,1)
*debug:      write(6,*)'nbas',nbas
*debug:      write(6,*)'zeta',zeta
*debug:      write(6,*)'ucont',ucont
      ljkoff = 0
      cnt_hatom = 0
c.. loop over shell types or atomic scf symmetries
      do i = 1,4
        type = i-1   ! type in NWChem terms
        sz_ints = (type+1)*(type+2)/2
        sz_ints = sz_ints*sz_ints  ! size of integral block
        lj = ljkoff

c.. loop over primitives for <bra|
        do j = 1,nbas(i)
          lj = lj + 1
          expa = zeta(lj)
          coefa = ucont(lj)
          lk = ljkoff
c.. loop over primitives for |ket>
          do k = 1,j
            lk=lk+1
            expb = zeta(lk)
            coefb = ucont(lk)
            lscr_guess = lscr
*debug:            write(6,*)' lscr ',lscr

              call ecp_integral(
     &            C,expa,coefa,1,1,type,i_c_A,
     &            C,expb,coefb,1,1,type,i_c_B,
     &            C,zeta_c,coef_c,n_prim_C,n_coef_c,
     &            ind_c, ind_c, atscf_n_zetac, atscf_n_zetac,
     &            l_c,i_cent_C,1,lmax_all,0,
     &            dbl_mb(k_ecp_c2s),mem_c2s,
     &            ecp_ints,sz_ints,1,  ! nblk = 1 for ecp integrals
     &            .true.,
     &            scr,lscr_guess,
     &            0)            
c.. make sure scratch array is big enough
            if (lscr_guess.gt.lscr) then
              write(luout,*)' lscr_guess =',lscr_guess
              write(luout,*)' lscr       =',lscr
              write(luout,*)' contact nwchem-support@emsl.pnl.gov'
              call errquit('atscf_ecp_build_hatom: fatal error',911,
     &       BASIS_ERR)
            else
c.. compute ecp integrals <bra|ecp|ket>
              call ecp_integral(
     &            C,expa,coefa,1,1,type,i_c_A,
     &            C,expb,coefb,1,1,type,i_c_B,
     &            C,zeta_c,coef_c,n_prim_C,n_coef_c,
     &            ind_c, ind_c, atscf_n_zetac, atscf_n_zetac,
     &            l_c,i_cent_C,1,lmax_all,0,
     &            dbl_mb(k_ecp_c2s),mem_c2s,
     &            ecp_ints,sz_ints,1,      ! nblk = 1 for ecp integrals
     &            .false.,
     &            scr,lscr,
     &            0)            
*              write(6,*)' ecp integrals, cart'
*              call output(ecp_ints,1,sz_ints,1,1,sz_ints,1,1)
              nbf_x = (type+1)*(type+2)/2
              nbf_s = 2*type+1
              
c.. transform cartesian block to spherical since the atomic scf 
*             is a spherical code.  This makes the integral block
*             diagonal with the same values for all 2*l+1 components.

              call spcart_tran1e(ecp_ints,scr,
     &            nbf_x,nbf_x,type,1,
     &            nbf_s,nbf_s,type,1,
     &            .false.)
*              write(6,*)' ecp integrals, spherical'
*              call output(ecp_ints,1,(nbf_s*nbf_s),
*     &            1,1,(nbf_s*nbf_s),1,1)
            endif

c.. compute index into hatom array
            cnt_hatom = cnt_hatom + 1
c...  use the zero component of the (-l:l) square spherical integral 
*     block.  nbf_x is now the index into odd rank square matrix 
*             (e.g., find the center element)
            nbf_x = (nbf_s*nbf_s - 1)/2 + 1 
            hatom(cnt_hatom) = ecp_ints(nbf_x)
          enddo
        enddo

c.. compute offset into scalar atomic scf zeta and ucont arrays
        ljkoff = ljkoff + nbas(i)
      enddo
      end
