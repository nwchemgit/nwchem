* $Id: atscf_ecp.F,v 1.1 1996-12-13 22:30:14 d3e129 Exp $
      subroutine atscf_ecp(geom,basis,tag,hatom,nbas,
     &    zeta,cont,ucont,nz)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
*****************************************#include "ecp_nwc.fh"
*
*.... computes the ecp contribution to the atomic scf
*
      integer geom  ! [input] geometry handle
      integer basis ! [input] basis set handle
      character*16  tag ! [input] tag for atom
      double precision hatom(*) ! [atom] ecp contribution to the potential
      integer nbas(4)
      integer nz
      double precision zeta(nz)
      double precision cont(nz)
      double precision ucont(nz)
c
      integer bgeom ! geometry handle that basis was loaded with
      integer nat   ! number of atoms
      integer icent ! center index
      character*16 gtag ! dummy geometry tag
      character*16 btag ! dummy basis tag
      integer clo   ! first non-uniuqe contraction on atom
      integer chi   ! last non-uniuqe contraction on atom
      integer i,j
      integer sz_hatom
      integer ecpis
c
      gtag = ' '
      btag = ' '
c
c compare geometry handles
c
      if (.not.bas_geom(basis,bgeom)) call errquit
     &    ('atscf_ecp: bas_geom failed ',911)
      if (geom.ne.bgeom) call errquit
     &    ('atscf_ecp: geom .ne. bgeom Fatal mismatch error',911)
c
c determine the center index
c
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('atscf_ecp: geom_ncent failed ',911)
c
      do icent = 1,nat
        if (.not.geom_cent_tag(geom, icent, gtag)) call errquit
     &      ('atscf_ecp: geom_cent_tag failed ',911)
        if (gtag.eq.tag) goto 00001
      enddo
00001 continue
*debug:      write(6,*)'tags 1', tag, gtag
c
c determine size of hatom
c
      sz_hatom = 0
      do i = 1,4
        j = nbas(i)
        sz_hatom = sz_hatom + j*(j+1)/2
      enddo
*debug:  write(6,*)' atscf_ecp: size of hatom is ',sz_hatom
      call dfill(sz_hatom,0.0d00,hatom,1)
c
c check if tag is an ecp center ?
c
      if (.not.geom_ecp_get(geom,icent)) return
c
c get ecpis .. ecp basis set handle
c
      if (.not.bas_get_ecp_handle(basis,ecpis)) call errquit
     &    ('atscf_ecp: bas_get_ecp_handle failed ',911)
c
c gather info for call to ecp_integral
c
      if (.not.bas_ce2cnr(ecpis,icent,clo,chi))
     &    call errquit('atscf_ecp:bas_ce2cnr failed',911)
c get tag from contraction      
      if (.not.bas_cont_tag(ecpis,clo,btag)) call errquit
     &    ('atscf_ecp:bas_cont_tag failed ',911)
      if (btag.ne.gtag) then
        write(luout,*)' atscf_ecp: tag mismatch '
        write(luout,*)' geometry tag: ',gtag
        write(luout,*)' basis    tag: ',btag
        call errquit('atscf_ecp: tag mismatch fatal error',911)
      endif
      call atscf_ecp_setup_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent)
      end
      subroutine atscf_ecp_setup_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
c::-passed
      integer geom
      integer basis
      integer ecpis
      integer nbas(4)
      integer nz
      double precision zeta(nz)
      double precision cont(nz)
      double precision ucont(nz)
      integer sz_hatom
      double precision hatom(sz_hatom)
      integer icent
c::-local
      integer sz_zeta
      integer lmaxbs
      integer lmax_all_ecp
      integer lscr, lbuf
      integer h_zetac,  k_zetac
      integer h_coefc,  k_coefc
      integer h_nprimc, k_nprimc
      integer h_ncoefc, k_ncoefc
      integer h_indc,   k_indc
      integer h_scr,    k_scr
      integer h_buf,    k_buf
c
c...  allocate space for info on ecp functions
c.. n_zeta_c number of exponents, need a max
c.1. zeta_c = exponents of ecp      (flat for all ecp info)
c.2. coef_c = contraction coeefs    (flat for all ecp info)
c.3. n_prim_c = info about primitives rval, type and center (0:2,-1:lmax,necp=1)
c.4. n_coef_c = info about prims summed over rval just type and center (-1:lmax,necp=1)
c.5. ind_c    = first coeff/exponent for each type (-1:lmax,necp=1)
c.. n_tot_C  = number of coeffs for each ecp center (necp = 1)
c.. n_zeta_C = total number ecp exponents summed over centers
c.. l_C      = max projector on center C (necp=1)
c.. n_c      = number of ecp centers necp=1=n_c
c.. l_ecp_max = max l val of projector on any center (lmax)
c
c.. only need 5 arrays use arbitrary length for zeta_c and coef_c with max for test
c.. get lmax in ecp basis and use this.  consistent array demension may hose it.
c.. 
c get ecp primitive exponents in ma
c
      sz_zeta = 500 ! should be okay for most atoms
      
c..   allocate memory for exponents
      if (.not.
     &    ma_alloc_get(mt_dbl,sz_zeta,
     &    'atscf ecp exponents',
     &    h_zetac,k_zetac)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
c..   allocate memory for coeffs
      if (.not.
     &    ma_alloc_get(mt_dbl,sz_zeta,
     &    'atscf ecp coeffs',
     &    h_coefc,k_coefc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
c.. get lmax for all ecps
      if (.not.bas_high_angular(ecpis,lmax_all_ecp)) call errquit
     &    ('atscf_ecp_setup_hatom: bas_high_angular failed',911)
c.. get lmax for all basis set
      if (.not.bas_high_angular(basis,lmaxbs)) call errquit
     &    ('atscf_ecp_setup_hatom: bas_high_angular failed',911)
      if (.not.
     &    ma_alloc_get(mt_int,(3*(lmax_all_ecp+2)),
     &    'atscf ecp nprimc',
     &    h_nprimc,k_nprimc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
      if (.not.
     &    ma_alloc_get(mt_int,(lmax_all_ecp+2),
     &    'atscf ecp nprimc',
     &    h_ncoefc,k_ncoefc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
      if (.not.
     &    ma_alloc_get(mt_int,(lmax_all_ecp+2),
     &    'atscf ecp nprimc',
     &    h_indc,k_indc)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
c..   allocate memory for scratch
      lbuf = (lmaxbs+1)*(lmaxbs+2)/2
      lbuf = lbuf*lbuf
      if (.not.
     &    ma_alloc_get(mt_dbl,lbuf,
     &    'atscf ecp integrals',
     &    h_buf,k_buf)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
c..   allocate memory for scratch
      lscr = max(lbuf*10,20000)
      if (.not.
     &    ma_alloc_get(mt_dbl,lscr,
     &    'atscf ecp scr',
     &    h_scr,k_scr)) call errquit
     &    ('atscf_ecp_setup_hatom: ma_alloc get failed',911)
c
      call atscf_ecp_build_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent,
     &    sz_zeta,dbl_mb(k_zetac),dbl_mb(k_coefc),
     &    lmax_all_ecp,
     &    int_mb(k_nprimc),int_mb(k_ncoefc),int_mb(k_indc),
     &    dbl_mb(k_scr),lscr,dbl_mb(k_buf),lbuf)

      if (.not.ma_free_heap(h_zetac)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_zetac',911)
      if (.not.ma_free_heap(h_coefc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_coefc',911)
      if (.not.ma_free_heap(h_nprimc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_nprimc',911)
      if (.not.ma_free_heap(h_ncoefc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_ncoefc',911)
      if (.not.ma_free_heap(h_indc)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_indc',911)
      if (.not.ma_free_heap(h_buf)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_buf',911)
      if (.not.ma_free_heap(h_scr)) call errquit
     &    ('atscf_ecp: ma_free_heap failed for h_scr',911)
      end
      subroutine atscf_ecp_build_hatom(geom,basis,ecpis,
     &    hatom,sz_hatom,nbas,zeta,cont,ucont,nz,icent,
     &    sz_zetac,zeta_c,coef_c,lmax_all,n_prim_c,
     &    n_coef_c, ind_c, scr,lscr, ecp_ints, l_ecp_ints)
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "basP.fh"
#include "bas_ibs_dec.fh"
#include "bas_exndcf_dec.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "ecp_nwc.fh"
c
c::-passed
      integer geom
      integer basis
      integer ecpis
      integer sz_hatom
      double precision hatom(sz_hatom)
      integer nbas(4)
      integer nz
      double precision zeta(nz)
      double precision cont(nz)
      double precision ucont(nz)
      integer icent
      integer sz_zetac
      double precision zeta_c(sz_zetac)
      double precision coef_c(sz_zetac)
      integer lmax_all
      integer n_prim_c(0:2,-1:lmax_all)
      integer n_coef_c(-1:lmax_all)
      integer ind_c(-1:lmax_all)
      integer lscr
      double precision scr(lscr)
      integer l_ecp_ints
      double precision ecp_ints(l_ecp_ints)
c::-local
      double precision C(3)
      integer type
      integer ncoef
      integer nprim
      integer ucent
      integer nn, f_cont, l_cont, nn_off
      integer ecp
      integer atscf_n_zetac
      integer n0, n1, n2
      integer iexp,icfp,irexp
      integer n_tot_c
      integer l_c
      integer i, j, k, lj, lk, ljkoff
      double precision expa, coefa, expb, coefb
      integer sz_ints, nbf_s, nbf_x
      integer lscr_guess
      integer cnt_hatom
c::-statement functions
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      ecp = ecpis + Basis_handle_offset
      ucent = sf_ibs_ce2uce(icent,ecp)
*debug:      write(6,*)' icent/ucent ',icent,'/',ucent
c      
c.. get number of contractions/primitives on tag
      atscf_n_zetac = infbs_tags(Tag_Nprim,ucent,ecp)
      if (atscf_n_zetac.gt.sz_zetac) then
        write(luout,*)' nprimc/ncoefc array too small '
        write(luout,*)'atscf_n_zetac  = ',atscf_n_zetac
        write(luout,*)'sz_zetac = ',sz_zetac
        write(luout,*)'contact nwchem-support@emsl.pnl.gov'
        call errquit
     &      ('atscf_ecp_build_hatom: fatal error: ',911)
      endif
      f_cont = infbs_tags(Tag_Fcont,ucent,ecp)
      l_cont = infbs_tags(Tag_Lcont,ucent,ecp)
*debug:      write(6,*)' f/l cont',f_cont,'/',l_cont
c
      l_c = -1000
      n_tot_c = 0
      nn_off = 1
      do nn = f_cont,l_cont
c
        type = infbs_cont(Cont_Type,nn,ecp)
        nprim = infbs_cont(Cont_Nprim,nn,ecp)
        ncoef = nprim*infbs_cont(Cont_Ngen,nn,ecp)
        if (nprim.ne.ncoef) then
          write(luout,*)'general contraction ecp basis are invalid now'
          call errquit('atscf_ecp_build_hatom: error',911)
        endif
        iexp  = infbs_cont(Cont_Iexp,nn,ecp)
        icfp  = infbs_cont(Cont_Icfp,nn,ecp)
        irexp = infbs_cont(Cont_Irexp,nn,ecp)
        if ((nn_off+nprim-1).gt.sz_zetac) call errquit
     &  ('atscf_ecp_build_hatom: too many exponents/coefficents',911)
        call dcopy(nprim,dbl_mb(mb_exndcf(iexp,ecp)),1,
     &      zeta_c(nn_off),1)
        call dcopy(nprim,dbl_mb(mb_exndcf(icfp,ecp)),1,
     &      coef_c(nn_off),1)
        call ecp_get_n3(dbl_mb(mb_exndcf(irexp,ecp)),nprim,n0,n1,n2)
        n_prim_c(0,type) = n0
        n_prim_c(1,type) = n1
        n_prim_c(2,type) = n2
        ind_c(type)      = nn_off
        n_coef_c(type)   = nprim
        n_tot_c          = n_tot_c + nprim
        l_c              = max(type,l_c)
        nn_off = nn_off + nprim
      enddo
c
*debug:      write(6,*)' exponent list for ecp '
*debug:      call output(zeta_c,1,atscf_n_zetac,1,1,atscf_n_zetac,1,1)
*debug:      write(6,*)' contraction list for ecp '
*debug:      call output(coef_c,1,atscf_n_zetac,1,1,atscf_n_zetac,1,1)
c      
      call dfill(3,0.0d00,C,1)
*debug:      write(6,*)'nbas',nbas
*debug:      write(6,*)'zeta',zeta
*debug:      write(6,*)'ucont',ucont
      ljkoff = 0
      cnt_hatom = 0
      do i = 1,4
        type = i-1
        sz_ints = (type+1)*(type+2)/2
        sz_ints = sz_ints*sz_ints
        lj = ljkoff
        do j = 1,nbas(i)
          lj = lj + 1
          expa = zeta(lj)
          coefa = ucont(lj)
          lk = ljkoff
          do k = 1,j
            lk=lk+1
            expb = zeta(lk)
            coefb = ucont(lk)
            lscr_guess = lscr
*debug:            write(6,*)' lscr ',lscr
            call ecp_integral(
     &          C,expa,coefa,1,1,type,
     &          C,expb,coefb,1,1,type,
     &          C,zeta_c,coef_c,n_prim_C,n_coef_c,
     &          ind_c, n_tot_c,
     &          atscf_n_zetac,
     &          l_c,1,lmax_all,0,
     &          dbl_mb(k_ecp_c2s),mem_c2s,
     &          ecp_ints,sz_ints,
     &          .true.,
     &          scr,lscr_guess,
     &          0)            
            if (lscr_guess.gt.lscr) then
              write(luout,*)' lscr_guess =',lscr_guess
              write(luout,*)' lscr       =',lscr
              write(luout,*)' contact nwchem-support@emsl.pnl.gov'
              call errquit('atscf_ecp_build_hatom: fatal error',911)
            else
              call ecp_integral(
     &            C,expa,coefa,1,1,type,
     &            C,expb,coefb,1,1,type,
     &            C,zeta_c,coef_c,n_prim_C,n_coef_c,
     &            ind_c, n_tot_c,
     &            atscf_n_zetac,
     &            l_c,1,lmax_all,0,
     &            dbl_mb(k_ecp_c2s),mem_c2s,
     &            ecp_ints,sz_ints,
     &            .false.,
     &            scr,lscr,
     &            0)            
*              write(6,*)' ecp integrals, cart'
*              call output(ecp_ints,1,sz_ints,1,1,sz_ints,1,1)
              nbf_x = (type+1)*(type+2)/2
              nbf_s = 2*type+1
              call spcart_tran1e(ecp_ints,scr,
     &            nbf_x,nbf_x,type,1,
     &            nbf_s,nbf_s,type,1,
     &            .false.)
*              write(6,*)' ecp integrals, spherical'
*              call output(ecp_ints,1,(nbf_s*nbf_s),
*     &            1,1,(nbf_s*nbf_s),1,1)
            endif
            cnt_hatom = cnt_hatom + 1
c...  use the zero component of the square spherical integral block
            nbf_x = (nbf_s*nbf_s - 1)/2 + 1 ! index into odd rank square matrix
            hatom(cnt_hatom) = ecp_ints(nbf_x)
          enddo
        enddo
        ljkoff = ljkoff + nbas(i)
      enddo
      end
