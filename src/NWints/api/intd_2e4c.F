C> \ingroup nwint
C> @{
C>
C> \brief Compute the 4-center 2-electron integral derivatives
C>
C> Compute the 4-center 2-electron integral derivatives as given by
C> \f{eqnarray*}{
C> \frac{\partial({\mu}{\rho}|{\nu}{\lambda})}{\partial X_x} = \int_{-\infty}^{\infty} \frac{\partial g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})}{\partial X_x}dr_{1}dr_{2}
C> \f}
C> The integral derivatives are stored in an order that is consistent with
C> the declaration `ERI(nint,ncoord,natom)`, where `nint` is the number of
C> integrals in the shell pair, `ncoord` is the number of Cartesian coordinates
C> and `natom` equals 4 for the 4 basis set expansion centers. The actual
C> lexical indeces of the atoms on which the shells `ish`, `jsh`, `ksh` and
C> `lsh` are centered are returned in `idatom`.
C>
      subroutine intd_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri, idatom)
c $Id$
      implicit none
c
c basic api routine to generate 4 center two electron integral derivatives 
c 
#include "stdio.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
c
c::external subroutines used
c errquit
c::functions
      integer int_nint_cart
      external int_nint_cart
ckw
      integer int_nint
      external int_nint
c
      logical cando_nw
      logical cando_txs
      logical cando_sim
      external cando_nw
      external cando_txs
      external cando_sim
ckw
c::passed
      integer brain !< [Input] basis set handle for bra basis
      integer ish   !< [Input] lexical contraction index 
      integer jsh   !< [Input] lexical contraction index 
      integer ketin !< [Input] basis set handle for ket basis
      integer ksh   !< [Input] lexical contraction index 
      integer lsh   !< [Input] lexical contraction index 
      integer lscr  !< [Input] length of scratch array
      integer leri  !< [Input] length of eri array
      double precision scr(lscr) !< [Scratch] scratch array for integral code.
      double precision eri(leri) !< [Output]  array for two electron integral derivatives.
c NOTE: length of idatom is always 4 because there can be at most 4 centers involved
      integer idatom(4)          !< [Output]  array identifying centers for derivatives
c
c                       ! e.g., the first  nint*3 derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
c                       !       the third  nint*3 derivatives go to center idatom(3)
c                       !       the fourth nint*3 derivatives go to center idatom(4)
c
c Order is...   nint*3*4 (3=> xyz, 4=atoms)
c
c  /                   |   
c | nint,  d <ij|kl>   |
c |      --------------|
c  \     d[idatom(1),x]|
c                          |
c       nint,  d <ij|kl>   |
c            --------------|
c            d[idatom(1),y]|
c                              |
c           nint,  d <ij|kl>   |
c                --------------|
c                d[idatom(1),z]|
c                                  |
c               nint,  d <ij|kl>   |
c                    --------------|
c                    d[idatom(2),x]|
c                                      |
c                   nint,  d <ij|kl>   |
c                        --------------|
c                        d[idatom(2),y]|
c                                           |
c                       nint,  d <ij|kl>    |
c                            -------------- |
c                            d[idatom(2),z] |
c                                              |   
c                           nint,  d <ij|kl>   |
c                                --------------|
c                                d[idatom(3),x]|
c                                                  |
c                               nint,  d <ij|kl>   |
c                                    --------------|
c                                    d[idatom(3),y]|
c                                                      |
c                                   nint,  d <ij|kl>   |
c                                        --------------|
c                                        d[idatom(3),z]|
c                                                          |
c                                       nint,  d <ij|kl>   |
c                                            --------------|
c                                            d[idatom(4),x]|
c                                                              |
c                                           nint,  d <ij|kl>   |
c                                                --------------|
c                                                d[idatom(4),y]|
c                                                                   \
c                                               nint,  d <ij|kl>     |
c                                                    --------------  |
c                                                    d[idatom(4),z] / 
c                     
c::local
      integer nint_ck, ucont
      integer bra, ket, ab_geom, cd_geom
      integer inp, igen, iexp, icf, icfs, itype, iatom
      integer jnp, jgen, jexp, jcf, jcfs, jtype, jatom
      integer knp, kgen, kexp, kcf, kcfs, ktype, katom
      integer lnp, lgen, lexp, lcf, lcfs, ltype, latom
c.rel-dmd
      logical status_rel, bra_rel, ket_rel
      logical i_rel, j_rel, k_rel, l_rel
      integer sbas, abas, bras, kets
ckw
      double precision roff(3)
      integer txs_i, txs_j, txs_k, txs_l
      logical status_nw, status_txs, status_sim
      logical dum_log
      integer nintzero, num_quart, dummy_lab
      double precision q4
      integer momentum
      integer scr_ptr,ii
ckw
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
#include "int_nbf.fh"
c
      i_rel = .false.
      j_rel = .false.
      k_rel = .false.
      l_rel = .false.
      nint_ck = int_nint_cart(brain,ish,brain,jsh,ketin,ksh,ketin,lsh)
      if (nint_ck*3*4.gt.leri) then
        write(luout,*) 'nint*3*4 = ',nint_ck*3*4
        write(luout,*) 'leri     = ',leri
        call errquit('intd_2e4c: nint>leri error',911, INT_ERR)
      endif
c
c  check if spherical/gencon/sp shell
c
      call int_nogencont_check(brain,'intd_2e4c:bra')
      call int_nogencont_check(ketin,'intd_2e4c:ket')
      call int_nospshell_check(brain,'intd_2e4c:bra')
      call int_nospshell_check(ketin,'intd_2e4c:ket')
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      bras = bra
      kets = ket
c
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      if (ab_geom.ne.cd_geom) then
        write(luout,*)
     &      'intd_2e4c.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_2e4c: geom error ',911, INT_ERR)
      endif

c... stat rel
      status_rel = dyall_mod_dir .and. .not.nesc_1e_approx
     &    .and. (brain .eq. ketin) .and. (brain .eq. ao_bsh)
      if (status_rel) then
c
c     get basis set handles; relativistic integral option valid
c     if bra or ket are the ao basis and bra and ket have both
c     functions relativistic
c
        bra_rel = .false.
        ket_rel = .false.
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        bras = sbas
        kets = sbas
        bra_rel = bra .eq. abas
        if (bra_rel) then
          ucont = sf_ibs_cn2ucn(ish,bra)
          i_rel = infbs_cont(CONT_RELLS ,ucont,bra) .ne. 0
          ucont = sf_ibs_cn2ucn(jsh,bra)
          j_rel = infbs_cont(CONT_RELLS ,ucont,bra) .ne. 0
          bra_rel = bra_rel .and. i_rel .and. j_rel
        end if
        ket_rel = ket .eq. abas
        if (ket_rel) then
          ucont = sf_ibs_cn2ucn(ksh,ket)
          k_rel = infbs_cont(CONT_RELLS ,ucont,ket) .ne. 0
          ucont = sf_ibs_cn2ucn(lsh,ket)
          l_rel = infbs_cont(CONT_RELLS ,ucont,ket) .ne. 0
          ket_rel = ket_rel .and. k_rel .and. l_rel
        end if
        status_rel = status_rel .and. (bra_rel .or. ket_rel)
      end if
c
      ucont = (sf_ibs_cn2ucn(ish,bra))
      inp   = infbs_cont(CONT_NPRIM,ucont,bra)
      igen  = infbs_cont(CONT_NGEN,ucont,bra)
      iexp  = infbs_cont(CONT_IEXP,ucont,bra)
      icf   = infbs_cont(CONT_ICFP,ucont,bra)
      itype = infbs_cont(CONT_TYPE,ucont,bra) 
      iatom = (sf_ibs_cn2ce(ish,bra))
      if (i_rel) ucont = ao_to_ls(ucont)
      icfs = infbs_cont(CONT_ICFP,ucont,bras)
c
      ucont = (sf_ibs_cn2ucn(jsh,bra))
      jnp   = infbs_cont(CONT_NPRIM,ucont,bra)
      jgen  = infbs_cont(CONT_NGEN,ucont,bra)
      jexp  = infbs_cont(CONT_IEXP,ucont,bra)
      jcf   = infbs_cont(CONT_ICFP,ucont,bra)
      jtype = infbs_cont(CONT_TYPE,ucont,bra) 
      jatom = (sf_ibs_cn2ce(jsh,bra))
      if (j_rel) ucont = ao_to_ls(ucont)
      jcfs = infbs_cont(CONT_ICFP,ucont,bras)
c
      ucont = (sf_ibs_cn2ucn(ksh,ket))
      knp   = infbs_cont(CONT_NPRIM,ucont,ket)
      kgen  = infbs_cont(CONT_NGEN,ucont,ket)
      kexp  = infbs_cont(CONT_IEXP,ucont,ket)
      kcf   = infbs_cont(CONT_ICFP,ucont,ket)
      ktype = infbs_cont(CONT_TYPE,ucont,ket) 
      katom = (sf_ibs_cn2ce(ksh,ket))
      if (k_rel) ucont = ao_to_ls(ucont)
      kcfs = infbs_cont(CONT_ICFP,ucont,kets)
c
      ucont = (sf_ibs_cn2ucn(lsh,ket))
      lnp   = infbs_cont(CONT_NPRIM,ucont,ket)
      lgen  = infbs_cont(CONT_NGEN,ucont,ket)
      lexp  = infbs_cont(CONT_IEXP,ucont,ket)
      lcf   = infbs_cont(CONT_ICFP,ucont,ket)
      ltype = infbs_cont(CONT_TYPE,ucont,ket) 
      latom = (sf_ibs_cn2ce(lsh,ket))
      if (l_rel) ucont = ao_to_ls(ucont)
      lcfs = infbs_cont(CONT_ICFP,ucont,kets)
c
c... new logic
      if (iatom.eq.jatom.and.jatom.eq.katom.and.katom.eq.latom) then
        call dcopy((nint_ck*3*4),0.0d00,0,eri,1)
        call ifill(4,-1,idatom,1)
        return
      endif
c
ckw--------
ckw
      momentum=abs(itype)+abs(jtype)+abs(ktype)+abs(ltype)
      status_txs =.true.
      if( momentum.lt.4) then 
        status_txs =.false.
      endif
      status_txs = status_txs .and.
     &    cando_txs(brain,ish,jsh).and.cando_txs(ketin,ksh,lsh)
      status_nw  =
     &    cando_nw(brain,ish,jsh).and.cando_nw(ketin,ksh,lsh) 
      if (lgen.gt.1.or.kgen.gt.1.or.jgen.gt.1.or.igen.gt.1) then
        status_txs = status_txs .and. .true.
        status_nw = .false.
      endif
      status_sim = cando_sim(brain,ish,jsh).and.cando_sim(ketin,ksh,lsh)
*************************************************************************
* texas one at a time is broke for now:  temporary fix
*************************************************************************
      status_txs = .false.
*************************************************************************
c
      if (status_txs .and. .not.status_rel) then
        call dcopy(3,0.0d00,0,roff,1)
        q4 = 1.0d00
        txs_i = ish
        txs_j = jsh
        txs_k = ksh
        txs_l = lsh
        num_quart=1
        dum_log=.false.
c
        call texas_hf2_m(
     &      brain,txs_i,txs_j,
     &      ketin,txs_k,txs_l,num_quart,
     &      q4,.false.,
c...............................use roff set false
     &      roff,roff,roff,roff,.false.,
     &      eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
c...............gen labs .. more_integrals
     &      nint_ck, .false., dum_log, scr, lscr, 0.0d0,'der1_int')
        if (nint_ck .eq. 0) then
          nintzero = int_nint(brain,ish,brain,jsh,ketin,ksh,ketin,lsh)
          nintzero = nintzero*12
          call dcopy(nintzero, 0.0d0, 0, eri, 1)
        endif
c
        if (nint_ck*12.gt.lscr) call errquit
     $      ('intd_2e4c: lscr is too small for texas derivatives',
     $      911, INT_ERR)
        scr_ptr = lscr - nint_ck*12 - 1 
        call dcopy((nint_ck*12),eri,1,scr(scr_ptr),1)
        call intd_texas_grad_switch(nint_ck,eri,scr(scr_ptr))
        call intd_sum(eri,nint_ck,idatom,iatom,jatom,katom,latom)
c
      else if (status_nw) then
        if (status_rel) then
          call rel_2e4cd_sf (
     &        coords(1,iatom,ab_geom),dbl_mb(mb_exndcf(iexp,bra)),
     &        dbl_mb(mb_exndcf(icf,bra)),dbl_mb(mb_exndcf(icfs,bras)),
     &        inp,igen,itype,iatom,
c     
     &        coords(1,jatom,ab_geom),dbl_mb(mb_exndcf(jexp,bra)),
     &        dbl_mb(mb_exndcf(jcf,bra)),dbl_mb(mb_exndcf(jcfs,bras)),
     &        jnp,jgen,jtype,jatom,
c
     &        coords(1,katom,cd_geom),dbl_mb(mb_exndcf(kexp,ket)),
     &        dbl_mb(mb_exndcf(kcf,ket)),dbl_mb(mb_exndcf(kcfs,kets)),
     &        knp,kgen,ktype,katom,
c
     &        coords(1,latom,cd_geom),dbl_mb(mb_exndcf(lexp,ket)),
     &        dbl_mb(mb_exndcf(lcf,ket)),dbl_mb(mb_exndcf(lcfs,kets)),
     &        lnp,lgen,ltype,latom,
c
     &        eri,nint_ck,.false.,.false.,.false.,.false.,
     &        scr,lscr,bra_rel,ket_rel,ss_one_cent,do_ssss,rel_dbg)
        else
          call hf2d(
     &      coords(1,iatom,ab_geom),dbl_mb(mb_exndcf(iexp,bra)),
     &      dbl_mb(mb_exndcf(icf,bra)),inp,igen,itype,iatom,
c
     &      coords(1,jatom,ab_geom),dbl_mb(mb_exndcf(jexp,bra)),
     &      dbl_mb(mb_exndcf(jcf,bra)),jnp,jgen,jtype,jatom,
c
     &      coords(1,katom,cd_geom),dbl_mb(mb_exndcf(kexp,ket)),
     &      dbl_mb(mb_exndcf(kcf,ket)),knp,kgen,ktype,katom,
c
     &      coords(1,latom,cd_geom),dbl_mb(mb_exndcf(lexp,ket)),
     &      dbl_mb(mb_exndcf(lcf,ket)),lnp,lgen,ltype,latom,
c
     &      eri,nint_ck,.false.,.false.,.false.,.false.,
     &      scr,lscr)
        end if
c
        call intd_sum(eri,nint_ck,idatom,iatom,jatom,katom,latom)
        call intd_2ec2s(eri,nint_ck,scr,lscr,
     &      itype,jtype,ktype,ltype,igen,jgen,kgen,lgen,
     &      bas_spherical(bra),bas_spherical(ket),idatom,
     T       .false.)
c
      else if (status_sim) then
#ifdef SIMINT_GRADIENT
        call nwcsim_hf2d(
     &        bra,ish,jsh,
     &        ket,ksh,lsh,
     &        nint_ck, eri, leri, scr, lscr)
c
        if(nint_ck.ne.0) then
           nint_ck=        
     i          int_nbf_x(itype)*int_nbf_x(jtype)*
     k          int_nbf_x(ktype)*int_nbf_x(ltype)

c     skip for sp
           if((bas_spherical(bra).or.bas_spherical(ket)).and.
     A      (itype.gt.1.or.jtype.gt.1.or.ktype.gt.1.or.ltype.gt.1))then
#if 1
           call intd_sum_sim(eri,nint_ck,idatom,iatom,jatom,katom,latom)
#else
              idatom(1)=iatom
              idatom(2)=jatom
              idatom(3)=katom
              idatom(4)=latom
#endif
           call intd_2ec2s(eri,nint_ck,scr,lscr,
     &          itype,jtype,ktype,ltype,igen,jgen,kgen,lgen,
     &          bas_spherical(bra),bas_spherical(ket),idatom,
     T          .true.)
           endif
        endif
#else
         call errquit(' simint 2e4c derivatives not ready yet',0,0)
#endif
      else
        write(luout,*)'intd_2e4c: could not use either texas or nwchem'
        write(luout,*)'           integral derivatives'
        write(luout,*)' Please notify nwchem-support@emsl.pnl.gov with'
        write(luout,*)'        the input deck and output available'
        call errquit('intd_2e4c: fatal error',911, INT_ERR)
      endif
c
      end
C>
C> \brief Based on translational invariance combine integral derivatives
C>
      subroutine intd_sum(eri,nint_ck,idatom,iatom,jatom,katom,latom)
      implicit none
c
      integer nint_ck
      double precision eri(nint_ck,3,4)
      integer idatom(4)
      integer iatom,jatom,katom,latom
c
      integer iduse
c
      call intd_logic_atom(idatom,iatom,jatom,katom,latom)
c
      do 00100 iduse = 2,4
        if (idatom(iduse).gt.0) then
          continue
        else
          call daxpy(nint_ck*3,1.0d00,eri(1,1,iduse),1,
     &           eri(1,1,abs(idatom(iduse))),1)
        endif
00100 continue
c
      end
      subroutine intd_sum_sim(eri,nint_ck,
     I     idatom,iatom,jatom,katom,latom)
      implicit none
c
      integer nint_ck
      double precision eri(3,4,nint_ck)
      integer idatom(4)
      integer iatom,jatom,katom,latom
c
      integer iduse,ii
c
      call intd_logic_atom(idatom,iatom,jatom,katom,latom)
c
      do iduse = 2,4
        if (idatom(iduse).gt.0) then
          continue
        else
           do ii=1,nint_ck
c          call yaxpy(nint_ck*3,1.0d00,eri(1,1,iduse),1,
c     &           eri(1,1,abs(idatom(iduse))),1)
              eri(1,abs(idatom(iduse)),ii)=
     E             eri(1,abs(idatom(iduse)),ii)+
     +             eri(1,iduse,ii)
              eri(1,iduse,ii)=0d0
              eri(2,abs(idatom(iduse)),ii)=
     E             eri(2,abs(idatom(iduse)),ii)+
     +             eri(2,iduse,ii)
              eri(2,iduse,ii)=0d0
              eri(3,abs(idatom(iduse)),ii)=
     E             eri(3,abs(idatom(iduse)),ii)+
     +             eri(3,iduse,ii)
              eri(3,iduse,ii)=0d0
           enddo
        endif
      enddo
c
      end
C>
C> \brief Work out aspects of the translational invariance
C>
      subroutine intd_logic_atom(idat,iat,jat,kat,lat)
      implicit none
      integer iat,jat,kat,lat
      integer idat(4)
c
      idat(1) = iat
      idat(2) = jat
      idat(3) = kat
      idat(4) = lat
      if (iat.eq.jat) idat(2) = -1
      if (iat.eq.kat) idat(3) = -1
      if (iat.eq.lat) idat(4) = -1
      if (jat.eq.kat) then
        if(idat(2).gt.0) then
          idat(3) = -2
        else
          idat(3) = idat(2)
        endif
      endif
      if (jat.eq.lat) then
        if(idat(2).gt.0) then
          idat(4) = -2
        else
          idat(4) = idat(2)
        endif
      endif
      if (kat.eq.lat) then
        if(idat(3).gt.0) then
          idat(4) = -3
        else
          idat(4) = idat(3)
        endif
      endif
      end
C>
C> \brief Transform 4-center 2-electron integral derivatives from 
C> Cartesian to spherical harmonic basis functions
C>
      subroutine intd_2ec2s(eri,nint_ck,scr,lscr,
     &    it,jt,kt,lt,igin,jgin,kgin,lgin,
     &    bra_sph,ket_sph,idatom,ltransp)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "int_nbf.fh"
c::passed
      integer nint_ck
      integer lscr
      double precision eri(nint_ck,3,4)
      double precision scr(lscr)
      integer it,jt,kt,lt
      integer igin,jgin,kgin,lgin
      logical bra_sph, ket_sph
      integer idatom(4)
c::local      
      integer ig,jg,kg,lg
      integer nint_x, nint_s
      integer i_nbf,j_nbf,k_nbf,l_nbf
      integer i_nbf_s,j_nbf_s,k_nbf_s,l_nbf_s
      integer zatom
      integer zyx, lda, ldb
      logical oprint,ltransp
c
      oprint=.false.
c
      ig = igin
      jg = jgin
      kg = kgin
      lg = lgin
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
      if (it.eq.-1) ig = 1
      if (jt.eq.-1) jg = 1
      if (kt.eq.-1) kg = 1
      if (lt.eq.-1) lg = 1
      if (bra_sph.or.ket_sph) then
        i_nbf = int_nbf_x(It)
        j_nbf = int_nbf_x(Jt)
        k_nbf = int_nbf_x(Kt)
        l_nbf = int_nbf_x(Lt)
        nint_x = i_nbf*j_nbf*k_nbf*l_nbf
C should be
C        nint_x = i_nbf*j_nbf*k_nbf*l_nbf*ig*jg*kg*lg
C ?
        if (nint_ck.ne.nint_x) call errquit
     &      ('intd_2ec2s: nint_ck.ne.nint_x diff=',(nint_ck-nint_x),
     &       INT_ERR)
c
c     transpose integrals  X(12,nint_x) --> Y(nint_x,12)
c
        if(ltransp) then
           lda=12
           ldb=nint_x
           if(lscr.lt.nint_x) call errquit(' lscr small ',0,0)
           call dcopy((nint_x*lda),eri,1,scr,1)
           call trspmo_block(scr,lda,  eri,ldb)
        endif
        if(bra_sph) then
          i_nbf_s = int_nbf_s(It)
          j_nbf_s = int_nbf_s(Jt)
          do zatom = 1,4
            if (idatom(zatom).gt.0) then
              do zyx = 1,3
                call spcart_bra2etran(eri(1,zyx,zatom),scr,
     &              j_nbf,i_nbf,j_nbf_s,i_nbf_s,
     &              Jt, It, jg, ig,
     &              (k_nbf*l_nbf),oprint)
C should be
C     &              (k_nbf*l_nbf*kg*lg),.true.)
C ?
              enddo
            endif
          enddo
          i_nbf = i_nbf_s
          j_nbf = j_nbf_s
        endif
        if(ket_sph) then
          k_nbf_s = int_nbf_s(Kt)
          l_nbf_s = int_nbf_s(Lt)
          do zatom = 1,4
            if (idatom(zatom).gt.0) then
              do zyx = 1,3
                call spcart_ket2etran(eri(1,zyx,zatom),scr,
     &              l_nbf,k_nbf,l_nbf_s,k_nbf_s,
     &              Lt, Kt, lg, kg,
     &              (i_nbf*j_nbf),oprint)
C should be
C     &              (i_nbf*j_nbf*ig*jg),.true.)
C ?
              enddo
            endif
          enddo
          k_nbf = k_nbf_s
          l_nbf = l_nbf_s
        endif
        nint_s = i_nbf*j_nbf*k_nbf*l_nbf
        if(ltransp) then
c
c     transpose back integrals  X(12,nint_x) --> Y(nint_x,12)
c
           call dcopy((nint_x*lda),eri,1,scr,1)
           call trspmo_block(scr,ldb,  eri,lda)
           return
        endif
        if (nint_s.gt.nint_x) then
          call errquit
     &      ('intd_2ec2s: nint_s >.nint_x diff=',(nint_s-nint_x),
     &       INT_ERR)
        elseif (nint_s.eq.nint_x) then
          return
        else
           call int_c2s_mv
     &          (eri,nint_x,nint_s,(3*4),scr,lscr,'intd_2e4c')
        endif
      endif
      end
      subroutine int_c2s_mv
     &    (int_buf,nint_x,nint_s,nblocks,scr,lscr,ctine)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c::passed
      integer lscr
      integer nint_x
      integer nint_s
      integer nblocks
      double precision int_buf(*)
      double precision scr(lscr)
      character*(*) ctine
c::local
      integer i
      integer z
      integer offset_x, offset_s
*
#if defined(VECTOR_MODE)
      if (nint_x*nblocks.gt.lscr) then
        write(luout,*)' calling routine: ',ctine
        call errquit
     &      ('int_c2s_mv: lscr to small by ',((nint_x*nblocks)-lscr),
     &         INT_ERR)
      endif
#endif
      if (nint_s.gt.nint_x) then
        call errquit
     &      ('int_c2s_mv: nint_s >.nint_x diff=',(nint_s-nint_x),
     &        INT_ERR)
      elseif (nint_s.eq.nint_x) then
        return
      endif
#if defined(VECTOR_MODE)
      call dcopy((nint_x*nblocks),int_buf,1,scr,1)
      call dcopy((nint_x*nblocks),0.0d00,0,int_buf,1)
      do z = 1,nblocks
        offset_x = (z-1)*nint_x + 1
        offset_s = (z-1)*nint_s + 1
        call dcopy(nint_s,scr(offset_x),1,int_buf(offset_s),1)
      enddo
#else
c** scalar  (in place)
      do z = 2,nblocks
        offset_x = (z-1)*nint_x 
        offset_s = (z-1)*nint_s 
        do i = 1,nint_s
          int_buf(offset_s+i) = int_buf(offset_x+i)
        enddo
      enddo
#endif
      end
c=====================================================
C>
C> \brief Reorder integrals derivatives from Texas to NWChem ordering
C>
      subroutine intd_texas_grad_switch(nint_ck,eri,scr)
      implicit none
      integer nint_ck
      double precision eri(nint_ck,3,4)
      double precision scr(12,nint_ck)
c
      double precision xa, xb, xc, xd
      double precision ya, yb, yc, yd
      double precision za, zb, zc, zd
      integer i_int
c
      do i_int = 1,nint_ck
         xa = scr(1,i_int)
         ya = scr(2,i_int)
         za = scr(3,i_int)
         xb = scr(4,i_int)
         yb = scr(5,i_int)
         zb = scr(6,i_int)
         xc = scr(7,i_int)
         yc = scr(8,i_int)
         zc = scr(9,i_int)
         xd = scr(10,i_int)
         yd = scr(11,i_int)
         zd = scr(12,i_int)
         eri(i_int,1,1) = xa
         eri(i_int,2,1) = ya
         eri(i_int,3,1) = za
         eri(i_int,1,2) = xb
         eri(i_int,2,2) = yb
         eri(i_int,3,2) = zb
         eri(i_int,1,3) = xc
         eri(i_int,2,3) = yc
         eri(i_int,3,3) = zc
         eri(i_int,1,4) = xd 
         eri(i_int,2,4) = yd 
         eri(i_int,3,4) = zd 
      enddo
      end
      subroutine trspmo_block(in,ld_in,out,ld_out)
#define CHUNK 16
      implicit none
      integer ld_in,ld_out
      double precision out(ld_out,ld_in)
      double precision in(ld_in,ld_out)
c
      integer i,j,i1,j1
c
      do i1=1,ld_in,CHUNK
         do j1=1,ld_out,CHUNK
            do i=i1,min(i1+CHUNK-1,ld_in)
!DEC$ LOOP COUNT AVG=CHUNK
cc!deC$ SIMD
               do j=j1,min(j1+CHUNK-1,ld_out)
                  out(j,i)=in(i,j)
               enddo
            enddo
         enddo
      enddo
      return
      end
c=====================================================
C> @}
