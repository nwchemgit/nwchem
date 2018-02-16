c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Computes 3 center 2-electron integrals
C>
C> Computes 3 center 2-electron integrals of the following kind:
C> \f{eqnarray*}{
C> ({\mu}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}
C> g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_2e3c}
c:tex-this routine computes the 3 center 2 electron integrals: 
c:tex-\begin{eqnarray*}
c:tex-({\mu}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}
c:tex-g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_2e3c(brain, ish, ketin, jsh, ksh,
     &       lscr, scr, leri, eri)
c:tex-\end{verbatim}
      implicit none
c
c basic api routine to generate a block of 3 center two electron integrals
c eri = <bra_g(ish)|ket_g(jsh).ket_g(ksh)>
c
#include "apiP.fh"
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
#include "int_nbf.fh"
#include "stdio.fh"
#include "rel_nwc.fh"
#include "util.fh"
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical cando_sp
      logical cando_sim
      logical int_chk_sh
      logical int_chk_init
      external cando_nw
      external cando_sp
      external cando_sim
      external int_chk_sh
      external int_chk_init
*----------------------------#define USE_TEXAS
cedo#if defined(USE_TEXAS_BROKE)
      logical cando_txs
      external cando_txs
cedo#endif
c::functions
      integer int_nint_cart
      external int_nint_cart
c:: passed
c:tex-\begin{verbatim}
      integer brain !< [Input] bra basis set handle
      integer ish   !< [Input] shell/contraction index
      integer ketin !< [Input] ket basis set handle
      integer jsh   !< [Input] shell/contraction index
      integer ksh   !< [Input] shell/contraction index
      integer lscr  !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] array
      integer leri  !< [Input] length of integral array
      double precision eri(leri) !< [Output] 2e3c integrals
c:tex-\end{verbatim}
c:: local
      logical shells_ok
      integer bra, ket 
      integer p_geom, cd_geom, ucont
      integer Lp, p_prim, p_gen, p_iexp, p_icfp, p_cent
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent
      logical status_nw, status_sp, status_sim
cedo#if defined(USE_TEXAS_BROKE)
      logical status_gen
      logical status_txs
      integer num_quart, dummy_lab(2)
      double precision roff(3), q4 
      integer txs_i,txs_j,txs_k,txs_d
      logical dum_log
cedo#endif
      integer nint
      logical OFALSE
      logical any_spherical
      logical used_txs
c.rel-dmd
      logical status_rel, ket_rel, j_rel, k_rel
      integer sbas, abas, kets, c_icfps, d_icfps
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      OFALSE = .false.
      used_txs=.false.
c
c check if ERI is big enough - TLW
c
      nint = int_nint_cart(brain,ish,ketin,jsh,ketin,ksh,ketin,0)
      if (nint.gt.leri) then
        write(luout,*) 'nint = ',nint
        write(luout,*) 'leri     = ',leri
        call errquit('int_2e3c: nint>leri error',911, INT_ERR)
      endiF
c
c check initialization
c
      if (.not.int_chk_init('int_2e3c'))
     &       call errquit('int_2e3c: int_init was not called' ,0,
     &       INT_ERR)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      if (.not.shells_ok)
     &       call errquit('int_2e3c: invalid contraction/shell',0,
     &       BASIS_ERR)
cedo#ifndef USE_TEXAS_BROKE
      call int_nogencont_check(brain,'int_2e3c:bra')
      call int_nogencont_check(ketin,'int_2e3c:ket')
chvd  I currently find no evidence that the code cannot handle SP-shells
c     call int_nospshell_check(brain,'int_2e3c:bra')
c     call int_nospshell_check(ketin,'int_2e3c:ket')
cedo#endif
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      kets = ket
      p_geom  = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
c
c  check if spherical
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
c
      if (p_geom.ne.cd_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_2e3c: WARNING: possible geometry inconsistency'
        write(luout,*)'bra geometry handle:',p_geom
        write(luout,*)'ket geometry handle:',cd_geom
        WarnP = 1
      endif
c
      p_cent  = (sf_ibs_cn2ce(ish,bra))
      c_cent  = (sf_ibs_cn2ce(jsh,ket))
      d_cent  = (sf_ibs_cn2ce(ksh,ket))
c
      ucont   = (sf_ibs_cn2ucn(ish,bra))
      Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
      p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
c
      ucont   = (sf_ibs_cn2ucn(jsh,ket))
      Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
      c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
c
      ucont   = (sf_ibs_cn2ucn(ksh,ket))
      Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
      d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)

*
c
c set integral code status
c
      status_nw = cando_nw(brain,ish,0).and.cando_nw(ketin,jsh,ksh)
      status_sp = cando_sp(brain,ish,0).and.cando_sp(ketin,jsh,ksh)
      status_sim = cando_sim(brain,ish,0).and.cando_sim(ketin,jsh,ksh)
cedo#if defined(USE_TEXAS_BROKE)
      status_txs = cando_txs(brain,ish,0).and.cando_txs(ketin,jsh,ksh)
      status_gen = (max(p_gen,c_gen,d_gen)) .gt. 1  ! if general contraction texas is only option

*                                (p|ff) or greater do texas
      if (.not.status_gen)
     &    status_txs = status_txs .and. ((abs(Lp)+abs(Lc)+abs(Ld)).ge.7)
cedo#endif
*
      status_rel = dyall_mod_dir .and. .not.nesc_1e_approx
     &    .and. (ketin .eq. ao_bsh)
      k_rel=.false.
      j_rel=.false.
      if (status_rel) then
c
c     get basis set handles; relativistic integral option valid
c     if bra or ket are the ao basis and bra and ket have both
c     functions relativistic
c
        ket_rel = .false.
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        kets = sbas
        ket_rel = ket .eq. abas
        if (ket_rel) then
          ucont = sf_ibs_cn2ucn(jsh,ket)
          j_rel = infbs_cont(CONT_RELLS ,ucont,ket) .ne. 0
          ucont = sf_ibs_cn2ucn(ksh,ket)
          k_rel = infbs_cont(CONT_RELLS ,ucont,ket) .ne. 0
          ket_rel = ket_rel .and. j_rel .and. k_rel
        end if
        status_rel = status_rel .and. ket_rel
      end if
c
      if (status_sp .and. .not.status_rel) then
        call genr70(
     &         brain,ish,coords(1,p_cent,p_geom),
     &                 0,coords(1,p_cent,p_geom),
     &         ketin,jsh,coords(1,c_cent,cd_geom),
     &               ksh,coords(1,d_cent,cd_geom),
     &         eri)
cedo#if defined(USE_TEXAS_BROKE)
      else if (status_txs .and. .not.status_rel) then
        num_quart = 1
        dummy_lab(1) = 0
        dummy_lab(2) = 0
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        txs_i = ish
        txs_j = jsh
        txs_k = ksh
        txs_d = 0
        dum_log = OFALSE
        q4 = 1.0d00
        nint = 0
        call texas_hf2_m(
     &        brain,txs_i,txs_d,
     &        ketin,txs_j,txs_k,num_quart,
     &        q4,OFALSE,roff,roff,roff,roff,OFALSE,
     &        eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
     &        nint, OFALSE, dum_log, scr, lscr, 0.0d0,'scfd_int')
        used_txs=.true.
cedo#endif
      else if(status_nw) then
c
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
        p_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        p_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        p_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(jsh,ket))
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        if (j_rel) ucont = ao_to_ls(ucont)
        c_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        ucont   = (sf_ibs_cn2ucn(ksh,ket))
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        if (k_rel) ucont = ao_to_ls(ucont)
        d_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c      
        if (status_rel) then
          call rel_2e4c_sf (
     &        coords(1,p_cent,p_geom), dbl_mb(mb_exndcf(p_iexp,bra)),
     &        dbl_mb(mb_exndcf(p_icfp,bra)), 
     &        dbl_mb(mb_exndcf(p_icfp,bra)),p_prim,p_gen,Lp,p_cent,
     &        coords(1,p_cent,p_geom),DCexp,DCcoeff,
     &        DCcoeff,1,1,0,p_cent,
     &        coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &        dbl_mb(mb_exndcf(c_icfp,ket)), 
     &        dbl_mb(mb_exndcf(c_icfps,kets)), c_prim,c_gen,Lc,c_cent,
     &        coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &        dbl_mb(mb_exndcf(d_icfp,ket)), 
     &        dbl_mb(mb_exndcf(d_icfps,kets)), d_prim,d_gen,Ld,d_cent,
c...................... canAB   canCD   canPQ   DryRun
     &        eri,leri,.false.,.false.,.false.,.false.,scr,lscr,
c............  bra_rel                    do_ssss
     &        .false.,ket_rel,ss_one_cent,.false.,rel_dbg)
        else
          call hf2(
     &        coords(1,p_cent,p_geom), dbl_mb(mb_exndcf(p_iexp,bra)),
     &        dbl_mb(mb_exndcf(p_icfp,bra)), p_prim, p_gen, Lp,
     &        coords(1,p_cent,p_geom), DCexp, 
     &        DCcoeff, 1, 1, 0,
     &        coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &        dbl_mb(mb_exndcf(c_icfp,ket)), c_prim, c_gen, Lc,
     &        coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &        dbl_mb(mb_exndcf(d_icfp,ket)), d_prim,d_gen,Ld,
c......................... canAB    canCD    canPQ
     &        eri, leri, OFALSE, OFALSE, OFALSE,
c............ dryrun
     &        OFALSE, scr, lscr)
        end if
c
      elseif(status_sim .and. .not.status_rel)  then
        call nwcsim_hf2_3c(
     &        bra,ish,
     &        ket,jsh,ksh,
     &        nint, eri, leri, scr, lscr)
      else
        write(luout,*)'int_2e3c: could not do nw integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        write(luout,*)' ksh   :',ksh
        call errquit('int_2e3c: fatal error ',0, INT_ERR)
      endif
      if (any_spherical.and.(.not.used_txs))then
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
          if (Lp.eq.-1) p_gen = 1
          if (Lc.eq.-1) c_gen = 1 
          if (Ld.eq.-1) d_gen = 1 
          call spcart_3ctran(eri,scr,lscr,
     &        int_nbf_x(Lp),int_nbf_s(Lp),Lp,p_gen,bas_spherical(bra),
     &        int_nbf_x(Lc),int_nbf_s(Lc),Lc,c_gen,bas_spherical(ket),
     &        int_nbf_x(Ld),int_nbf_s(Ld),Ld,d_gen,bas_spherical(ket),
     &        OFALSE)
        endif
      end
C> @}
