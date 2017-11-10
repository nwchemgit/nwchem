c $Id$
*
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_2e4c}
c:tex-this routine computes the 4 center (traditional) 2 electron integrals: 
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
C> \ingroup nwint
C> @{
c
C> \brief Calculate a block of 4 center 2-electron integrals
c
C> The routine takes 4 shells of basis functions. The first pair of
C> shells generates the "bra" charge distribution and the second pair
C> the "ket" charge distribution. This routine evaluates the integral
C> of the Coulomb operator for the repulsion of these charge
C> distributions. In summary the following expression is evaluated:
C> \f{eqnarray*}{
C>    ({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty}
C>    g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}
C>    g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C>
C> The order of the integrals returned can be specified considering
C> that the rank of the basis functions in the shells are given by:
C> - `ibflo` is the lowest rank in `ish`, and `ibfhi` is the highest
C> - `jbflo` is the lowest rank in `jsh`, and `jbfhi` is the highest
C> - `kbflo` is the lowest rank in `ksh`, and `kbfhi` is the highest
C> - `lbflo` is the lowest rank in `lsh`, and `lbfhi` is the highest
C>
C> The number of basis functions in each shell is given by:
C> - `inbf = ibfhi - ibflo + 1`
C> - `jnbf = jbfhi - jbflo + 1`
C> - `knbf = kbfhi - kbflo + 1`
C> - `lnbf = lbfhi - lbflo + 1`
C>
C> and further
C> - `klnbf = knbf * lnbf`
C> - `jklnbf = jnbf * klnbf`
C> - `ijklnbf = ijnbf * jklnbf`
C>
C> then 
C> \f$(ij|kl) = \mathrm{eri}((i-ibflo)*jklnbf+(j-jbflo)*klnbf+
C> (k-kbflo)*lnbf+(l-lbflo)+1)\f$
C>
C> At present this routine provides the following algorithms:
C>
C> - SP-rotated axis 
C>
C> - McMurchie-Davidson
C>
C> - Obara-Saika (Texas)
C>
C> Not available are:
C>
C> - Rys quadrature (Hondo)
C>
C> The choice of algorithm is controlled by capability checks and user
C> inputs (see the various cando functions for details).
c
      subroutine int_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri)
c:tex-\end{verbatim}
      implicit none
c
c basic api routine to generate a block of two electron integrals
c eri = <bra_g(ish).bra_g(jsh) | ket_g(ksh).ket_g(lsh)>
c
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "apiP.fh"
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
#include "global.fh"
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical cando_sp
      logical cando_txs
      logical cando_sim
      external cando_nw
      external cando_sp
      external cando_txs
      external cando_sim
#define PERF_NORMAL
#if defined(PERF_NORMAL)
      logical int_chk_sh
      logical int_chk_init
      external int_chk_sh
      external int_chk_init
#endif
c:: passed
c:tex-\begin{verbatim}
      integer brain              !< [Input] the bra basis set handle
      integer ish                !< [Input] the shell index
      integer jsh                !< [Input] the shell index
      integer ketin              !< [Input] the ket basis set handle
      integer ksh                !< [Input] the shell index
      integer lsh                !< [Input] the shell index
      integer lscr               !< [Input] the length of scratch array
      double precision scr(lscr) !< [Scratch] array
      integer leri               !< [Input] the length of integral array
      double precision eri(leri) !< [Output] the 2-electron integrals
c:tex-\end{verbatim}
c:: local
      integer bra, ket 
      integer ab_geom, cd_geom, ucont
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_cent, a_icfps
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_cent, b_icfps
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent, c_icfps
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent, d_icfps
c.txs
      integer num_quart, dummy_lab(2)
      integer txs_i,txs_j,txs_k,txs_l
      double precision roff(3)
      double precision q4
      integer nint
      logical dum_log
      logical status_sp, status_nw, status_txs, status_gen,
     L     status_sim
      integer texas_ang_limit
c.rel-dmd
      logical status_rel, bra_rel, ket_rel
      logical i_rel, j_rel, k_rel, l_rel
      integer sbas, abas, bras, kets
c
      logical any_spherical
      integer a_nbf,   b_nbf,   c_nbf,   d_nbf
      integer a_nbf_s, b_nbf_s, c_nbf_s, d_nbf_s
      integer ab_gen, ab_cmp, cd_gen, cd_cmp
c
#if defined(PERF_NORMAL)
      logical shells_ok
      logical inline_chk_sh
c
      integer WarnP
      save WarnP
      data WarnP /0/
c... statement function for int_chk_sh
      inline_chk_sh(bra,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(bra)))
#endif
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
#if defined(PERF_NORMAL)
c
c check initialization
c
      if (.not.int_chk_init('int_2e4c'))
     &       call errquit('int_2e4c: int_init was not called' ,0,
     &       INT_ERR)
c
c     define center information required 
c
      if (.not. bas_check_handle(brain,'int_2e4c:bra')) call errquit
     $     ('int_2e4c: invalid bra basis handle', brain, BASIS_ERR)
      if (.not. bas_check_handle(ketin,'int_2e4c:ket')) call errquit
     $     ('int_2e4c: invalid ket basis handle', ketin, BASIS_ERR)
#endif
*
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      bras = bra
      kets = ket
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
#if defined(PERF_NORMAL)
c
c     check input shell ranges
c
      shells_ok = inline_chk_sh(bra,ish)
      shells_ok = shells_ok .and. inline_chk_sh(bra,jsh)
      shells_ok = shells_ok .and. inline_chk_sh(ket,ksh)
      shells_ok = shells_ok .and. inline_chk_sh(ket,lsh)
      if (.not.shells_ok)
     &       call errquit('int_2e4c: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      if (ab_geom.ne.cd_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_2e4c: WARNING: possible geometry inconsistency'
        write(luout,*)'bra geometry handle:',ab_geom
        write(luout,*)'ket geometry handle:',cd_geom
        WarnP = 1
      endif
c
#endif
      a_cent  = (sf_ibs_cn2ce(ish,bra))
      b_cent  = (sf_ibs_cn2ce(jsh,bra))
      c_cent  = (sf_ibs_cn2ce(ksh,ket))
      d_cent  = (sf_ibs_cn2ce(lsh,ket))
c
c     used ony for nw and sim
      any_spherical = .false.
c
c
c... stat sp
      status_sp = cando_sp(brain,ish,jsh).and.cando_sp(ketin,ksh,lsh)

c... stat texas
      texas_ang_limit = 8   ! (ff|pp) or greater integral
      ucont   = (sf_ibs_cn2ucn(ish,bra))
      La      = infbs_cont(CONT_TYPE ,ucont,bra)
      a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
      ucont   = (sf_ibs_cn2ucn(jsh,bra))
      Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
      b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
      ucont   = (sf_ibs_cn2ucn(ksh,ket))
      Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
      c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
      ucont   = (sf_ibs_cn2ucn(lsh,ket))
      Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
      d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
      status_txs = cando_txs(brain,ish,jsh).and.cando_txs(ketin,ksh,lsh)
      status_gen = max(a_gen,b_gen,c_gen,d_gen).gt.1
      if (.not.status_gen) then
        status_txs = status_txs .and.
     &        ( (La+Lb+Lc+Ld).gt.texas_ang_limit  )
      endif
c... stat nw
      status_nw  = cando_nw(brain,ish,jsh).and.cando_nw(ketin,ksh,lsh)
c... stat rel
      status_rel = dyall_mod_dir .and. .not.nesc_1e_approx
     &    .and. (brain .eq. ketin) .and. (brain .eq. ao_bsh)
c... stat sim
      status_sim = cando_sim(brain,ish,jsh).and.cando_sim(ketin,ksh,lsh)
c
c  The following 4 lines are to take care of compiler warnings.
c
      i_rel = .false.
      j_rel = .false.
      k_rel = .false.
      l_rel = .false.
      if (status_rel) then
c
c     get basis set handles; relativistic integral option valid
c     if bra or ket are the ao basis and bra and ket have both
c     functions relativistic
c
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
      if (status_sp) then
c
        call genr70(
     &         brain,ish,coords(1,a_cent,ab_geom),
     &               jsh,coords(1,b_cent,ab_geom),
     &         ketin,ksh,coords(1,c_cent,cd_geom),
     &               lsh,coords(1,d_cent,cd_geom),
     &         eri)
c        
      elseif(status_txs .and. .not.status_rel)  then
        num_quart = 1
        dummy_lab(1) = 0
        dummy_lab(2) = 0
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        dum_log = .false.
*...check symmetry here
        q4 = 1.0d00
        txs_i = ish
        txs_j = jsh
        txs_k = ksh
        txs_l = lsh
        call texas_hf2_m(
     &        brain,txs_i,txs_j,
     &        ketin,txs_k,txs_l,num_quart,
     &        q4,.false.,
c............................... use roff set false
     &        roff,roff,roff,roff,.false.,
     &        eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
c............... gen labs .. more_integrals
     &        nint, .false., dum_log, scr, lscr, 0.0d0,'scfd_int')
        if (nint .eq. 0) then
           a_nbf = sf_ibs_cn2bfr(2,ish,bra) - sf_ibs_cn2bfr(1,ish,bra)+1
           b_nbf = sf_ibs_cn2bfr(2,jsh,bra) - sf_ibs_cn2bfr(1,jsh,bra)+1
           c_nbf = sf_ibs_cn2bfr(2,ksh,ket) - sf_ibs_cn2bfr(1,ksh,ket)+1
           d_nbf = sf_ibs_cn2bfr(2,lsh,ket) - sf_ibs_cn2bfr(1,lsh,ket)+1
           call dcopy(a_nbf*b_nbf*c_nbf*d_nbf, 0.0d0, 0, eri, 1)
        endif

      elseif(status_nw) then
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        La      = infbs_cont(CONT_TYPE ,ucont,bra)
        a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
        if (i_rel) ucont = ao_to_ls(ucont)
        a_icfps = infbs_cont(CONT_ICFP ,ucont,bras)
c
        ucont   = (sf_ibs_cn2ucn(jsh,bra))
        Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
        b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
        if (j_rel) ucont = ao_to_ls(ucont)
        b_icfps = infbs_cont(CONT_ICFP ,ucont,bras)
c      
        ucont   = (sf_ibs_cn2ucn(ksh,ket))
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        if (k_rel) ucont = ao_to_ls(ucont)
        c_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        ucont   = (sf_ibs_cn2ucn(lsh,ket))
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        if (l_rel) ucont = ao_to_ls(ucont)
        d_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        if (status_rel) then
          call rel_2e4c_sf (
     &        coords(1,a_cent,ab_geom), dbl_mb(mb_exndcf(a_iexp,bra)),
     &        dbl_mb(mb_exndcf(a_icfp,bra)),
     &        dbl_mb(mb_exndcf(a_icfps,bras)), a_prim,a_gen,La,a_cent,
     &        coords(1,b_cent,ab_geom), dbl_mb(mb_exndcf(b_iexp,bra)),
     &        dbl_mb(mb_exndcf(b_icfp,bra)), 
     &        dbl_mb(mb_exndcf(b_icfps,bras)), b_prim,b_gen,Lb,b_cent,
     &        coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &        dbl_mb(mb_exndcf(c_icfp,ket)), 
     &        dbl_mb(mb_exndcf(c_icfps,kets)), c_prim,c_gen,Lc,c_cent,
     &        coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &        dbl_mb(mb_exndcf(d_icfp,ket)), 
     &        dbl_mb(mb_exndcf(d_icfps,kets)), d_prim,d_gen,Ld,d_cent,
c...................... canAB   canCD   canPQ   DryRun
     &        eri,leri,.false.,.false.,.false.,.false.,scr,lscr,
     &        bra_rel,ket_rel,ss_one_cent,do_ssss,rel_dbg)
        else
          call hf2(
     &        coords(1,a_cent,ab_geom), dbl_mb(mb_exndcf(a_iexp,bra)),
     &        dbl_mb(mb_exndcf(a_icfp,bra)), a_prim, a_gen, La,
     &        coords(1,b_cent,ab_geom), dbl_mb(mb_exndcf(b_iexp,bra)),
     &        dbl_mb(mb_exndcf(b_icfp,bra)), b_prim, b_gen, Lb,
     &        coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &        dbl_mb(mb_exndcf(c_icfp,ket)), c_prim, c_gen, Lc,
     &        coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &        dbl_mb(mb_exndcf(d_icfp,ket)), d_prim,d_gen,Ld,
c....................... canAB    canCD   canPQ   dryrun
     &        eri, leri, .false.,.false.,.false.,.false., scr, lscr)
        end if
      elseif(status_sim .and. .not.status_rel)  then
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)

        call nwcsim_hf2(
     &        bra,ish,jsh,
     &        ket,ksh,lsh,
     &        nint, eri, leri, scr, lscr)
        if (nint .eq. 0) then
           a_nbf = sf_ibs_cn2bfr(2,ish,bra) - sf_ibs_cn2bfr(1,ish,bra)+1
           b_nbf = sf_ibs_cn2bfr(2,jsh,bra) - sf_ibs_cn2bfr(1,jsh,bra)+1
           c_nbf = sf_ibs_cn2bfr(2,ksh,ket) - sf_ibs_cn2bfr(1,ksh,ket)+1
           d_nbf = sf_ibs_cn2bfr(2,lsh,ket) - sf_ibs_cn2bfr(1,lsh,ket)+1
           call dcopy(a_nbf*b_nbf*c_nbf*d_nbf, 0.0d0, 0, eri, 1)
        endif

      else
        write(luout,*)'int_2e4c: could not do sp or nw integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        write(luout,*)' ksh   :',ksh
        write(luout,*)' lsh   :',lsh
        call errquit('int_2e4c: fatal error ',0, INT_ERR)
      endif
c
c
c eri has cartesian block of integrals (llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
c
        if (any_spherical)then
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
c
          if (la.eq.-1) a_gen = 1
          if (lb.eq.-1) b_gen = 1
          if (lc.eq.-1) c_gen = 1
          if (ld.eq.-1) d_gen = 1
          
          a_nbf = int_nbf_x(La)
          b_nbf = int_nbf_x(Lb)
          c_nbf = int_nbf_x(Lc)
          d_nbf = int_nbf_x(Ld)
          cd_gen = c_gen*d_gen
          ab_gen = a_gen*b_gen
          ab_cmp = a_nbf*b_nbf
          cd_cmp = c_nbf*d_nbf
          if(bas_spherical(bra)) then
            a_nbf_s = int_nbf_s(La)
            b_nbf_s = int_nbf_s(Lb)
            call spcart_bra2etran(eri,scr,
     &            b_nbf,a_nbf,b_nbf_s,a_nbf_s,
     &            Lb, La, b_gen, a_gen,
     &            cd_gen*cd_cmp,.false.)
            a_nbf = a_nbf_s
            b_nbf = b_nbf_s
            ab_cmp = a_nbf*b_nbf
          endif
          if(bas_spherical(ket)) then
            c_nbf_s = int_nbf_s(Lc)
            d_nbf_s = int_nbf_s(Ld)
            call spcart_ket2etran(eri,scr,
     &            d_nbf,c_nbf,d_nbf_s,c_nbf_s,
     &            Ld, Lc, d_gen, c_gen,
     &            ab_gen*ab_cmp,.false.)
            c_nbf = c_nbf_s
            d_nbf = d_nbf_s
          endif
        endif
      end
C> @}
