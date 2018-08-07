c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute 2 center 2-electron integrals.
C>
C> Computes 2 center 2-electron integrals of the form:
C> \f{eqnarray*}{
C> ({\mu}|{\nu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})dr_{1}dr_{2}
C> \f}
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_2e2c}
c:tex-this routine computes the 2 center 2 electron integrals: 
c:tex-\begin{eqnarray*}
c:tex-({\mu}|{\nu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_2e2c(brain, ish, ketin, jsh, 
     &       lscr, scr, leri, eri)
c:tex-\end{verbatim}
      implicit none
c
c basic api routine to generate a block of 2 center two electron integrals
c eri = <bra_g(ish)|ket_g(jsh)>
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
#define USE_TEXAS
cedo#if defined(USE_TEXAS_BROKE)
      logical cando_txs
      external cando_txs
cedo#endif
c:: passed
c:tex-\begin{verbatim}
      integer brain !< [Input] bra basis set handle
      integer ish   !< [Input] shell/contraction index
      integer ketin !< [Input] ket basis set handle
      integer jsh   !< [Input] shell/contraction index
      integer lscr  !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] array
      integer leri  !< [Input] length of integral array
      double precision eri(leri) !< [Output] 2e2c integrals
c:tex-\end{verbatim}
c:: local
      logical shells_ok
      integer bra, ket 
      integer p_geom, q_geom, ucont
      integer Lp, p_prim, p_gen, p_iexp, p_icfp, p_cent
      integer Lq, q_prim, q_gen, q_iexp, q_icfp, q_cent
      logical status_nw, status_sp,status_sim
cedo#if defined(USE_TEXAS_BROKE)
      logical status_gen
      logical status_txs
      integer num_quart, dummy_lab(2)
      double precision roff(3), q4 
      integer txs_i,txs_j,txs_d1,txs_d2
      integer nint
      logical dum_log
cedo#endif
      logical any_spherical
      logical used_txs
      logical OFALSE
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c
c
      OFALSE = .false.
      used_txs=.false.
c
c check initialization
c
      if (.not.int_chk_init('int_2e2c'))
     &       call errquit('int_2e2c: int_init was not called' ,0,0)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
c
      if (.not.shells_ok)
     &       call errquit('int_2e2c: invalid contraction/shell',0,0)
cedo#ifndef USE_TEXAS_BROKE
      call int_nogencont_check(brain,'int_2e2c:bra')
      call int_nogencont_check(ketin,'int_2e2c:ket')
chvd  I currently find no evidence that the code cannot handle SP-shells
c     call int_nospshell_check(brain,'int_2e2c:bra')
c     call int_nospshell_check(ketin,'int_2e2c:ket')
cedo#endif
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      p_geom = ibs_geom(bra)
      q_geom = ibs_geom(ket)
c
c  check if spherical
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
c
      if (p_geom.ne.q_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_2e2c: WARNING: possible geometry inconsistency'
        write(luout,*)'bra geometry handle:',p_geom
        write(luout,*)'ket geometry handle:',q_geom
        WarnP = 1
      endif
c
      p_cent  = (sf_ibs_cn2ce(ish,bra))
      q_cent  = (sf_ibs_cn2ce(jsh,ket))
c
      ucont   = (sf_ibs_cn2ucn(ish,bra))
      Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
      p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
c
      ucont   = (sf_ibs_cn2ucn(jsh,ket))
      Lq      = infbs_cont(CONT_TYPE ,ucont,ket)
      q_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
c
c set status variables
c
      status_nw = cando_nw(brain,ish,0).and.cando_nw(ketin,jsh,0)
      status_sp = cando_sp(brain,ish,0).and.cando_sp(ketin,jsh,0)
      status_sim = cando_sim(brain,ish,0).and.cando_sim(ketin,jsh,0)
cedo#if defined(USE_TEXAS_BROKE)
      status_txs = cando_txs(brain,ish,0).and.cando_txs(ketin,jsh,0)
      status_gen = (max(p_gen,q_gen)) .gt. 1 ! if general contraction texas is only option
*        greater than (f|f) do texas
      if (.not.status_gen)
     &    status_txs = status_txs .and. ((Lp+Lq).gt.6)
cedo#endif

c
      if (status_sp) then
        call genr70(
     &         brain,ish,coords(1,p_cent,p_geom),
     &                 0,coords(1,p_cent,p_geom),
     &         ketin,jsh,coords(1,q_cent,q_geom),
     &                 0,coords(1,q_cent,q_geom),
     &         eri)
cedo#if defined(USE_TEXAS_BROKE)
      else if (status_txs) then
        num_quart = 1
        dummy_lab(1) = 0
        dummy_lab(2) = 0
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        dum_log = OFALSE
        q4 = 1.0d00
        nint = 0
        txs_i = ish
        txs_j = jsh
        txs_d1 = 0
        txs_d2 = 0
        call texas_hf2_m(
     &      brain,txs_i,txs_d1,
     &      ketin,txs_j,txs_d2,num_quart,
     &      q4,OFALSE,roff,roff,roff,roff,OFALSE,
     &      eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
     &      nint, OFALSE, dum_log, scr, lscr, 0.0d0,'scfd_int')
        used_txs=.true.

cedo#endif
c
      else if(status_nw) then
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
        p_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        p_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        p_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(jsh,ket))
        Lq      = infbs_cont(CONT_TYPE ,ucont,ket)
        q_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        q_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        q_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        q_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c
        call hf2(
     &         coords(1,p_cent,p_geom), dbl_mb(mb_exndcf(p_iexp,bra)),
     &         dbl_mb(mb_exndcf(p_icfp,bra)), p_prim, p_gen, Lp,
     &         coords(1,p_cent,p_geom), DCexp,
     &         DCcoeff           ,      1,     1, 0,
     &         coords(1,q_cent,q_geom), dbl_mb(mb_exndcf(q_iexp,ket)),
     &         dbl_mb(mb_exndcf(q_icfp,ket)), q_prim, q_gen, Lq,
     &         coords(1,q_cent,q_geom), DCexp,
     &         DCcoeff           ,      1,     1, 0,
c......................... canAB    canCD    canPQ
     &         eri, leri, OFALSE, OFALSE, OFALSE,
c............. dryrun
     &         OFALSE, scr, lscr)
c
      elseif(status_sim)  then
        call nwcsim_hf2_2c(
     &        bra,ish,
     &        ket,jsh,
     &        nint, eri, leri, scr, lscr)

      else
        write(luout,*)'int_2e2c: could not do nw integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        call errquit('int_2e2c: fatal error ',0,0)
      endif
      if (any_spherical.and.(.not.used_txs))then
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
          if (Lp.eq.-1) p_gen = 1
          if (Lq.eq.-1) q_gen = 1 
          call spcart_2ctran(eri,scr,lscr,
     &        int_nbf_x(Lp),int_nbf_s(Lp),Lp,p_gen,bas_spherical(bra),
     &        int_nbf_x(Lq),int_nbf_s(Lq),Lq,q_gen,bas_spherical(ket),
     &        OFALSE)
       endif
      end
C> @}
