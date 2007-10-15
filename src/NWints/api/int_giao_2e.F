c $Id: int_giao_2e.F,v 1.2 2007-10-15 23:36:52 bert Exp $
*
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_2e4c}
c:tex-this routine computes the 4 center (traditional) 2 electron integrals: 
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2}g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_giao_2e(brain, ish, jsh, ketin, ksh, lsh,
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
      common/testdata/timing(20),irepeat
      double precision timing
      integer irepeat
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical cando_sp
      logical cando_txs
      external cando_nw
      external cando_sp
      external cando_txs
c:: passed
c:tex-\begin{verbatim}
      integer brain ! [input] bra basis set handle
      integer ish   ! [input] shell/contraction index
      integer jsh   ! [input] shell/contraction index
      integer ketin ! [input] ket basis set handle
      integer ksh   ! [input] shell/contraction index
      integer lsh   ! [input] shell/contraction index
      integer lscr  ! [input] length of scratch array
      double precision scr(lscr) ! [scratch] array
      integer leri  ! [input|output] length of integral array
      double precision eri(*) ! [output] 2e4c integrals
c:tex-\end{verbatim}
c:: local
      integer bra, ket 
      integer ab_geom, cd_geom, ucont, ityp
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_cent, a_icfps
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_cent, b_icfps
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent, c_icfps
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent, d_icfps
c
      double precision roff(3)
      double precision q4
      integer nint
      logical dum_log, do_bra, do_ket
      logical status_sp, status_nw, status_txs, status_gen
      integer texas_ang_limit
      logical ieqj, keql
      integer sbas, abas, bras, kets
c
      logical any_spherical
      integer a_nbf,   b_nbf,   c_nbf,   d_nbf
      integer a_nbf_s, b_nbf_s, c_nbf_s, d_nbf_s
      integer ab_gen, ab_cmp, cd_gen, cd_cmp,i,j
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c     timing(10)=timing(10)-util_wallsec()
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      bras = bra
      kets = ket
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      a_cent  = (sf_ibs_cn2ce(ish,bra))
      b_cent  = (sf_ibs_cn2ce(jsh,bra))
      c_cent  = (sf_ibs_cn2ce(ksh,ket))
      d_cent  = (sf_ibs_cn2ce(lsh,ket))
c
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
c
      status_sp = cando_sp(brain,ish,jsh).and.cando_sp(ketin,ksh,lsh)
      if (.not.status_sp) then
c
        ieqj = ish.eq.jsh
        keql = ksh.eq.lsh
c
        ucont   = sf_ibs_cn2ucn(ish,bra)
        La      = infbs_cont(CONT_TYPE ,ucont,bra)
        a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
        a_icfps = infbs_cont(CONT_ICFP ,ucont,bras)
c
        ucont   = sf_ibs_cn2ucn(jsh,bra)
        Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
        b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
        b_icfps = infbs_cont(CONT_ICFP ,ucont,bras)
c      
        ucont   = sf_ibs_cn2ucn(ksh,ket)
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        c_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        ucont   = sf_ibs_cn2ucn(lsh,ket)
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        d_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        a_nbf = int_nbf_x(La)
        b_nbf = int_nbf_x(Lb)
        c_nbf = int_nbf_x(Lc)
        d_nbf = int_nbf_x(Ld)
        leri=a_nbf*b_nbf*c_nbf*d_nbf*a_gen*b_gen*c_gen*d_gen
        call hnd_giahnd(
     &      coords(1,a_cent,ab_geom), dbl_mb(mb_exndcf(a_iexp,bra)),
     &      dbl_mb(mb_exndcf(a_icfp,bra)), a_prim, a_gen, La,
     &      coords(1,b_cent,ab_geom), dbl_mb(mb_exndcf(b_iexp,bra)),
     &      dbl_mb(mb_exndcf(b_icfp,bra)), b_prim, b_gen, Lb,
     &      coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &      dbl_mb(mb_exndcf(c_icfp,ket)), c_prim, c_gen, Lc,
     &      coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &      dbl_mb(mb_exndcf(d_icfp,ket)), d_prim,d_gen,Ld,
     &      ieqj, keql, eri, leri, scr, lscr)
c
c eri has cartesian block of integrals (llo:lhi,klo:khi,jlo:jhi,ilo:ihi,6)
c
        if (any_spherical) then
          a_nbf_s = int_nbf_s(La)
          b_nbf_s = int_nbf_s(Lb)
          c_nbf_s = int_nbf_s(Lc)
          d_nbf_s = int_nbf_s(Ld)
          cd_gen = c_gen*d_gen
          ab_gen = a_gen*b_gen
          ab_cmp = a_nbf*b_nbf
          cd_cmp = c_nbf*d_nbf
          do_bra=bas_spherical(bra)
          do_ket=bas_spherical(ket)
          call giao_to_sph(eri,leri,scr,La,Lb,Lc,Ld,a_nbf,b_nbf,c_nbf,
     &                     d_nbf,a_nbf_s,b_nbf_s,c_nbf_s,d_nbf_s,
     &                     a_gen,b_gen,c_gen,d_gen,ab_gen,ab_cmp,
     &                     cd_gen,cd_cmp,do_bra,do_ket)
        endif
      else
        write(luout,*)'int_giao_2e: cannot do sp integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        write(luout,*)' ksh   :',ksh
        write(luout,*)' lsh   :',lsh
        call errquit('int_giao_2e: fatal error ',0, INT_ERR)
      endif
      end
c
      subroutine giao_to_sph(eri,leri,scr,La,Lb,Lc,Ld,a_nbf,b_nbf,c_nbf,
     &                     d_nbf,a_nbf_s,b_nbf_s,c_nbf_s,d_nbf_s,
     &                     a_gen,b_gen,c_gen,d_gen,ab_gen,ab_cmp,
     &                     cd_gen,cd_cmp,do_bra,do_ket)
      implicit none
c
      integer leri
      double precision eri(leri,6), scr(*)
      integer La,Lb,Lc,Ld,a_nbf,b_nbf,c_nbf,d_nbf
      integer a_nbf_s,b_nbf_s,c_nbf_s,d_nbf_s,a_gen,b_gen,c_gen,d_gen
      integer ab_gen,ab_cmp,cd_gen,cd_cmp
      logical do_bra,do_ket
c
      integer ityp
c
      ab_cmp = a_nbf_s*b_nbf_s
      do ityp = 1, 6
         call spcart_bra2etran(eri(1,ityp),scr,
     &        b_nbf,a_nbf,b_nbf_s,a_nbf_s,
     &        Lb, La, b_gen, a_gen,
     &        cd_gen*cd_cmp,.false.)
         call spcart_ket2etran(eri(1,ityp),scr,
     &        d_nbf,c_nbf,d_nbf_s,c_nbf_s,
     &        Ld, Lc, d_gen, c_gen,
     &        ab_gen*ab_cmp,.false.)
      enddo
      a_nbf = a_nbf_s
      b_nbf = b_nbf_s
      c_nbf = c_nbf_s
      d_nbf = d_nbf_s
      end
c
