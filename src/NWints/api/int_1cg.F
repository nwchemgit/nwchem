c $Id$
*
**********************************************************************
C> \ingroup nwint
C> @{
C>
C> \brief Integrate Gaussians over all space
C>
C> This routine computes the integral of the one center gaussian.
C> \f{eqnarray*}{
C> ({\mu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})dr_{1}
C> \f}
C>
c:tex-% This is part of the API Standard Integral routines
c:tex-\subsection{int\_1cg}
c:tex-This routine computes the integral of the one center gaussian.
c:tex-\begin{eqnarray*}
c:tex-({\mu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1cg(i_basis,ish,lscr,scr,l1g,G1e)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "util.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle
      integer ish     !< [Input] lexical shell/contraction index
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr)  !< [Scratch] scratch space 
      integer l1g     !< [Input] length of integral buffer array
      double precision G1e(l1g)   !< [Output] one center gaussian integrals
c:tex-\end{verbatim}
c::local
      logical any_spherical
      logical shells_ok
      integer i_geom, ibas, ucont
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer i_nbf_x, i_nbf_s
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1cg'))
     &       call errquit('int_1cg: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      if (.not.shells_ok)
     &       call errquit('int_1cg: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
c
      if (cando_nw_1e(i_basis,ish,0)) then
        call int_hf1sp(
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,i_cent,i_geom), DCexp, DCcoeff,
     &        1, 1, 0, i_cent,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
c............................. doS    doT     doV     canonical
     &        G1e,scr,scr,l1g,.true.,.false.,.false.,.false.,
c.............. dryrun
     &        .false.,scr,lscr,'int_1cg')
      elseif (cando_nw(i_basis,ish,0)) then
        call hf1(
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li,
     &        coords(1,i_cent,i_geom),DCexp,
     &        DCcoeff, 1, 1, 0,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
c............................. doS    doT     doV     canonical
     &         G1e,scr,scr,l1g,.true.,.false.,.false.,.false.,
c.............. dryrun
     &         .false.,scr,lscr)
      else
        call errquit('int_1cg: could not do sp or nw integrals',0,
     &       INT_ERR)
      endif
c
*     G1e now has the cartesian integral block  (ilo:ihi)
*      
      any_spherical = bas_spherical(ibas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
c ... since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
c
*.. transform on i component
      i_nbf_x = int_nbf_x(Li)
      i_nbf_s = int_nbf_s(Li)
      call spcart_tran1e(G1e,scr,
     &    1,i_nbf_x,0,1,
     &    1,i_nbf_s,Li,i_gen,
     &    .false.)
c
      end
C> @}
