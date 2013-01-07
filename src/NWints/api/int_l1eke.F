c $Id$
*
C>
C> \ingroup nwint
C> @{
C>
C> \brief Compute the non-zero 1 electron kinetic energy integrals
C>
C> This routine computes the 1 electron kinetic energy integrals, 
C> (\f$T\f$),
C> with labels and zero valued integrals removed:
C> 
C> \f{eqnarray*}{
C> T & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l1eke}
c:tex-This routine computes the 1 electron kinetic energy integrals, ($T$).
c:tex-with labels and ``zero'' integrals removed:
c:tex-
c:tex-\begin{eqnarray*}
c:tex-T & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l1eke(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lke, Ke, lscr, scr, numgen) 
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
#include "util.fh"
c
c::external subroutines used
c errquit
c::function
      logical int_chk_sh
      integer int_nint_cart
      external int_chk_sh, int_nint_cart
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      double precision zerotol   !< [Input] zero threshold
      integer lke                !< [Input] length of 1e buffers
      integer numgen             !< [Output] number of Ke integrals
      integer ilab(lke)          !< [Output] i bas fun labels array
      integer jlab(lke)          !< [Output] j bas fun labels array
      double precision Ke(lke)   !< [Output] kinetic energy integrals
c:tex-\end{verbatim}
c::local
      integer ibas, jbas, icount, i, j
      logical ijbas
      logical shells_ok
      logical nonzero
      integer numint, newlscr
c
#include "bas_ibs_sfn.fh"
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_l1eke: invalid contraction/shell',0,
     &        BASIS_ERR)
c
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
      numint =   int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      newlscr = lscr - numint
      call int_1eke(i_basis, ish, j_basis, jsh,
     &    newlscr,scr(numint+1),numint,scr)
c
      numgen = 0
      icount = 0
      do 00100 i = (sf_ibs_cn2bfr(1,ish,ibas)),
     &      (sf_ibs_cn2bfr(2,ish,ibas))
        do 00200 j = (sf_ibs_cn2bfr(1,jsh,jbas)),
     &        (sf_ibs_cn2bfr(2,jsh,jbas))
          icount = icount + 1
          nonzero = abs(scr(icount)).ge.zerotol
          if (nonzero) then
            if ((.not.ijbas).or.i.ge.j) then
              numgen = numgen + 1
              Ke(numgen) = scr(icount)
              ilab(numgen) = i
              jlab(numgen) = j
            endif
          endif
c
00200   continue
00100 continue
c
      end
C> @}
