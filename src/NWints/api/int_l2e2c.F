c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute the non-zero 2 center 2-electron integrals.
C>
C> Computes 2 center 2-electron integrals of the form:
C> \f{eqnarray*}{
C> ({\mu}|{\nu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})dr_{1}dr_{2}
C> \f}
C> The zero valued integrals are removed and the remaining integrals
C> identified with appropriate labels.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l2e2c}
c:tex-this routine computes the 2 center 2 electron integrals
c:tex-with labels and ``zero'' integrals removed: 
c:tex-\begin{eqnarray*}
c:tex-({\mu}|{\nu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l2e2c(brain, ish, ketin, jsh, 
     &       zerotol, leri, eri, nint, ilab, jlab, 
     &       lscr, scr)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
#include "util.fh"
c
c::external subroutines used
c errquit
c::functions
      integer int_nint_cart
      external int_nint_cart
      logical int_chk_sh
      external int_chk_sh
c::passed
c:tex-\begin{verbatim}
      integer brain !< [Input] bra basis set handle
      integer ish   !< [Input] shell/contraction index
      integer ketin !< [Input] ket basis set handle
      integer jsh   !< [Input] shell/contraction index
      integer lscr  !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] array
      double precision zerotol   !< [Input] zero threshold
      integer leri  !< [Input] length of integral array
      integer nint  !< [Output] num of ints computed
      integer ilab(leri) !< [Output] i bas func label array
      integer jlab(leri) !< [Output] j bas func label array
      double precision eri(leri) !< [Output] 2e2c integrals
c:tex-\end{verbatim}
c::local
      logical shells_ok
      integer i,j,bra,ket,icount
      integer numint, newlscr
c
c:statement function
c
#include "bas_ibs_sfn.fh"
c
c check shells 
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      if (.not.shells_ok)
     &       call errquit('int_l2e2c: invalid contraction/shell',0,
     &         BASIS_ERR)
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
c compute eris into scr instead of copying them
c
c 
      numint = int_nint_cart(brain,ish,ketin,jsh,0,0,0,0)
      newlscr = lscr - numint
      call int_2e2c(brain,ish,ketin,jsh,
     &    newlscr,scr(numint+1),numint,scr)
c
      nint = 0
      icount = 0
      do 00100 i = (sf_ibs_cn2bfr(1,ish,bra)),
     &      (sf_ibs_cn2bfr(2,ish,bra))
        do 00200 j = (sf_ibs_cn2bfr(1,jsh,ket)),
     &        (sf_ibs_cn2bfr(2,jsh,ket))
          icount = icount + 1
          if (abs(scr(icount)).ge.zerotol) then
            if (bra.eq.ket) then
              if (i.ge.j) then
                nint = nint + 1
                eri(nint) = scr(icount)
                ilab(nint) = i
                jlab(nint) = j
              endif
            else
              nint = nint + 1
              eri(nint) = scr(icount)
              ilab(nint) = i
              jlab(nint) = j
            endif
          endif
00200   continue
00100 continue
c
      end
C> @}
