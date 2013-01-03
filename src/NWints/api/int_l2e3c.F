c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute non-zero 3 center 2-electron integrals
C>
C> Computes 3 center 2-electron integrals of the following kind:
C> \f{eqnarray*}{
C> ({\mu}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}
C> g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C> All zero valued integrals are removed and the remaining integrals
C> identified with appropriate labels.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l2e3c}
c:tex-this routine computes the 3 center 2 electron integrals
c:tex-with labels and ``zero'' integrals removed: 
c:tex-\begin{eqnarray*}
c:tex-({\mu}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{1}{r_{12}}
c:tex-g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l2e3c(brain, ish, ketin, jsh, ksh,
     &       zerotol, canket, leri, eri, nint, ilab, jlab, klab,
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
#include "stdio.fh"
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
      integer ksh   !< [Input] shell/contraction index
      integer lscr  !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] array
      double precision zerotol   !< [Input] zero threshold
      integer leri  !< [Input] length of integral array
      integer  nint !< [Output] number of integrals computed
      integer ilab(leri) !< [Output] i bas fun labels array
      integer jlab(leri) !< [Output] j bas fun labels array
      integer klab(leri) !< [Output] k bas fun labels array
      double precision eri(leri) !< [Output] 2e3c integrals
      logical canket !< [Input] canonicalize ket bas. fun. label pairs
c:tex-\end{verbatim}
c::local
      logical shells_ok
      integer i,j,k,bra,ket,icount
      integer numint, newlscr
c
#include "bas_ibs_sfn.fh"
c
c check shells 
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      if (.not.shells_ok)
     &       call errquit('int_l2e3c: invalid contraction/shell',0,
     &       BASIS_ERR)
c
c   check canonicalizations of input shells for canket
c      
      shells_ok = jsh.ge.ksh
      if (.not.shells_ok) then
        write(luout,*)'int_2e3c: shells not canonical on input '
        write(luout,*)'bra basis set handle:',brain
        write(luout,*)'ket basis set handle:',ketin
        write(luout,*)'                 ish:',ish
        write(luout,*)'                 jsh:',jsh
        write(luout,*)'                 ksh:',ksh
        call errquit('int_l2e3c: shells not in canonical order',0, 
     &     BASIS_ERR)
      endif
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
c
c compute eri (instead of copy) in front of scr()
c
      numint = int_nint_cart(brain,ish,ketin,jsh,ketin,ksh,0,0)
      newlscr = lscr - numint
      call int_2e3c(brain,ish,ketin,jsh,ksh,
     &    newlscr,scr(numint+1),numint,scr)
c
      nint = 0
      icount = 0
      do 00100 i = (sf_ibs_cn2bfr(1,ish,bra)),
     &      (sf_ibs_cn2bfr(2,ish,bra))
        do 00200 j = (sf_ibs_cn2bfr(1,jsh,ket)),
     &        (sf_ibs_cn2bfr(2,jsh,ket))
          do 00300 k = (sf_ibs_cn2bfr(1,ksh,ket)),
     &          (sf_ibs_cn2bfr(2,ksh,ket))
            icount = icount + 1
            if (abs(scr(icount)).ge.zerotol) then
              if (canket)then
                if(j.ge.k) then
                  nint = nint + 1
                  eri(nint) = scr(icount)
                  ilab(nint) = i
                  jlab(nint) = j
                  klab(nint) = k
                endif
              else
                nint = nint + 1
                eri(nint) = scr(icount)
                ilab(nint) = i
                jlab(nint) = j
                klab(nint) = k
              endif
            endif
00300     continue
00200   continue
00100 continue
c
      end
C> @}
