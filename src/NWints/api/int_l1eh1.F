c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute the 1 electron hamiltonian with zero integrals
C> removed
C>
C> This routine computes the 1 electron hamiltonian, (\f$H1\f$)
C> with labels and zero valued integrals removed. 
C> 
C> \f{eqnarray*}{
C> H1 & = & T + V      \\\\
C> T  & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\\\
C>    & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}
C> \nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> V  & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>    & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
C> {-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} 
C> \f}
C> 
C> If an ECP is defined then the ECP integral contributions are summed 
C> directly into the $H1$ integrals.  
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l1eh1}
c:tex-This routine computes the 1 electron hamiltonian, ($H1$)
c:tex-with labels and ``zero'' integrals removed. 
c:tex-
c:tex-\begin{eqnarray*}
c:tex-H1 & = & T + V      \\
c:tex-T  & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\
c:tex-   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}
c:tex-\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-V  & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
c:tex-{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} 
c:tex-\end{eqnarray*}
c:tex-
c:tex-If an ECP is defined then the ECP integral contributions are summed 
c:tex-directly into the $H1$ integrals.  
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l1eh1(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lh1, H1, lscr, scr, numgen) 
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "util.fh"
c
c::external subroutines used
c errquit
c::function
      logical int_chk_sh
      external int_chk_sh
      integer int_nint_cart
      external int_nint_cart
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision zerotol   !< [Input] zero threshold 
      double precision scr(lscr) !< [Scratch] scratch array
      integer lh1                !< [Input] length of 1e buffers.
      integer numgen             !< [Output] number of H1 integrals
      integer ilab(lh1)          !< [Output] i bas fun labels array
      integer jlab(lh1)          !< [Output] j bas fun labels array
      double precision H1(lh1)   !< [Output] 1e H1 integrals
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
     &       call errquit('int_l1eh1: invalid contraction/shell',0, 
     &       BASIS_ERR)
c      
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
      numint = int_nint_cart(i_basis,ish,j_basis,jsh, 0,0, 0,0)
      newlscr = lscr - numint
      call int_1eh1(i_basis, ish, j_basis, jsh,
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
              H1(numgen) = scr(icount)
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
