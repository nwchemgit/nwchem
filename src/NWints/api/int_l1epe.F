c $Id: int_l1epe.F,v 1.10 2003-10-17 22:53:42 carlfahl Exp $
*
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l1epe}
c:tex-This routine computes the 1 electron potential integrals, ($V$):
c:tex-If an ECP is defined then the ECP integral contributions are summed 
c:tex-directly into the V integrals.  Integrals are computed with labels 
c:tex-and ``zero'' integrals removed.  
c:tex-
c:tex-\begin{eqnarray*}
c:tex-V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
c:tex-{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} 
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l1epe(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lpe, Pe, lscr, scr, numgen) 
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
      external int_chk_sh
      integer int_nint_cart
      external int_nint_cart
c::passed
c:tex-\begin{verbatim}
      integer i_basis ! [input] basis set handle for ish
      integer ish     ! [input] i shell/contraction
      integer j_basis ! [input] basis set handle for jsh
      integer jsh     ! [input] j shell/contraction
      integer lscr    ! [input] length of scratch array
      double precision scr(lscr) ! [scratch] scratch array
      double precision zerotol   ! [input] zero integral threshold
      integer lpe                ! [input] length of potential buffer
      integer numgen             ! [output] number of integrals generated
      integer ilab(lpe)          ! [output] i bas fun labels array
      integer jlab(lpe)          ! [output] j bas fun labels array
      double precision Pe(lpe)   ! [output] potential integrals
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
     &       call errquit('int_l1epe: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
      numint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0, 0,0)
      newlscr = lscr - numint
* compute pe integrals
      call int_1epe(i_basis, ish, j_basis, jsh,
     &    newlscr, scr(numint+1), numint,scr)
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
              Pe(numgen) = scr(icount)
              ilab(numgen) = i
              jlab(numgen) = j
            endif
          endif
c
00200   continue
00100 continue
c
      end
