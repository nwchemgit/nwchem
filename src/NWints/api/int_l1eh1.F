c $Id: int_l1eh1.F,v 1.5 1996-10-07 16:50:33 d3e129 Exp $
*
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsubsection{int\_l1eh1}
c:tex-This routine computes the 1 electron hamiltonian, ($H1$)
c:tex-with labels and ``zero'' integrals removed. 
c:tex-\[
c:tex-H1 = T + V      
c:tex-\]
c:tex-\[
c:tex-T = ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\]
c:tex-\[
c:tex-V = ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} 
c:tex-\]
c:tex-If an ecp is defined then the integral contributions are summed 
c:tex-directly into the $H1$ integrals.  
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l1eh1(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lh1, H1, lscr, scr, numgen) 
c:tex-\end{verbatim}
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c
c::external subroutines used
c errquit
c::function
      logical int_chk_sh
      external int_chk_sh
c::passed
c:tex-\begin{verbatim}
      integer i_basis ! [input] basis set handle for ish
      integer ish     ! [input] i shell/contraction
      integer j_basis ! [input] basis set handle for jsh
      integer jsh     ! [input] j shell/contraction
      integer lscr    ! [input] length of scratch array
      double precision zerotol   ! [input] zero threshold 
      double precision scr(lscr) ! [scratch] scratch array
      integer lh1                ! [input] length of one electron buffer.
      integer numgen             ! [output] number of H1 integrals
      integer ilab(lh1)          ! [output] i bas fun labels array
      integer jlab(lh1)          ! [output] j bas fun labels array
      double precision H1(lh1)   ! [output] one electron hamiltonian integrals
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
     &       call errquit('int_l1eh1: invalid contraction/shell',0)
c
c  check if spherical
c
      call int_nospherical_check(i_basis,'int_l1eh1')
      call int_nospherical_check(j_basis,'int_l1eh1')
c      
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
      numint =         (sf_ibs_cn2bfr(2,ish,ibas))-
     &      (sf_ibs_cn2bfr(1,ish,ibas))+1
      numint = numint*((sf_ibs_cn2bfr(2,jsh,jbas))-
     &      (sf_ibs_cn2bfr(1,jsh,jbas))+1)
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
