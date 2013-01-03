c $Id$
*
C> \ingroup nwint
C> @{
C> \brief Compute any 1-electron integrals removing any zero integrals
C>
C> This routine computes the 1 electron integrals \f$S\f$, \f$T\f$, and \f$V\f$:
C> \f{eqnarray*}{
C> S & = & ({\mu}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> T & = & ({\mu}|-\frac{1}{2}\nabla^{2}|{\nu}) \\\\
C>   & = & -\frac{1}{2}\int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\nabla^{2}(r_{1})\\
C>         g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
C> {-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} 
C> \f}
C> All zero valued integrals are removed and the remaining integrals
C> identified with appropriate labels.
C>
C> If an ECP is defined then the ECP integral contributions are summed
C> directly into the V integrals.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l1eall}
c:tex-This routine computes the 1 electron integrals S, T, and V 
c:tex-with labels and it removes ``zero'' integrals:
c:tex-\begin{eqnarray*} 
c:tex-S & = & ({\mu}|{\nu}) \\
c:tex-  & = & \int_{{-}\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-T & = & ({\mu}|\frac{-1}{2}{\nabla^{2}}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}{\nabla^{2}}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-Z_{\alpha}}
c:tex-{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*} 
c:tex-
c:tex-If an ECP is defined then the ECP integral contributions are summed 
c:tex-directly into the V integrals.  
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l1eall(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lstv, S, T, V, lscr, scr, numstv) 
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
#include "util.fh"
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
      double precision zerotol !< [Input] zero threshold for integrals
      double precision scr(lscr) !< [Scratch] scratch array
      integer lstv               !< [Input] length of one electron buffers
      integer ilab(lstv)         !< [Output] i bas. fun. index array
      integer jlab(lstv)         !< [Output] j bas. fun. index array
      double precision T(lstv)   !< [Output] kinetic integral buffer
      double precision V(lstv)   !< [Output] potential integral buffer
      double precision S(lstv)   !< [Output] overlap integral buffer
      integer numstv             !< [Output] number of returned integrals
c:tex-\end{verbatim}
c::local
      integer ibas, jbas, iscrS, iscrT, iscrV, icount, i, j
      logical ijbas, shells_ok
      logical nonzero
      integer numint, newlscr
c
#include "bas_ibs_sfn.fh"
c
c check shells
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_l1eall: invalid contraction/shell',0,
     &       BASIS_ERR)
c
c check canoncialization of shells
c
      ijbas = i_basis .eq. j_basis
      if (.not.(ijbas.and.ish.ge.jsh)) then
        write(luout,*)' shells not in canonical order '
        write(luout,*)' ish = ',ish
        write(luout,*)' jsh = ',jsh
        call errquit('int_l1eall: canonical error',0, INT_ERR)
      endif
c
      if ((3*lstv).gt.lscr) then
        write(luout,*)'int_l1eall: need more scratch space '
        write(luout,*)'int_l1eall: have :',lscr
        write(luout,*)'int_l1eall: need :',(3*lstv)
        call errquit('int_l1eall: scratch space error ',0, INT_ERR)
      endif
c
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
c
      numint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0, 0,0)
      iscrS = 1
      iscrT = iscrS + numint
      iscrV = iscrT + numint
      newlscr = lscr - 3*numint
      call int_1eall(i_basis,ish,j_basis,jsh,
     &    newlscr,scr(3*numint+1),numint,
     &    scr(iscrS),scr(iscrT),scr(iscrV))
c
      numstv = 0
      icount = 0
      do 00100 i = (sf_ibs_cn2bfr(1,ish,ibas)),
     &      (sf_ibs_cn2bfr(2,ish,ibas))
        do 00200 j = (sf_ibs_cn2bfr(1,jsh,jbas)),
     &        (sf_ibs_cn2bfr(2,jsh,jbas))
          nonzero = abs(scr((icount+iscrS))).ge.zerotol
          nonzero = nonzero .and.(abs(scr((icount+iscrT))).ge.zerotol)
          nonzero = nonzero .and.(abs(scr((icount+iscrV))).ge.zerotol)
          if (nonzero) then
            if ((.not.ijbas).or.i.ge.j) then
              numstv = numstv + 1
              S(numstv) = scr((icount+iscrS))
              T(numstv) = scr((icount+iscrT))
              V(numstv) = scr((icount+iscrV))
              ilab(numstv) = i
              jlab(numstv) = j
            endif
          endif
c
          icount = icount + 1
00200   continue
00100 continue
c
      end
C> @}
