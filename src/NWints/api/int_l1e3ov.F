c $Id$
*
C> \ingroup nwint
C> @{
C> \brief Calculates 3-center overlap integrals with labels and without
C> zero valued integrals
C>
C> This routine computes the 3 center overlap integral:
C> \f{eqnarray*}{
C> ({\mu}{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})g_{\lambda}(X_{\lambda},r_{1})dr_{1}
C> \f}
C> Integrals that are zero are removed and the remaining integrals are
C> identified with appropriate labels.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l1e3ov}
c:tex-This routine computes the 3 center overlap integral 
c:tex-with labels and it removes ``zero'' integrals:
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1}g_{\lambda}(X_{\lambda},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l1e3ov(i_basis, ish, j_basis, jsh, k_basis, ksh,
     &       zerotol, lov3, OV3, ilab, jlab, klab, numov3, lscr, scr)
c:tex-\end{verbatim}
c
c routine to calculate labels for 3 center overlap integrals 
c
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
c..  errquit
c::function
      logical int_chk_sh
      external int_chk_sh
      integer int_nint_cart
      external int_nint_cart
c::passed
c:tex-\begin{verbatim}
      integer i_basis            !< [Input] basis set handle for ish
      integer ish                !< [Input] i shell/contraction
      integer j_basis            !< [Input] basis set handle for jsh
      integer jsh                !< [Input] j shell/contraction
      integer k_basis            !< [Input] basis set handle for ksh
      integer ksh                !< [Input] k shell/contraction
      double precision zerotol   !< [Input] "zero" tolerance threshold
      integer lov3               !< [Input] length of 3c overlap array
      double precision OV3(lov3) !< [Output] 3c overlap integrals
      integer ilab(lov3)         !< [Output] i labels for 3c overlap ints
      integer jlab(lov3)         !< [Output] j labels for 3c overlap ints
      integer klab(lov3)         !< [Output] k labels for 3c overlap ints
      integer numov3             !< [Output] number of integrals and
                                 !< labels generated and returned
      integer lscr               !< [Input] length of scratch array
      double precision scr(lscr) !< [Input] scratch array
c:tex-\end{verbatim}
c::local
      integer ibas, jbas, kbas  ! real basis set array index
      integer i, j, k, icount   ! and do loop variables counters
      logical shells_ok         ! all shells valid
      logical nonzero           ! for nonzero check
      integer numint, newlscr
c
      logical ijbas, jkbas, ikbas, allbas  ! check of basis set handles for canonicalization
c
#include "bas_ibs_sfn.fh"
c
c
c check shells
c
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      shells_ok = shells_ok .and. int_chk_sh(k_basis,ksh)
      if (.not.shells_ok)
     &       call errquit('int_l1e3ov: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      kbas = k_basis + BASIS_HANDLE_OFFSET
      ijbas  = ibas.eq.jbas
      jkbas  = jbas.eq.kbas
      ikbas  = ibas.eq.kbas
      allbas = ijbas .and. jkbas
c
      numint = int_nint_cart(i_basis,ish,j_basis,jsh,k_basis,ksh,0,0)
      newlscr = lscr-numint
      call int_1e3ov(i_basis, ish, j_basis, jsh, k_basis, ksh,
     &    newlscr,scr(numint+1),numint,scr)
c
      numov3 = 0
      icount = 0
c
      do 00100 i = (sf_ibs_cn2bfr(1,ish,ibas)),
     &      (sf_ibs_cn2bfr(2,ish,ibas))
        do 00200 j = (sf_ibs_cn2bfr(1,jsh,jbas)),
     &        (sf_ibs_cn2bfr(2,jsh,jbas))
          do 00300 k = (sf_ibs_cn2bfr(1,ksh,kbas)),
     &          (sf_ibs_cn2bfr(2,ksh,kbas))
            icount = icount + 1
            nonzero = abs(scr(icount)).ge.zerotol
            if (nonzero) then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .!  keep the integral
              if     (allbas) then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .!  i>=j>=k
                if(i.ge.j.and.j.ge.k) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              elseif (ijbas)  then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .! i>=j ; free k
                if (i.ge.j) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              elseif (jkbas)  then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .! j>=k ; free i
                if (j.ge.k) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              elseif (ikbas)  then
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . .! i>=k ; free j  (least likely!!)
                if (i.ge.k) then
                  numov3 = numov3 + 1
                  OV3(numov3)  = scr(icount)
                  ilab(numov3) = i
                  jlab(numov3) = j
                  klab(numov3) = k
                endif
              else
                numov3 = numov3 + 1
                OV3(numov3)  = scr(icount)
                ilab(numov3) = i
                jlab(numov3) = j
                klab(numov3) = k
              endif
            endif
00300     continue
00200   continue
00100 continue
c      
      end
C> @}
