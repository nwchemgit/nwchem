c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Calculate a block of non-zero 4 center 2-electron integrals
C>
C> The routine takes 4 shells of basis functions. The first pair of
C> shells generates the "bra" charge distribution and the second pair
C> the "ket" charge distribution. This routine evaluates the integral
C> of the Coulomb operator for the repulsion of these charge
C> distributions. In summary the following expression is evaluated:
C> \f{eqnarray*}{
C>    ({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty}
C>    g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}
C>    g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C> All zero valued integrals are removed and the remaining integrals
C> identified with appropriate labels.
C>
C> At present this routine provides the following algorithms:
C>
C> - SP-rotated axis 
C>
C> - McMurchie-Davidson
C>
C> - Obara-Saika (Texas)
C>
C> Not available are:
C>
C> - Rys quadrature (Hondo)
C>
C> The choice of algorithm is controlled by capability checks and user
C> inputs (see the various cando functions for details).
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_l2e4c}
c:tex-this routine computes the 4 center (traditional) 2 electron integrals
c:tex-with labels and ``zero'' integrals removed: 
c:tex-
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_l2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       zerotol, canonicalize, leri, eri, nint, ilab, jlab, klab,
     &       llab, lscr, scr)
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
      logical  cando_sp, cando_nw, cando_txs
      logical  int_chk_sh
      integer int_nint_cart
      external int_chk_sh
      external cando_sp, cando_nw, cando_txs
      external int_nint_cart
c::passed
c:tex-\begin{verbatim}
      integer brain !< [Input] bra basis set handle
      integer ish   !< [Input] shell/contraction index
      integer jsh   !< [Input] shell/contraction index
      integer ketin !< [Input] ket basis set handle
      integer ksh   !< [Input] shell/contraction index
      integer lsh   !< [Input] shell/contraction index
      double precision zerotol !< [Input] zero threshold
      integer lscr  !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] array
      integer leri  !< [Input] length of integral array
      integer nint  !< [Output] number of integrals computed
      integer ilab(leri) !< [Output] i bas fun label array 
      integer jlab(leri) !< [Output] j bas fun label array 
      integer klab(leri) !< [Output] k bas fun label array 
      integer llab(leri) !< [Output] l bas fun label array 
      double precision eri(leri) !< [Output] 2e4c integrals
      logical canonicalize !< [Input] canonicalize labels
c:tex-\end{verbatim}
c::local
      logical status_sp, status_txs
      logical braket, shells_ok
#if defined(INT_LIMIT_CHECK)
      logical allshl, ignore_ij_kl
#endif
*       integer ijsh, klsh
      integer isym2,iscr1
      integer i,j,k,l,bra,ket,icount,ij,kl
      integer numint, newlscr
      logical dum_log_f, dum_log_t
      integer num_quart, intp 
      integer txs_i,txs_j,txs_k,txs_l
      double precision roff(3), val, q4

c
c:statement function
c
#include "bas_ibs_sfn.fh"
      isym2(i,j)=(max(i,j)*(max(i,j)-1))/2 + min(i,j)
c
c check initialization and shells 
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(brain,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,lsh)
      if (.not.shells_ok)
     &       call errquit('int_l2e4c: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      status_sp  = cando_sp(brain, ish, jsh)  .and.
     &    cando_sp(ketin, ksh, lsh)
      status_txs = cando_txs(brain, ish, jsh) .and.
     &    cando_txs(ketin, ksh, lsh)          .and.
     &    brain.eq.ketin                                 ! right now texas uses only one basis
c
      if (status_txs .and. (.not.status_sp)) then
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        num_quart = 1
        dum_log_f = .false.
        dum_log_t = .true.
        q4 = 1.0d00
        txs_i = ish
        txs_j = jsh
        txs_k = ksh
        txs_l = lsh
        call texas_hf2_m (
     &      brain, txs_i, txs_j,
     &      ketin, txs_k, txs_l, num_quart,
     &      q4, dum_log_f,
     &      roff, roff, roff, roff, dum_log_f,
     &      eri, leri,
     &      ilab, jlab, klab, llab, nint, dum_log_t,
     &      dum_log_f, scr, lscr, 0.0d0,'scfd_int')
#if defined(INT_LIMIT_CHECK)
        braket = brain.eq.ketin
        allshl = isym2(ish,jsh) .eq. isym2(ksh,lsh)
        ignore_ij_kl = (.not.braket) .or. (.not.allshl)
        icount = 0
        do 10100 intp = 1,nint
          i = ilab(intp)
          j = jlab(intp)
          k = klab(intp)
          l = llab(intp)
          ij = isym2(i,j)
          kl = isym2(k,l)
          val = eri(intp)
          if (abs(val).ge.zerotol) then
            if (i.ge.j.and.k.ge.l.and.
     &          (ignore_ij_kl.or.(ij.ge.kl)))then
              icount = icount + 1
              if (icount.lt.intp) then
                eri(icount)  = val
                if (canonicalize.and.braket) then
                  call int_canon(i,j,k,l,
     &                ilab(icount), jlab(icount),
     &                klab(icount), llab(icount))
                else
                  ilab(icount) = i
                  jlab(icount) = j
                  klab(icount) = k
                  llab(icount) = l
                endif
              else if (canonicalize.and.braket) then
                call int_canon(i,j,k,l,
     &              ilab(intp), jlab(intp),
     &              klab(intp), llab(intp))
              endif
            endif
          endif
10100   continue
        nint = icount
#else
        braket = brain.eq.ketin
        icount = 0
        do 10100 intp = 1,nint
          i = ilab(intp)
          j = jlab(intp)
          k = klab(intp)
          l = llab(intp)
          ij = isym2(i,j)
          kl = isym2(k,l)
          val = eri(intp)
          if (abs(val).ge.zerotol) then
            icount = icount + 1
            eri(icount)  = val
            if (canonicalize.and.braket) then
              call int_canon(i,j,k,l,
     &            ilab(icount), jlab(icount),
     &            klab(icount), llab(icount))
            else if (canonicalize) then
              if (i.lt.j) then
                ilab(icount) = j
                ilab(icount) = i
              else
                ilab(icount) = i
                jlab(icount) = j
              endif
              if (k.lt.l) then
                klab(icount) = l
                llab(icount) = k
              else
                klab(icount) = k
                llab(icount) = l
              endif
            else
              ilab(icount) = i
              jlab(icount) = j
              klab(icount) = k
              llab(icount) = l
            endif
          endif
10100   continue
        nint = icount
#endif
      else
        
c
c compute eri (instead of copy) in front of scr()
c
c
        bra = brain + BASIS_HANDLE_OFFSET
        ket = ketin + BASIS_HANDLE_OFFSET
        
        numint = int_nint_cart( brain, ish, brain, jsh,
     &       ketin, ksh, ketin, lsh)
        iscr1=numint
#ifdef USE_SIMINT
c     this avoids segv in simint SIMD
#define SIMINT_SIMD_LEN 8
        call util_align(iscr1,SIMINT_SIMD_LEN)
#endif
        newlscr = lscr - iscr1
        if(newlscr.lt.0) call errquit(' intl2e4c: insuff scratch',
     c       newlscr,0)
        call int_2e4c(brain,ish,jsh,ketin,ksh,lsh,
c     &      newlscr,scr(numint+1),numint,scr)
csimint    two args    needs to be  aligned with SIMINT_SIMD_LEN
     &      newlscr,scr(iscr1+1),numint,scr)
c
#if defined(INT_LIMIT_CHECK)
c
        braket = bra .eq. ket
        allshl = isym2(ish,jsh) .eq. isym2(ksh,lsh)
        ignore_ij_kl = (.not.braket) .or. (.not.allshl)
c
        nint = 0
        icount = 0
        do 00100 i = (sf_ibs_cn2bfr(1,ish,bra)),
     &        (sf_ibs_cn2bfr(2,ish,bra))
          do 00200 j = (sf_ibs_cn2bfr(1,jsh,bra)),
     &          (sf_ibs_cn2bfr(2,jsh,bra))
            ij = isym2(i,j)
            do 00300 k = (sf_ibs_cn2bfr(1,ksh,ket)),
              (sf_ibs_cn2bfr(2,ksh,ket))
              do 00400 l = (sf_ibs_cn2bfr(1,lsh,ket)),
     &              (sf_ibs_cn2bfr(2,lsh,ket))
                kl = isym2(k,l)
                icount = icount + 1
                if (abs(scr(icount)).ge.zerotol) then
                  if (i.ge.j.and.k.ge.l.and.
     &                (ignore_ij_kl.or.(ij.ge.kl)))then
                    nint = nint + 1
                    eri(nint) = scr(icount)
                    if (canonicalize.and.braket) then
                      call int_canon(i,j,k,l,
     &                    ilab(nint),jlab(nint),
     &                    klab(nint),llab(nint))
                    else
                      ilab(nint) = i
                      jlab(nint) = j
                      klab(nint) = k
                      llab(nint) = l
                    endif
                  endif
                endif
00400         continue
00300       continue
00200     continue
00100   continue
#else
        braket = bra .eq. ket
        nint = 0
        icount = 0
        do 00100 i = (sf_ibs_cn2bfr(1,ish,bra)),
     &        (sf_ibs_cn2bfr(2,ish,bra))
          do 00200 j = (sf_ibs_cn2bfr(1,jsh,bra)),
     &          (sf_ibs_cn2bfr(2,jsh,bra))
            do 00300 k = (sf_ibs_cn2bfr(1,ksh,ket)),
     &            (sf_ibs_cn2bfr(2,ksh,ket))
              do 00400 l = (sf_ibs_cn2bfr(1,lsh,ket)),
     &              (sf_ibs_cn2bfr(2,lsh,ket))
                icount = icount + 1
                if (abs(scr(icount)).ge.zerotol) then
                  nint = nint + 1
                  eri(nint) = scr(icount)
                  if (canonicalize.and.braket) then
                      call int_canon(i,j,k,l,
     &                    ilab(nint),jlab(nint),
     &                    klab(nint),llab(nint))
                  elseif (canonicalize) then
                    if (i.lt.j) then
                      ilab(nint) = j
                      ilab(nint) = i
                    else
                      ilab(nint) = i
                      jlab(nint) = j
                    endif
                    if (k.lt.l) then
                      klab(nint) = l
                      llab(nint) = k
                    else
                      klab(nint) = k
                      llab(nint) = l
                    endif
                  else
                    ilab(nint) = i
                    jlab(nint) = j
                    klab(nint) = k
                    llab(nint) = l
                  endif
                endif
00400         continue
00300       continue
00200     continue
00100   continue
#endif
c
      endif
      end
C> @}
