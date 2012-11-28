*     $Id$
C> \ingroup nwint
C> @{
C>
C> \brief Compute the integrals of a block of shell quartets and return
C> the non-zero integrals and corresponding labels
C>
C> Compute and return the non-zero integrals and corresponding labels
C> based on the input of a group of shell quartets (a block). This 
C> function cannot be without a call to `intb_init4c` with the same block
C> information. The data returned is compressed in the sense that all
C> zero valued integrals have been removed.
C>
C> The integrals evaluated are defined as
C> \f{eqnarray*}{
C> ({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C> The following code excerpt describes the proper use of the blocking
C> API routines.
C> 
C> * begin atom/shell loops
C> \code
C>        call collect_group_of_shells()
C>        okay = intb_init4c($\cdots$) ! with group of shells info     
C>        if (.not.okay) stop ' error setting up blocking interface '
C>  00001 continue
C>        more = intb_2e4c($\cdots$) ! with group of shells info
C>        call use_integrals_in_some_way()
C>        if (more) goto 00001
C> \endcode
C> 
C> * end atom/shell loops
C>
C> \return Returns .true. if not all integrals fitted in the buffer
C> provided and an additional call to `intb_2e4c` is required, and
C> returns .false. if all integrals were retrieved.
c:tex-
c:tex-\subsection{intb\_2e4c}
c:tex-This logical function returns the integrals and labels based on the
c:tex-input of a group of shell quartets (a block).  This function cannot be
c:tex-called without a call to intb\_init4c using the same block information.
c:tex-``zero'' integrals are removed.
c:tex-
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2})g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Return Values:}
c:tex-\begin{tabular}{|c|p{60mm}|} \hline
c:tex- .true.  & more integrals from this block another call to intb\_2e4c required \\ \hline
c:tex- .false. & all integrals from this block computed\\ \hline
c:tex-\end{tabular}
c:tex-
c:tex-The following code excerpt describes the proper use of the blocking
c:tex-API routines.\\
c:tex-{\it Pseudo Code:}
c:tex-\begin{verbatim}
c:tex-
c:tex-* begin atom/shell loops
c:tex-
c:tex-      call collect_group_of_shells()
c:tex-      okay = intb_init4c($\cdots$) ! with group of shells info     
c:tex-      if (.not.okay) stop ' error setting up blocking interface '
c:tex-00001 continue
c:tex-      more = intb_2e4c($\cdots$) ! with group of shells info
c:tex-      call use_integrals_in_some_way()
c:tex-      if (more) goto 00001
c:tex-
c:tex-* end atom/shell loops
c:tex-
c:tex-\end{verbatim}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function intb_2e4c(brain, icl, jcl, ketin, kcl, lcl,
     $     num_q, q4, use_q4, zerotol, canonicalize,
     $     ilab, jlab, klab, llab, eri,
     $     l_erilab, nints, lscr, scr)
c:tex-\end{verbatim}
      implicit none
#include "numb_qP.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
#include "util.fh"
#include "mafdecls.fh"
c::   functions
      logical intb_nw_2e4c
      external intb_nw_2e4c
c::   passed
c:tex-\begin{verbatim}
      integer brain              !< [Input]  basis handle for bra 
      integer ketin              !< [Input]  basis handle for ket 
      integer num_q              !< [Input]  number of quartets input
      integer icl(num_q)         !< [Input]  i-contraction quartet labels
      integer jcl(num_q)         !< [Input]  j-contraction quartet labels
      integer kcl(num_q)         !< [Input]  k-contraction quartet labels
      integer lcl(num_q)         !< [Input]  l-contraction quartet labels
      logical use_q4             !< [Input]  logical for use of q4
      double precision q4(num_q) !< [Input] symmetry prefactors for ints
*     
      integer l_erilab           !< [Input]  length of eri and label arrays
      integer ilab(l_erilab)     !< [Output] integral labels for ``i'' 
      integer jlab(l_erilab)     !< [Output] integral labels for ``j'' 
      integer klab(l_erilab)     !< [Output] integral labels for ``k'' 
      integer llab(l_erilab)     !< [Output] integral labels for ``l'' 
      integer nints              !< [Output] number of integrals returned
      double precision eri(l_erilab) !< [Output] integral values
      integer lscr               !< [Input] length of scratch array
      double precision scr(lscr) !< [Input] scratch array
      double precision zerotol   !< [Input] zero integral threshold
      logical canonicalize       !< [Input] Canonicalize integral labels?
c:tex-\end{verbatim}
c::   local
      double precision roff(3)
c     
c     Cannot do texas.
c     
      intb_2e4c = .false.
      if (num_q.eq.0) then
        nints = 0
        return
      endif
      if (brain.ne.ketin .or. (.not.txs_b2e4c)) then
         intb_2e4c = intb_nw_2e4c(
     $        brain, icl, jcl, ketin, kcl, lcl,
     $        num_q, q4, use_q4, zerotol, canonicalize,
     $        ilab, jlab, klab, llab, eri,
     $        l_erilab, nints, lscr, scr)
         return
      endif
c     
c     If faster to use SP than texas
c     
      if (nsp .gt. 0) then
         intb_2e4c = intb_nw_2e4c(
     $        brain, icl, jcl, ketin, kcl, lcl,
     $        nsp, q4, use_q4, zerotol, canonicalize,
     $        ilab, jlab, klab, llab, eri,
     $        l_erilab, nints, lscr, scr)
         if ((nints.eq.0) .or. (.not. intb_2e4c)) then
c     
c     Must skip past those already done with SP
c     
            num_q_start = nsp + 1
            nsp = 0
            if (num_q_start .gt. num_q) then
               intb_2e4c = .false.
               num_q_start = -99999999
            else
               intb_2e4c = .true.
            endif
            return
         else
            return
         endif
      endif
c     
      intb_2e4c = .true.
c     
      call dcopy(3,0.0d00,0,roff,1)
      call texas_hf2_m(
     $     brain, icl(num_q_start), jcl(num_q_start),
     $     ketin, kcl(num_q_start), lcl(num_q_start), 
     $     num_q-num_q_start+1,
     $     q4(num_q_start), use_q4,
     $     roff,roff,roff,roff,.false.,
     $     eri,l_erilab, ilab, jlab, klab, llab, nints, .true.,
     $     more_ints_txs, scr, lscr, zerotol,'scfd_int')
      if (nints.gt.l_erilab) then
         write(luout,*)' number of integrals from texas   :',nints
         write(luout,*)' size of integral and label buffer:',l_erilab
         write(luout,*)' you must increase eri buffer/label size '
         call errquit(' integral buffer overflow ',911, MEM_ERR)
      endif
      intb_2e4c = more_ints_txs
c     
      end
C>
C> \brief Compute the integrals of a block of shell quartets and return
C> the non-zero integrals and corresponding labels using the NWChem 
C> McMurchie-Davidson code.
C>
C> Compute and return the non-zero integrals and corresponding labels
C> based on the input of a group of shell quartets (a block). This 
C> routine explicitly uses the McMurchie-Davidson code, it should \b NEVER
C> be called directly by an application module! Please call `intb_2e4c`
C> instead.
c:tex-
c:tex-\subsection{intb\_nw\_2e4c}
c:tex-This logical function returns the integrals and labels based on the
c:tex-input of a group of shell quartets (a block).  This interfaces to
c:tex-the NWChem McMurchie-Davidson code.  This routine should {\it never}
c:tex-be called directly by an application module!!
c:tex-``zero'' integrals are removed.
c:tex-
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}g_{\nu}(X_{\nu},r_{2}g_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Return Values:}
c:tex-\begin{tabular}{|c|p{60mm}|} \hline
c:tex- .true.  & more integrals from this block another call to intb\_2e4c required \\ \hline
c:tex- .false. & all integrals from this block computed\\ \hline
c:tex-\end{tabular}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function intb_nw_2e4c(brain, icl, jcl, ketin, kcl, lcl,
     $     num_q, q4, use_q4, zerotol, canonicalize,
     $     ilab, jlab, klab, llab, eri,
     $     l_erilab, nints, lscr, scr)
c:tex-\end{verbatim}
      implicit none
#include "numb_qP.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::   functions
      integer int_nint_cart
      external int_nint_cart
c::   passed
c:tex-\begin{verbatim}
      integer brain              !< [Input]  basis set handle for bra 
      integer ketin              !< [Input]  basis set handle for ket 
      integer num_q              !< [Input]  number of quartets input
      integer icl(num_q)         !< [Input]  i-contraction quartet labels
      integer jcl(num_q)         !< [Input]  j-contraction quartet labels
      integer kcl(num_q)         !< [Input]  k-contraction quartet labels
      integer lcl(num_q)         !< [input]  l-contraction quartet labels
      logical use_q4             !< [Input]  logical for use of q4
      double precision q4(num_q) !< [Input] symmetry prefactors for ints
*     
      integer l_erilab           !< [Input]  length of eri and label arrays
      integer ilab(l_erilab)     !< [Output] integral labels for ``i'' 
      integer jlab(l_erilab)     !< [Output] integral labels for ``j'' 
      integer klab(l_erilab)     !< [Output] integral labels for ``k'' 
      integer llab(l_erilab)     !< [Output] integral labels for ``l'' 
      integer nints              !< [Output] number of integrals returned
      double precision eri(l_erilab) !< [Output] integral values
      integer lscr               !< [Input] length of scratch array
      double precision scr(lscr) !< [Input] scratch array
      double precision zerotol   !< [Input] zero integral threshold
      logical canonicalize       !< [Input] Canonicalize integral labels?
c:tex-\end{verbatim}
c::   local
      integer ip_int, leri_left
      integer numint, nint_for_q
#if defined(DEBUG_TXS)
      integer numinti, numintj, numintk, numintl
#endif
      integer i_q, i_q_done
c     
#include "bas_ibs_sfn.fh"
c     
      intb_nw_2e4c = .true.
c     
      ip_int = 1
      leri_left = l_erilab
      i_q_done = num_q_start - 1
#if defined(DEBUG_TXS)
      write(luout,*)' num_q_start ',num_q_start
      write(luout,*)' num_q_last  ',num_q_last
#endif
      if (num_q_start.eq.num_q_last) then
         write(luout,*)' l_erilab: too small ',l_erilab
         call errquit(' intb_nw_2e4c: infinite blocking loop ',911,
     &          INT_ERR)
      endif
      num_q_last = num_q_start
      nints = 0
      do 00100 i_q = num_q_start, num_q
         numint=int_nint_cart(
     &       brain,icl(i_q),brain,jcl(i_q),
     &       ketin,kcl(i_q),ketin,lcl(i_q))
         if (numint.gt.leri_left) then
#if defined(DEBUG_TXS)
           numinti=int_nint_cart(brain,icl(i_q),0,0,0,0,0,0)
           numintj=int_nint_cart(brain,jcl(i_q),0,0,0,0,0,0)
           numintk=int_nint_cart(ketin,kcl(i_q),0,0,0,0,0,0)
           numintl=int_nint_cart(ketin,lcl(i_q),0,0,0,0,0,0)
           write(luout,*)' shells ',icl(i_q),jcl(i_q),kcl(i_q),lcl(i_q)
           write(luout,*)' numinti  ',numinti
           write(luout,*)' numintj  ',numintj
           write(luout,*)' numintk  ',numintk
           write(luout,*)' numintl  ',numintl
           write(luout,*)' numint   ',numint
           write(luout,*)' leri_left',leri_left
#endif
            goto 00002
         else
            call int_l2e4c(brain, icl(i_q), jcl(i_q),
     $           ketin, kcl(i_q), lcl(i_q),
     $           zerotol, canonicalize, leri_left, eri(ip_int),
     $           nint_for_q,
     $           ilab(ip_int), jlab(ip_int), klab(ip_int), llab(ip_int),
     $           lscr, scr)
            if (use_q4 .and. nint_for_q.gt.0)
     $           call dscal(nint_for_q,q4(i_q),eri(ip_int),1)
            i_q_done = i_q_done + 1
#if defined(DEBUG_TXS)
            write(luout,*)' numinti  ',numinti
            write(luout,*)' numintj  ',numintj
            write(luout,*)' numintk  ',numintk
            write(luout,*)' numintl  ',numintl
            write(luout,*)' i_q_done = ',i_q_done,
     &          '   shells ',icl(i_q),jcl(i_q),kcl(i_q),lcl(i_q),
     &          '   numint ',numint,' nint_for_q ',nint_for_q
#endif
            ip_int = ip_int + nint_for_q
            leri_left = leri_left - nint_for_q
            nints = nints + nint_for_q
         endif
00100 continue
      intb_nw_2e4c = .false.
00002 continue
      num_q_start = i_q_done + 1
#if defined(DEBUG_TXS)
      write(luout,*)' num_q_start on next pass is ',num_q_start
#endif
c     
      end
C> @}
