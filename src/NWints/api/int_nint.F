c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Calculate the number of integrals for a given
C> shell quartet
C>
C> This is an internal routine which should never be called by an
C> NWChem application module directly.
C>
C> This routine computes the number of integrals for a given 
C> shell/contraction grouping; if an input shell is zero then 
C> the routine ignores this shell.  This routine works
C> for both cartesian and spherical basis sets.
C>
c:tex-% this is an internal API routine
c:tex-\subsection{int\_nint}
c:tex-This routine computes the number of integrals for a given 
c:tex-shell/contraction grouping; if an input shell is zero then 
c:tex-the routine ignores this shell.  This routine will work 
c:tex-for both cartesian and spherical basis sets.
c:tex-This routine should never be called by an NWChem 
c:tex-application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      integer function int_nint(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
c:tex-\end{verbatim}
      implicit none
c
c
#include "bas.fh"
#include "errquit.fh"
c::passed
c:tex-\begin{verbatim}
      integer ibasin   !< [Input] basis set handle for icnt
      integer icnt     !< [Input] contraction index (e.g., ish)
      integer jbasin   !< [Input] basis set handle for jcnt
      integer jcnt     !< [Input] contraction index (e.g., jsh)
      integer kbasin   !< [Input] basis set handle for kcnt
      integer kcnt     !< [Input] contraction index (e.g., ksh)
      integer lbasin   !< [Input] basis set handle for lcnt
      integer lcnt     !< [Input] contraction index (e.g., lsh)
c:tex-\end{verbatim}
c::local
      integer hi,lo,nbf
c
      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_nint = 0
        return
      endif
c
      int_nint = 1
c
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_cn2bfr(ibasin,icnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911, 
     &           INT_ERR)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_cn2bfr(jbasin,jcnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911,
     &         INT_ERR)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_cn2bfr(kbasin,kcnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911, INT_ERR)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_cn2bfr(lbasin,lcnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911, INT_ERR)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
      end
C>
C> \brief Calculate the number of integrals for a given
C> quartet of unique shells
C>
C> This is an internal routine which should never be called by an
C> NWChem application module directly.
C>
C> This routine computes the number of integrals for a given 
C> shell/contraction grouping; if an input shell is zero then 
C> the routine ignores this shell. The input shell must be a
C> unique shell in the sense of the basis set API. This routine
C> works for both cartesian and spherical basis sets.
C>
c:tex-% this is an internal API routine
c:tex-\subsection{int\_unint}
c:tex-This routine computes the number of integrals for a given 
c:tex-shell/contraction grouping; if an input shell is zero then 
c:tex-the routine ignores this shell.  The input shell must be a 
c:tex-unique shell in the sense of the basis set API.  This 
c:tex-routine will work for both cartesian and spherical basis sets.
c:tex-This routine should never be called by an NWChem 
c:tex-application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      integer function int_unint(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
c:tex-\end{verbatim}
      implicit none
c
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
c::passed
c:tex-\begin{verbatim}
      integer ibasin   !< [Input] basis set handle for icnt
      integer icnt     !< [Input] unique contraction index (e.g., ish)
      integer jbasin   !< [Input] basis set handle for jcnt
      integer jcnt     !< [Input] unique contraction index (e.g., jsh)
      integer kbasin   !< [Input] basis set handle for kcnt
      integer kcnt     !< [Input] unique contraction index (e.g., ksh)
      integer lbasin   !< [Input] basis set handle for lcnt
      integer lcnt     !< [Input] unique contraction index (e.g., lsh)
c:tex-\end{verbatim}
c::local
      integer nbf
c
      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_unint = 0
        return
      endif
c
      int_unint = 1
c
c..icnt/ish
      if (icnt.gt.0) then
        nbf = nbf_from_ucont(icnt,ibasin)
        int_unint = int_unint*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        nbf = nbf_from_ucont(jcnt,jbasin)
        int_unint = int_unint*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        nbf = nbf_from_ucont(kcnt,kbasin)
        int_unint = int_unint*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        nbf = nbf_from_ucont(lcnt,lbasin)
        int_unint = int_unint*nbf
      endif
      end
C>
C> \brief Calculate the number of cartesian integrals for a given
C> shell quartet
C>
C> This is an internal routine which should never be called by an
C> NWChem application module directly.
C>
C> This routine computes the number of integrals for a given 
C> shell/contraction grouping; if an input shell is zero then 
C> the routine ignores this shell.  This routine will work 
C> for both cartesian and spherical basis sets. However, in either
C> case it will calculate the number of integrals as if the basis sets
C> were cartesian. The reason for this is that in the case of
C> spherical harmonic basis sets the integral codes calculate the
C> cartesian integrals first and transform them to spherical harmonics.
C> Because of this all buffers have to be large enough to hold the 
C> larger sets of cartesian integrals that are calculated as 
C> intermediates.
C>
c:tex-% this is an internal API routine
c:tex-\subsection{int\_nint\_cart}
c:tex-This routine computes the number of integrals for a given 
c:tex-shell/contraction grouping; if an input shell is zero then 
c:tex-the routine ignores this shell.  This 
c:tex-routine will work for both cartesian and spherical basis 
c:tex-sets, but {\it returns the cartesian size} (this is how the 
c:tex-integrals are computed!).
c:tex-This routine should never be called by an NWChem 
c:tex-application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      integer function int_nint_cart(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
c
c::passed
c:tex-\begin{verbatim}
      integer ibasin   !< [Input] basis set handle for icnt
      integer icnt     !< [Input] contraction index (e.g., ish)
      integer jbasin   !< [Input] basis set handle for jcnt
      integer jcnt     !< [Input] contraction index (e.g., jsh)
      integer kbasin   !< [Input] basis set handle for kcnt
      integer kcnt     !< [Input] contraction index (e.g., ksh)
      integer lbasin   !< [Input] basis set handle for lcnt
      integer lcnt     !< [Input] contraction index (e.g., lsh)
c:tex-\end{verbatim}
c:local
      integer type, nprim, ngen, spcart
      integer nbf

      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_nint_cart = 0
        return
      endif

      int_nint_cart = 1

      
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_continfo(ibasin,icnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911, 
     &              INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_continfo(jbasin,jcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911, 
     &               INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_continfo(kbasin,kcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911,
     &              INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_continfo(lbasin,lcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911, 
     &                 INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
      end
C>
C> \brief Calculate the number of cartesian integrals for a given
C> quartet of unique shells
C>
C> This is an internal routine which should never be called by an
C> NWChem application module directly.
C>
C> This routine computes the number of integrals for a given 
C> shell/contraction grouping; if an input shell is zero then 
C> the routine ignores this shell. The shells have to be unique
C> shells in the sense of the basis set API. This routine will work 
C> for both cartesian and spherical basis sets. However, in either
C> case it will calculate the number of integrals as if the basis sets
C> were cartesian. The reason for this is that in the case of
C> spherical harmonic basis sets the integral codes calculate the
C> cartesian integrals first and transform them to spherical harmonics.
C> Because of this all buffers have to be large enough to hold the 
C> larger sets of cartesian integrals that are calculated as 
C> intermediates.
C>
c:tex-% this is an internal API routine
c:tex-\subsection{int\_unint\_cart}
c:tex-This routine computes the number of integrals for a given 
c:tex-shell/contraction grouping; if an input shell is zero then 
c:tex-the routine ignores this shell.  The input shell must be a 
c:tex-unique shell in the sense of the basis set API.  This 
c:tex-routine will work for both cartesian and spherical basis 
c:tex-sets, but {\it returns the cartesian size} (this is how the 
c:tex-integrals are computed!).
c:tex-This routine should never be called by an NWChem 
c:tex-application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      integer function int_unint_cart(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
c
c::passed
c:tex-\begin{verbatim}
      integer ibasin   !< [Input] basis set handle for icnt
      integer icnt     !< [Input] unique contraction index (e.g., ish)
      integer jbasin   !< [Input] basis set handle for jcnt
      integer jcnt     !< [Input] unique contraction index (e.g., jsh)
      integer kbasin   !< [Input] basis set handle for kcnt
      integer kcnt     !< [Input] unique contraction index (e.g., ksh)
      integer lbasin   !< [Input] basis set handle for lcnt
      integer lcnt     !< [Input] unique contraction index (e.g., lsh)
c:tex-\end{verbatim}
c::local
      integer type, nprim, ngen, spcart
      integer nbf
c::functions
      logical  bas_ucontinfo
      external bas_ucontinfo

      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_unint_cart = 0
        return
      endif

      int_unint_cart = 1

      
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_ucontinfo(ibasin,icnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911,
     &           INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_ucontinfo(jbasin,jcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911,
     &            INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_ucontinfo(kbasin,kcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911,
     &            INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_ucontinfo(lbasin,lcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911,
     &          INT_ERR)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
      end
* not used :      integer function int_nint_raw(ibasin,icnt,jbasin,jcnt,
* not used :     &       kbasin,kcnt,lbasin,lcnt)
* not used :      implicit none
* not used :#include "bas.fh"
* not used :#include "nwc_const.fh"
* not used :#include "int_nbf.fh"
* not used :c
* not used :c::passed
* not used :      integer ibasin   ! [input] basis set handle for icnt
* not used :      integer icnt     ! [input] contraction index (e.g., ish)
* not used :      integer jbasin   ! [input] basis set handle for jcnt
* not used :      integer jcnt     ! [input] contraction index (e.g., jsh)
* not used :      integer kbasin   ! [input] basis set handle for kcnt
* not used :      integer kcnt     ! [input] contraction index (e.g., ksh)
* not used :      integer lbasin   ! [input] basis set handle for lcnt
* not used :      integer lcnt     ! [input] contraction index (e.g., lsh)
* not used :c:local
* not used :      integer type, nprim, ngen, spcart
* not used :      integer nbf
* not used :
* not used :      if (
* not used :     &    icnt.eq.jcnt.and.
* not used :     &    jcnt.eq.kcnt.and.
* not used :     &    kcnt.eq.lcnt.and.
* not used :     &    icnt.eq.0
* not used :     &   ) then
* not used :        int_nint_raw = 0
* not used :        return
* not used :      endif
* not used :
* not used :      int_nint_raw = 1
* not used :
* not used :      
* not used :c..icnt/ish
* not used :      if (icnt.gt.0) then
* not used :        if (.not.bas_continfo(ibasin,icnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_nint_raw: bas_continfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_nint_raw = int_nint_raw*nbf
* not used :      endif
* not used :c
* not used :c..jcnt/jsh
* not used :      if (jcnt.gt.0) then
* not used :        if (.not.bas_continfo(jbasin,jcnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_nint_raw: bas_continfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_nint_raw = int_nint_raw*nbf
* not used :      endif
* not used :c
* not used :c..kcnt/ksh
* not used :      if (kcnt.gt.0) then
* not used :        if (.not.bas_continfo(kbasin,kcnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_nint_raw: bas_continfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_nint_raw = int_nint_raw*nbf
* not used :      endif
* not used :c
* not used :c..lcnt/lsh
* not used :      if (lcnt.gt.0) then
* not used :        if (.not.bas_continfo(lbasin,lcnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_nint_raw: bas_continfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_nint_raw = int_nint_raw*nbf
* not used :      endif
* not used :      end
* not used :      integer function int_unint_raw(ibasin,icnt,jbasin,jcnt,
* not used :     &       kbasin,kcnt,lbasin,lcnt)
* not used :      implicit none
* not used :#include "nwc_const.fh"
* not used :#include "int_nbf.fh"
* not used :c
* not used :c::passed
* not used :      integer ibasin   ! [input] basis set handle for icnt
* not used :      integer icnt     ! [input] contraction index (e.g., ish)
* not used :      integer jbasin   ! [input] basis set handle for jcnt
* not used :      integer jcnt     ! [input] contraction index (e.g., jsh)
* not used :      integer kbasin   ! [input] basis set handle for kcnt
* not used :      integer kcnt     ! [input] contraction index (e.g., ksh)
* not used :      integer lbasin   ! [input] basis set handle for lcnt
* not used :      integer lcnt     ! [input] contraction index (e.g., lsh)
* not used :c::local
* not used :      integer type, nprim, ngen, spcart
* not used :      integer nbf
* not used :c::functions
* not used :      logical  bas_ucontinfo
* not used :      external bas_ucontinfo
* not used :
* not used :      if (
* not used :     &    icnt.eq.jcnt.and.
* not used :     &    jcnt.eq.kcnt.and.
* not used :     &    kcnt.eq.lcnt.and.
* not used :     &    icnt.eq.0
* not used :     &   ) then
* not used :        int_unint_raw = 0
* not used :        return
* not used :      endif
* not used :
* not used :      int_unint_raw = 1
* not used :
* not used :      
* not used :c..icnt/ish
* not used :      if (icnt.gt.0) then
* not used :        if (.not.bas_ucontinfo(ibasin,icnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_unint_raw = int_unint_raw*nbf
* not used :      endif
* not used :c
* not used :c..jcnt/jsh
* not used :      if (jcnt.gt.0) then
* not used :        if (.not.bas_ucontinfo(jbasin,jcnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_unint_raw = int_unint_raw*nbf
* not used :      endif
* not used :c
* not used :c..kcnt/ksh
* not used :      if (kcnt.gt.0) then
* not used :        if (.not.bas_ucontinfo(kbasin,kcnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_unint_raw = int_unint_raw*nbf
* not used :      endif
* not used :c
* not used :c..lcnt/lsh
* not used :      if (lcnt.gt.0) then
* not used :        if (.not.bas_ucontinfo(lbasin,lcnt,type,nprim,ngen,spcart))
* not used :     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
* not used :        nbf = int_nbf_x(type)
* not used :        int_unint_raw = int_unint_raw*nbf
* not used :      endif
* not used :      end
C> @}
