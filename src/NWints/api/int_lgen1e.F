c $Id$
*
C> \ingroup nwint
C>
C> \brief Generate labels for general 2 index 1 electron integrals
C>
C> This routine generates labels for general 2 index one 
C> electron integrals.  This is mostly unused since the other 
C> integral type specific label routines are now used.  This 
C> routine requires that the integral block be computed prior 
C> to the label call.  Other routines now integrate label 
C> generation with computation.  
C>
c:tex-% this is part of the API Standard Integral routines
c:tex-\subsection{int\_l1gen1e}
c:tex-This routine generates labels for general 2 index one 
c:tex-electron integrals.  This is mostly unused since the other 
c:tex-integral type specific label routines are now used.  This 
c:tex-routine requires that the integral block be computed prior 
c:tex-to the label call.  Other routines now integrate label 
c:tex-generation with computation.  
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_lgen1e(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, l1e, Gen1e, lscr, scr, numgen) 
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
c dcopy
c errquit
c::function
      logical int_chk_sh
      external int_chk_sh
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [input] bra basis set handle
      integer ish     !< [input] bra shell lexical index
      integer j_basis !< [input] ket basis set handle
      integer jsh     !< [input] ket shell lexical index
      double precision zerotol !< [input] zero threshold
      integer l1e       !< [input] length of buffers for integrals
      integer ilab(l1e) !< [output] i bas func labels array
      integer jlab(l1e) !< [output] j bas func labels array
      double precision Gen1e(l1e) !< [input/output] 1e integrals
      integer lscr !< [input] length of scratch array
      double precision scr(lscr) !< [scratch] array
      integer numgen !< [output] number of integrals 
                     !< saved and returned
c:tex-\end{verbatim}
c::local
      integer ibas, jbas, icount, i, j
      logical ijbas
      logical shells_ok
      logical nonzero
c
#include "bas_ibs_sfn.fh"
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_lgen1e: invalid contraction/shell',0, 
     &            BASIS_ERR)
c      
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
c copy Gen1e vector to scratch array
c      
      if (l1e.gt.lscr)
     & call errquit('int_lgen1e: scratch array too small. need:',l1e,
     &        INT_ERR)
      
      call dcopy(l1e,Gen1e,1,scr,1)
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
              Gen1e(numgen) = scr(icount)
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
