c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Check to see if a given shell is valid
C>
C> \return Return .true. if the shell is valid, and .false. otherwise
C>
C> This function checks to see that the given shell is valid
C> This subroutine call can be replaced by a statement 
C> function sequence:
C> \code
C> #include "basP.fh"
C> #include "geobasmapP.fh"
C>       logical inline_chk_sh
C>       inline_chk_sh(bra,ish) =
C>      $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(bra)))
C> \endcode
C> Where bra is the lexical basis index (not the handle).
C> or you could use the following with the handle.
C> \code
C>       inline_chk_sh(bra,ish) =
C>      $     ((ish.gt.0) .and. 
C>      $(ish.le.ncont_tot_gb(bra+Basis_Handle_Offset)))
C> \endcode
C>
c:tex-% part of the internal API routines
c:tex-\subsection{int\_chk\_sh}
c:tex-This function checks to see that the given shell is valid
c:tex-Returns .true. if so else returns .false. if not.
c:tex-This subroutine call can be replaced by a statement 
c:tex-function sequence:
c:tex-\begin{verbatim}
c:tex-#include "basP.fh"
c:tex-#include "geobasmapP.fh"
c:tex-      logical inline_chk_sh
c:tex-      inline_chk_sh(bra,ish) =
c:tex-     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(bra)))
c:tex-\end{verbatim}
c:tex-Where bra is the lexical basis index (not the handle).
c:tex-Or you could use the following with the handle:
c:tex-\begin{verbatim}
c:tex-      inline_chk_sh(bra,ish) =
c:tex-     $     ((ish.gt.0) .and. 
c:tex-     $(ish.le.ncont_tot_gb(bra+Basis_Handle_Offset)))
c:tex-\end{verbatim}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function int_chk_sh(basisin,shell)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin !< [Input] basis set handle
      integer shell   !< [Input] lexical shell index
c:tex-\end{verbatim}
c::local
      integer basis
c
      int_chk_sh = .true.
#ifdef API_DEBUG
      int_chk_sh = bas_check_handle(basisin,'int_chk_sh')
      if (.not.int_chk_sh) return
#endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
      int_chk_sh = shell.gt.0 .and.
     &             shell.le.ncont_tot_gb(basis)
      if (.not.int_chk_sh) then
        write(luout,*)'invalid contraction/shell index'
        write(luout,*)' shell range is: <1:',ncont_tot_gb(basis),'>'
        write(luout,*)' input shell is: ',shell
      endif
      end
C> @}
