c $Id: int_chk_init.F,v 1.2 1996-10-07 16:50:20 d3e129 Exp $
*
c:tex-% part of the internal API routines
c:tex-\subsubsection{int\_chk\_init}
c:tex-This function checks to see that the integral API is initialized.
c:tex-Returns .true. if initialized and .false. if not.
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function int_chk_init(msg)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c:: passed
c:tex-\begin{verbatim}
      character*(*) msg ! [input] usually indentfy calling routine
c:tex-\end{verbatim}
c
c check initialization
c
      int_chk_init = init_int.ne.0
      if (.not.int_chk_init) then
        write(6,*)' ERROR:: <',msg,'>'
        write(6,*)' int_init was not called '
      endif
      end
