c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Check that the integral API is initialized
C>
C> Check that the integral API is initialized. If it is not 
C> then a warning is constructed from the message information
C> in the argument and printed on the output. The message text
C> usually is chosen to identify the calling routine.
C>
C> \return Return .true. if the API is initialized, and .false. 
C> otherwise.
C>
c:tex-% part of the internal API routines
c:tex-\subsection{int\_chk\_init}
c:tex-This function checks to see that the integral API is initialized.
c:tex-Returns .true. if initialized and .false. if not.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function int_chk_init(msg)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "stdio.fh"
c:: passed
c:tex-\begin{verbatim}
      character*(*) msg !< [Input] The message text. 
c:tex-\end{verbatim}
c
c check initialization
c
      int_chk_init = init_int.ne.0
      if (.not.int_chk_init) then
        write(luout,*)' ERROR:: <',msg,'>'
        write(luout,*)' int_init was not called '
      endif
      end
C> @}
