C $Id$
*
C> \ingroup nwint
C> @{
C> 
C> \brief Activates the standard integral accuracy settings
C>
C> Different levels of integral accuracy are set up in 
C> `int_init`. This routine activates the accuracy settings for
C> "standard integral accuracy". The specifics of the integral
C> accuracy are picked up from a common block and then transferred
C> to the internals of the various integral packages.
C>
c:tex-% part of API Internal Routines
c:tex-\subsection{int\_acc\_std}
c:tex-This routine sets the integral threshold for radial cutoffs in all 
c:tex-integral codes used in the api via a parameter statement.  Other
c:tex-routines have access via the apiP.fh common blocks and the set/get API.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_acc_std()
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c:tex-The default ''standard'' value for the integral accuracy is:
c:tex-\begin{verbatim}
c      parameter(val_def = 1.0d-12)
c:tex-\end{verbatim}
c
      val_int_acc = val_int_acc_std
      call int_sp_acc(val_int_acc)
      call texas_set_accy(val_int_acc)
      call hnd_set_acc(val_int_acc)
      end
c......................................................................
C> 
C> \brief Activates the high integral accuracy settings
C>
C> Different levels of integral accuracy are set up in 
C> `int_init`. This routine activates the accuracy settings for
C> "high integral accuracy". The specifics of the integral
C> accuracy are picked up from a common block and then transferred
C> to the internals of the various integral packages.
C>
c:tex-% part of API Internal Routines
c:tex-\subsection{int\_acc\_high}
c:tex-This routine sets the integral threshold to ``high'' accuracy 
c:tex-for radial cutoffs in all integral codes used in the api via a 
c:tex-parameter statement.  Other routines have access via the apiP.fh 
c:tex-common blocks and the set/get API.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_acc_high()
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c:tex-The default ''high accuracy'' value for the integral API is:
c:tex-\begin{verbatim}
c      parameter (val_def_high = 1.0d-30)
c:tex-\end{verbatim}
c:tex-This is needed for certain algorithms within NWChem, e.g., 
c:tex-the Schwarz inequality screening.
c:tex-
      val_int_acc = val_int_acc_high
      call int_sp_acc(val_int_acc)
      call texas_set_accy(val_int_acc)
      call hnd_set_acc(val_int_acc)
      end
c......................................................................
C>
C> \brief Retrieve the current integral accuracy setting
C>
C> Retrieves the currently active integral accuracy setting and 
C> returns that value.
C>
c:tex-% part of API Internal Routines
c:tex-\subsection{int\_acc\_get}
c:tex-This routine returns the current integral threshold 
c:tex-for radial cutoffs in all integral codes used in the api via a 
c:tex-parameter statement.  
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_acc_get(retval)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c:tex-\begin{verbatim}
      double precision retval !< [Output] current threshold
c:tex-\end{verbatim}
      retval = val_int_acc
      end
c......................................................................
C>
C> \brief Sets the current integral accuracy setting
C>
C> Sets the currently active integral accuracy setting to the value 
C> specified by the argument.
C>
c:tex-% part of API Internal Routines
c:tex-\subsection{int\_acc\_set}
c:tex-This routine sets the current integral threshold 
c:tex-for radial cutoffs in all integral codes used in the api via a 
c:tex-parameter statement.  
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_acc_set(setval)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
c:tex-\begin{verbatim}
      double precision setval !< [Input] new threshold
c:tex-\end{verbatim}
      val_int_acc = setval
      call int_sp_acc(val_int_acc)
      call texas_set_accy(val_int_acc)
      call hnd_set_acc(val_int_acc)
      end
C> @}
