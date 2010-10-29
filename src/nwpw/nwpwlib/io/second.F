#define TCGMSG
*
* $Id$
*
#ifdef MPI
      SUBROUTINE current_second(T)
      implicit none
      real*8 T
      real*8 MPI_Wtime
      external MPI_Wtime

      T=MPI_Wtime()

      RETURN
      END
#endif

#ifdef TCGMSG
      SUBROUTINE current_second(T)
      implicit none
      real*8 T
#include "tcgmsg.fh"

      T=TCGTIME()

      RETURN
      END
#endif
