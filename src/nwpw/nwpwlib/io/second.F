#define TCGMSG
*
* $Id: second.F,v 1.1 2001-08-30 01:56:20 bylaska Exp $
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
