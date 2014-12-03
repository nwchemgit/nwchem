

*
*     ***********************************
*     *					*
*     *	        K1dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine K1dB_SumAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "K1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
    
      if (np_k.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_k(),
     >                       9+MSGDBL,sum,1,'+')
      end if

      return
      end


*     ***********************************
*     *					*
*     *	        K1dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine K1dB_ISumAll(sum)
c     implicit none
      integer  sum


#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "K1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k



      if (np_k.gt.1) then
        call GA_PGROUP_IGOP(Parallel3d_comm_k(),
     >                      9+MSGINT,sum,1,'+')
      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         K1dB_MaxAll             *
*     *                                 *
*     ***********************************
      subroutine K1dB_MaxAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "K1dB.fh"

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

      if (np_k.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_k(),
     >                       9+MSGDBL,sum,1,'max')
      end if
      return
      end

*     ***********************************
*     *                                 *
*     *         K1dB_MinAll             *
*     *                                 *
*     ***********************************
      subroutine K1dB_MinAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "K1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

      if (np_k.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_k(),
     >                       9+MSGDBL,sum,1,'min')
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        K1dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine K1dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "bafdecls.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"
#include "K1dB.fh"



*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k


      if (np_k.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_k(),
     >                       9+MSGDBL,sum,n,'+')
      end if
      return
      end


*     ***********************************
*     *					*
*     *	        K1dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine K1dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "K1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k


      if (np_k.gt.1) then
        call GA_PGROUP_IGOP(Parallel3d_comm_k(),
     >                      9+MSGINT,sum,n,'+')
      end if

      return
      end


c $Id$
