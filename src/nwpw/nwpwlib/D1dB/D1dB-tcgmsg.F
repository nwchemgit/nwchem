
#define NBLOCKS 4


*
*     ***********************************
*     *					*
*     *	        D1dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine D1dB_SumAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "D1dB.fh"


*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j
    

      if (np_j.gt.1) then
         call GA_PGROUP_DGOP(Parallel2d_comm_j(),
     >                       9+MSGDBL,sum,1,'+')
      end if

      return
      end


*     ***********************************
*     *					*
*     *	        D1dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine D1dB_ISumAll(sum)
c     implicit none
      integer  sum


#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "D1dB.fh"


*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j



      if (np_j.gt.1) then
        call GA_PGROUP_IGOP(Parallel2d_comm_j(),
     >                      9+MSGINT,sum,1,'+')
      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         D1dB_MaxAll             *
*     *                                 *
*     ***********************************
      subroutine D1dB_MaxAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "D1dB.fh"

      integer msglen,mpierr,np
      real*8  sumall

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      if (np_j.gt.1) then
         call GA_PGROUP_DGOP(Parallel2d_comm_j(),
     >                       9+MSGDBL,sum,1,'max')
      end if


      return
      end



*     ***********************************
*     *					*
*     *	        D1dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine D1dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "mafdecls.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"
#include "D1dB.fh"



*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j


      if (np_j.gt.1) then
         call GA_PGROUP_DGOP(Parallel2d_comm_j(),
     >                       9+MSGDBL,sum,n,'+')
      end if
      return
      end


*     ***********************************
*     *					*
*     *	        D1dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine D1dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "D1dB.fh"


*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j


      if (np_j.gt.1) then
        call GA_PGROUP_IGOP(Parallel2d_comm_j(),
     >                      9+MSGINT,sum,n,'+')
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      D1dB_Brdcst_values         *
*     *                                 *
*     ***********************************

      subroutine D1dB_Brdcst_values(psend,nsize,sum)
      implicit none
      integer psend,nsize
      real*8  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "D1dB.fh"

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      integer np

      call Parallel2d_np_j(np)
      if (np.gt.1) then
         call GA_PGROUP_BRDCST(Parallel2d_comm_j(),
     >                         9+MSGDBL,sum,mdtob(nsize),psend)
      end if

      return
      end

