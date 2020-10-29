
#define NBLOCKS 4


*
*     ***********************************
*     *					*
*     *	        C1dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine C1dB_SumAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "C1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j
    

      if (np_j.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_j(),
     >                       9+MSGDBL,sum,1,'+')
      end if

      return
      end


*     ***********************************
*     *					*
*     *	        C1dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine C1dB_ISumAll(sum)
c     implicit none
      integer  sum


#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "C1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j



      if (np_j.gt.1) then
        call GA_PGROUP_IGOP(Parallel3d_comm_j(),
     >                      9+MSGINT,sum,1,'+')
      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         C1dB_MaxAll             *
*     *                                 *
*     ***********************************
      subroutine C1dB_MaxAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "C1dB.fh"

      integer msglen,mpierr,np
      real*8  sumall

*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j

      if (np_j.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_j(),
     >                       9+MSGDBL,sum,1,'max')
      end if


      return
      end



*     ***********************************
*     *					*
*     *	        C1dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine C1dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "bafdecls.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"
#include "C1dB.fh"



*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j


      if (np_j.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_j(),
     >                       9+MSGDBL,sum,n,'+')
      end if
      return
      end


*     ***********************************
*     *					*
*     *	        C1dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine C1dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "C1dB.fh"


*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j


      if (np_j.gt.1) then
        call GA_PGROUP_IGOP(Parallel3d_comm_j(),
     >                      9+MSGINT,sum,n,'+')
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      C1dB_Brdcst_values         *
*     *                                 *
*     ***********************************

      subroutine C1dB_Brdcst_values(psend,nsize,sum)
      implicit none
      integer psend,nsize
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "C1dB.fh"

*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j

      integer np

      call Parallel3d_np_j(np)
      if (np.gt.1) then
         call GA_PGROUP_BRDCST(Parallel3d_comm_j(),
     >                         9+MSGDBL,sum,mdtob(nsize),psend)
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         C1dB_isendrecv          *
*     *                                 *
*     ***********************************
      subroutine C1dB_isendrecv(pto,  ssize,sdata,
     >                          pfrom,rsize,rdata,
     >                          request,reqcnt)
      implicit none
      integer pto,ssize
      real*8  sdata(*)
      integer pfrom,rsize
      real*8  rdata(*)
      integer request(*)
      integer reqcnt

#include "C1dB.fh"


*     **** local variables ****
      integer msgtype,mpierr

*     **** external functions ****
      integer  Parallel3d_comm_j
      external Parallel3d_comm_j


      call errquit(' C1dB_isendrecv:not implemented!',0,0)
      return
      end

*     ***********************************
*     *                                 *
*     *         C1dB_WaitAll            *
*     *                                 *
*     ***********************************
*
*  This routine waits for the sends and receives to 
* finish that were started with C1dB_isendrecv
* 
      subroutine C1dB_WaitAll(request,reqcnt)
      implicit none
      integer request(*)
      integer reqcnt

#include "bafdecls.fh"
#include "errquit.fh"

      call errquit(' C1dB_WaitAll:not implemented!',0,0)
      return
      end

c     ****************************************
c     *                                      *
c     *          C1dB_Brdcst_step            *
c     *                                      *
c     ****************************************
c
c  This routine performs step l of a butterfly Broadcast all algorithm. The step
c  l spans from 0..(Level-1) where the number of levels is Level = Log(np_j)/Log(2).
c
c   Entry - l: Butterfly step 0...(Level-1)
c           na: an array of length np_j containing the number of orbitals per taskid_j
c           blocks0: number of blocks to send size=blocks0, 
c                    the exceptions are:
c                    if blocks0==0: the block size is size=2**l.
c                    if blocks0==-1: block size is size=(np_j-2**Level)/2 + 1 for l==(Level-1),
c                                    blocksize is size=2**l otherwise
c           n2ft3d: leading size of psi_rep
c           psi_rep: data array
c   Exit -
c           psi_rep: modified data array 
c           requests,reqcnt: tags for asychronous message passing
c
      subroutine C1dB_Brdcst_step(l,na,blocks0,
     >                            n2ft3d,psi_rep,
     >                            requests,reqcnt)
      implicit none
      integer l,na(*),blocks0
      integer n2ft3d
      real*8 psi_rep(n2ft3d,*)
      integer requests(*),reqcnt

*     *** local variables ***
      integer taskid_j,np_j
      integer i,pr,ps,shift,size,Level
      integer pto,pfrom,rsize,ssize,rindx,sindx

*     *** local variables ***
      integer  Butter_levels
      external Butter_levels

      call errquit(' C1dB_Brdcst_step:not implemented!',0,0)
      return
      end


c     ****************************************
c     *                                      *
c     *           C1dB_Reduce_step           *
c     *                                      *
c     ****************************************
c
c  This routine performs step l of a butterfly Reduceall algorithm. The step
c  l spans from 0..(Level-1) where the number of levels is Level = Log(np_j)/Log(2).
c
c   Entry - l: Butterfly step 0...(Level-1)
c           na: an array of length np_j containing the number of orbitals per taskid_j
c           blocks0: number of blocks to send size=blocks0, 
c                    the exceptions are:
c                    if blocks0==0: the block size is size=2**l.
c                    if blocks0==-1: block size is size=(np_j-2**Level)/2 + 1 for l==(Level-1),
c                                    blocksize is size=2**l otherwise
c           n2ft3d: leading size of psi_rep
c           hpsi_rep: data array
c           tmp: tempory data array. Needs to be at least  n2ft3d*size
c
c   Exit - hpsi_rep: modified data array
c          requests,reqcnt: tags for asychronous message passing
c
      subroutine C1dB_Reduce_step(l,na,blocks0,
     >                            n2ft3d,hpsi_rep,tmp)
      implicit none
      integer l,na(*),blocks0
      integer n2ft3d
      real*8 hpsi_rep(n2ft3d,*)
      real*8 tmp(*)

*     *** local variables ***
      integer taskid_j,np_j
      integer i,pr,ps,size,shift,Level,pfrom,pto
      integer rsize,ssize,rindx,sindx
      integer requests(10),reqcnt

*     *** local variables ***
      integer  Butter_levels
      external Butter_levels

      call errquit(' C1dB_Reduce_step:not implemented!',0,0)
      return
      end
c $Id$
