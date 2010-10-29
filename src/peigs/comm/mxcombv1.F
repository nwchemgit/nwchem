*
* $Id$
*
*======================================================================
*
* DISCLAIMER
*
* This material was prepared as an account of work sponsored by an
* agency of the United States Government.  Neither the United States
* Government nor the United States Department of Energy, nor Battelle,
* nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
* ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
* COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
* SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
* INFRINGE PRIVATELY OWNED RIGHTS.
*
* ACKNOWLEDGMENT
*
* This software and its documentation were produced with Government
* support under Contract Number DE-AC06-76RLO-1830 awarded by the United
* States Department of Energy.  The Government retains a paid-up
* non-exclusive, irrevocable worldwide license to reproduce, prepare
* derivative works, perform publicly and display publicly by or for the
* Government, including the right to distribute to other Government
* contractors.
*
*======================================================================
*
*  -- PEIGS  routine (version 2.1) --
*     Pacific Northwest Laboratory
*     July 28, 1995
*
*======================================================================
	integer function mxlbuf( )
c
c     mxlbuf -- returns the value of PARAMETER BUFSIZ used in MXCOMBV1.
c               The user must pass mxcombv1 a data buffer, datbuf, of length at
c               least mxlbuf() = bufsiz bytes.  This routine allows one to
c               find out how much workspace is required.
c
#include "cmbbfr1.h"
        mxlbuf = bufsiz
        return
        end
	integer function mxcombv1(data,ifuncv,iteml,
     $       nitem,nprocs,nodel,
     &       itype, datbuf)
c
c     mxcombv -- implement combining on a specified set of processors
c                in a portable fashion, using a vectorizable combining
c                function.
c
c     Function:  This routine is identical to mxcomb, except for using
c                a vectorizable combining function (i.e., one that
c                combines two vectors instead of two scalars).
c
c     Arguments:
c         data        (typeless in/out vector) data to be combined.
c                        Length of the data vector is iteml*nitem bytes.
c         ifuncv      (external integer function) vector combining function:
c                        integer function ifuncv (adata,bdata,iteml,nitems)
c         iteml       (integer input) length, in bytes, of each data item.
c                        At present, this must be a multiple of 4 bytes
c                        and must be no more than 512 bytes.  (The latter
c                        restriction is controlled by the bufsiz
c                        compile-time parameter.)
c         nitem       (integer input) number of data items.
c         nprocs      (integer input) number of processors (length of nodel
c                        vector), or EXPRESS constant ialnod (from
c                        the /xpress/ common block)
c         nodel       (integer input vector) list of processors across
c                        which to combine.  ** This list must be IDENTICAL
c                        on all processors within the list.  Even the list
c                        ordering must be the same on all participating
c                        processors. **
c         itype       (integer input) message type to be used.
c
c         datbuf      (integer workspace) INTEGER workspace of
c                        length at least bufsiz/INTSIZ.  Should be aligned
c                        on an 8-byte boundary for best performance
c                        (see comments in cmbbfr1.h).
c
c
c     Returns:  For now, this version always returns "successful", using
c               the kxcomb convention for that value.  Future versions
c               should return failure codes for incorrect arguments and/or
c               communications errors, but it is unclear exactly when and
c               how such returns should be made.  Cavaet caller!
c
c     Limitations:
c        1. The nodel list must be IDENTICAL on all participating
c           processors.  Deadlock will result if it isn-t.
c
c     Method:  Collapse and multicast using a binary spanning tree
c              on the nodel list as provided by the caller.  No attempt
c              is made to optimize communications locality by reordering
c              the list.
c              
c              The collapse-and-multicast approach guarantees that all 
c              participating nodes will receive bit-wise identical 
c              results, even if -ifuncv- isn-t strictly commutative and
c              associative.  It can be less efficient than other methods,
c              depending on gory details of the interconnect, but is
c              easy to code and should be tolerable under most conditions.
c              In particular, the method is appropriate if the underlying
c              communication system behaves as a (virtual) complete
c              interconnect whose performance is dominated by message
c              startup.
c
c     Revision History:
c        1/10/91, design/skeleton by Rik Littlefield.
c        1/17/91, details/debug by Jack Cravatt & Rik Littlefield.
c        4/27/91, finite-buffer deadlock fix by Rik Littlefield.
c        8/02/94, convert datbuf into an argument by Dave Elwood
c       11/11/94, Modify to deal with INTEGER not of length 4 bytes.
c
c     The -data- argument should really be declared integer*4, since
c     its indexing assumes 4-byte integers, but that syntax is
c     rejected by some compilers (e.g., tfc).
c
       	integer data(*)
	external ifuncv
        integer iteml	
        integer nitem
        integer nprocs
        integer nodel(nprocs)
        integer itype
        integer i, item
c     
c     Local variables -- see code for usage.
c
#include "datasizes.h"
*
      integer mypnum, mylpos
      integer prtnr, prtlpos
      logical rcving, xmting
      integer iteml4
      integer mxfer, axfer
      integer xferlen
      integer ibuf, nbuf
      integer idat
      integer istep
      integer iretrn
      integer istat
	integer env(4)
#include "cmbbfr1.h"
	integer datbuf(*)
c
c     EXPRESS-environment common block that contains certain variables
c     and constants.  In this routine, we need -ialnod-.
c
      common /xpress/ nocare, norder, nonode, ihost, ialnod, ialprc
c
c     If nprocs = ialnod, then we need to know the actual number
c     of processors participating in the combine (-nppart-).
c
      istat = mxpara(env)
      if (nprocs.ne.ialnod) then
         nppart = nprocs
      else
         nppart = env(2)
      endif
c
c     Quick exit for special case of no combining required.
c
      if (nppart.lt.2) go to 60
c
c     Determine our processor number and where we are in the list of
c     processor nodes.  If we-re not in this list, then we-re not
c     participating, so immediately exit.  The basic output from this
c     code is -mylpos-, which is this processor-s position in the list
c     in terms of first element = position 0.  This makes the
c     calculations of partner simpler.
c
      mypnum = env(1)
      if (nprocs.ne.ialnod) then
         do 10 i = 1,nprocs
           if (mypnum .eq. nodel(i)) then
              mylpos = i-1
              go to 20 
           endif
10        continue
         go to 60 
      else
         mylpos = mypnum
      endif
c
c     OK, so this processor is participating.  We now need to do a series
c     of ceil(log2(nproc)) collapse-and-combine phases, to compute the
c     result, followed by a series of expand phases to rebroadcast it.
c
c     Determine the item length in terms of integer*4 data declarations.
c
 20	iteml4 = iteml/NBYTEI
c       mxfer = iteml4*(bufsiz/iteml)
	mxfer = iteml*nitem
c
c     Part one, collapse-and-combine by interpreting the node list as
c     a binary tree stored in nodel.  The variable -istep- goes as
c     1,2,4,8,... while istep < nppart-1
c
      istep = 1
30    if (istep .lt. nppart) then
c        Determing whether we are sending, receiving/combining, or
c        idling on this phase.  If sending or receiving, also determine
c        the list position of our communications partner.
c
         rcving = .false.
         xmting = .false.
c
         if (mod(mylpos,2*istep) .eq. 0) then
            prtlpos = mylpos + istep
            if (prtlpos .lt. nppart) rcving = .true.
         endif
c
         if (mod(mylpos,2*istep) .eq. istep) then
            prtlpos = mylpos - istep
            xmting = .true.
         endif
c
c        Now, if we are not idling, we need to spin down the data
c        buffers doing the communicating and combining.
c
         if (rcving .or. xmting) then
            if (nprocs.ne.ialnod) then
               prtnr = nodel(prtlpos+1)
            else
               prtnr = prtlpos
            endif
            idat = 1 - iteml4
            ibuf = bufsiz
            nbuf = 0
c
c           Explicit handshake with partner.  This is to keep
c           a node from getting inundated with messages that it
c           isn-t prepared to accept & dispose of, in case the
c           number of processors isn-t a power of 2.  E.g., if
c           there are 3 processors, we don-t want processor #2
c           getting ahead of the game and sending a bunch of data
c           to processor #0 at the same time that processor #1 is
c           trying to do the same.  What may happen in that case is
c           that the unaccepted messages from #2 fill up the queues,
c           so that the desired ones from #1 can-t get through, and
c           clunk! we-re deadlocked.  (This problem observed on
c           early release Delta software, underlying problem fixed
c           by subsequent release.)
c
c.            if (rcving) then
c.               istat = mxwrit (0,1,prtnr,itype)
c.            else
c.               istat = mxread (junk,1,prtnr,itype)
c.            endif
c
c..            do 40 item = 1,nitem
              item = 1
   39         if (item.le.nitem) then
c
c              Advance pointers to next data item.
c
               idat = idat + iteml4
               ibuf = ibuf + iteml4
c
c              If buffer empty, refill it. 
c
               if (ibuf.gt.nbuf) then
		   axfer = (nitem*iteml4 + 1) - idat
	           if(mxfer .le. axfer) then
                      xferlen = mxfer*NBYTEI
                   else
                      xferlen = axfer*NBYTEI
                   endif
                   if (rcving) then
                      istat = mxread(datbuf, xferlen, prtnr, itype)
                      ibuf = 1
                      nbuf = xferlen/NBYTEI
                   else
                      istat = mxwrit(data(idat), xferlen, prtnr, itype)
                      ibuf = 1
                      nbuf = xferlen/NBYTEI
                   endif
                endif
c
c              Call combining function on the current data item
c
               ntodo = nbuf / iteml4
c.	write (*,*) -me,ntodo, ibuf, idat = -,mypnum,ntodo,ibuf,idat
          if (rcving) then
	  iretrn = ifuncv(data(idat),datbuf(ibuf),iteml,ntodo)
          endif
          idat = idat + (ntodo-1)*iteml4
          ibuf = ibuf + (ntodo-1)*iteml4
          item = item + (ntodo-1)
c
c              End of do item = 1,nitem
              item = item + 1
              go to 39
              endif
40            continue
c
c        At this point, we-ve processed all the data for this level
c        in the tree, so advance one level and go do it again.
c
c        end if (rcving .or. xmting)
         endif
c
      istep = istep * 2
      go to 30 
c
c     end if (istep.lt.nppart)
      endif
c
c     At this point, we-ve completed the collapse.  The result is
c     sitting in the first processor in nodel, and istep has been
c     increased until it is larger than nppart.  To implement the
c     multicast, we essentially reverse the above communication
c     process.  This time, though, we don-t need to mess around
c     with buffering, and can just send the whole data vector in a
c     single comm call.
c
      istep = istep/2

50    if (istep.ge.1) then
c
c        The xmting/rcving flags are computed as above, but with their
c        sense reversed -- last time-s receiver is this time-s transmitter
c
         rcving = .false.
         xmting = .false.
c
         if (mod(mylpos,2*istep) .eq. 0) then
            prtlpos = mylpos + istep
            if (prtlpos .lt. nppart) xmting = .true.
         endif
c
         if (mod(mylpos,2*istep) .eq. istep) then
            prtlpos = mylpos -istep
            rcving = .true.
         endif
		

         if (nprocs.ne.ialnod) then
            prtnr = nodel(prtlpos+1)
         else
            prtnr = prtlpos
         endif

         if (xmting) then
		istat = mxwrit(data, nitem*iteml, prtnr, itype)
         else if (rcving) then
		istat = mxread(data, nitem*iteml, prtnr, itype)
         endif
c
c        End of one phase.
c
         istep = istep/2
         go to 50
      endif
c
c     All done -- exit "successful"
c
60    mxcombv1 = 0
      return
      end
