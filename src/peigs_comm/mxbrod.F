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
      integer function mxbrod (data,isrc,isize,nprocs,nodel,itype)
c     mxbcst -- broadcast (from isrc to all others in nodel).
c
c     Function:  This routine provides most of the functionality of
c                EXPRESS routine 'kxbrod', in a portable implementation.
c                The calling sequence and semantics are intended to be
c                identical to kxbrod, except for limitations described
c                below.
c
c     Arguments:
c         data		(typeless in/out vector) data to be broadcast.
c                          Length of the data vector is isize bytes.
c         isrc		(integer input) length, in bytes, of data to
c                          be broadcast.
c         nprocs      (integer input) number of processors (length of nodel
c                        vector), or EXPRESS constant 'ialnod' (from
c                        the /xpress/ common block)
c         nodel       (integer input vector) list of processors across
c                        which to combine.  ** This list must be IDENTICAL
c                        on all processors within the list.  Even the list
c                        ordering must be the same on all participating
c                        processors. **
c         itype       (integer input) message type to be used.
c
c     Returns:  For now, this version always returns "successful", using
c               the kxbrod convention for that value.  Future versions
c               should return failure codes for incorrect arguments and/or
c               communications errors, but it is unclear exactly when and
c               how such returns should be made.  Cavaet caller!
c
c     Limitations:
c        1. The nodel list must be IDENTICAL on all participating
c           processors.  Deadlock will result if it isnt.
c
c     Method:  Non-pipelined binary spanning tree on the nodel list 
c              as provided by the caller.
c
c     Revision History:
c        8/3/91, initial version by Rik Littlefield, loosely patterned
c                  after then-current mxcomb.
c
      integer*4 data(*)
      integer isrc
      integer isize
      integer nprocs
      integer nodel(nprocs)
      integer itype
c
c     Local variables -- see code for usage.
c
      integer srcpos
      integer dstsrc
      integer mypnum, mylpos
      integer prtnr, prtlpos
      logical rcving, xmting
      integer istep
      integer istat
      integer env(4)
      integer i
      integer nppart
c
c     EXPRESS-environment common block that contains certain variables
c     and constants.  In this routine, we need 'ialnod'.
c
      common /xpress/ nocare, norder, nonode, ihost, ialnod, ialprc
c
c     If nprocs = ialnod, then we need to know the actual number
c     of processors participating in the combine ('nppart').
c
      istat = mxpara(env)
      if (nprocs.ne.ialnod) then
         nppart = nprocs
      else
         nppart = env(2)
      endif
c
c     Quick exit for special case of no broadcast required.
c
      if (nppart.lt.2) go to 60
c
c     Determine our processor number and where we are in the list of
c     processor nodes.  If we're not in this list, then we're not
c     participating, so immediately exit.  The basic output from this
c     code is mylpos, which is this processors position in the list
c     in terms of first element = position 0, and srcpos, which is
c     the position of processor isrc.
c
      mypnum = env(1)
      if (nprocs.ne.ialnod) then
         mylpos = -1
         do 10 i = 1,nprocs
            if (mypnum .eq. nodel(i)) mylpos = i-1
            if (isrc .eq. nodel(i)) srcpos = i-1
10       continue
         if (mylpos.eq.-1) go to 60
      else
         mylpos = mypnum
         srcpos = isrc
      endif
c.      write (*,*) mxbrod: mypnum,mylpos,isrc,srcpos=,
c.     -            mypnum,mylpos,isrc,srcpos
c
c     OK, so this processor is participating.  We now need to do a series
c     of ceil(log2(nproc)) broadcast phases.  Who does what in each phase
c     is determined by each processors distance (in the list) from
c     the source processor, so well compute that now as dstsrc.
c
      dstsrc = mylpos - srcpos
      if (dstsrc.lt.0) dstsrc = dstsrc + nppart
c
      istep = 1
50    if (istep.lt.nppart) then
c
         xmting = .false.
         rcving = .false.

         if (dstsrc.lt.istep) then
            if (dstsrc+istep.lt.nppart) then
               xmting = .true.
               prtlpos = mylpos + istep
               if (prtlpos.ge.nppart) prtlpos = prtlpos - nppart
            endif
c
         else if (dstsrc.lt.2*istep) then
            rcving = .true.
            prtlpos = mylpos - istep
            if (prtlpos.lt.0) prtlpos = prtlpos + nppart
         endif

         if (rcving .or. xmting) then
            if (nprocs.ne.ialnod) then
               prtnr = nodel(prtlpos+1)
            else
               prtnr = prtlpos
            endif
         endif

c.      write (*,*) 'me,prtnr,rcving,xmting = ',mypnum,prtnr,rcving,xmting

         if (xmting) then
            istat = mxwrit(data, isize, prtnr, itype)
         else if (rcving) then
            istat = mxread(data, isize, prtnr, itype)
         endif

c         End of one phase

         istep = istep * 2
         go to 50
      endif
c
c     All done -- exit "successful"
c
60    mxbrod = 0
      return
      end
