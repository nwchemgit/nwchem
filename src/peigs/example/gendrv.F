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
      program gendriv
c
c  dummy driver for the c calling the general eigensystem solver
c     for the mxsubs.cpp calls
c
      integer ienv(4)
c     
c     
c     
      external mxinit, mxpara
      external mxpend
c
#ifdef TCGMSG
      call pbeginf
#endif
      call mxinit
c     
c     call mxpara(ienv)
c     
c     call setdbg(1)
c     
c
      call main1()
c     
      call mxpend ()
      stop
      end
