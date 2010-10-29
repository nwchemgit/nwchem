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
c PeIGS utility routine
c
c
      integer function menode(Npprocs, nproclist)
      implicit none
      integer Npprocs, nproclist(*)
c
c     Purpose:
c
c     returns the index of the first occurrance of
c     the current processor ( relative node number ) id
c     in a list of Npprocs processors
c     given in a list (nproclist)
c     
      integer mxmynd
      external mxmynd
      external xerbla
      integer indx, me
c     
c     
      me = mxmynd ()
      indx = 0
 10   continue
      indx = indx + 1
      if ( indx .gt. Npprocs ) goto 999
      if ( me .ne. nproclist(indx)) goto 10
      menode = indx - 1
      return
c     
  999 CONTINUE
* 999  write(*,*) 'ERROR: function menode, node: ',  
*     1     me, ' not in node list', Npprocs
c
c      call xerbla('menode', 1)
c
      menode = -1
      return
      end
      
