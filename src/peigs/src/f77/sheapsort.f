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
      subroutine sheapsort(n, ra)
c     
c     
c     
      integer n
      real ra(n)
c     
      integer l, ir, i,j
      real rra
c     
c     
      if ( n .eq. 0 ) return
      if ( n .eq. 1 ) return
      
c     
      l = n/2 + 1
      ir = n
c     
 10   continue
      if ( l .gt. 1 ) then
         l = l-1
         rra = ra(l)
      else
         rra = ra(ir)
         ra(ir) = ra(1)
         ir = ir - 1
         if ( ir .eq. 1 ) then
            ra(1) = rra
            return
         endif
      endif
c     
      if ( ir .eq. 1 ) goto 30
c     
      i = l
      j = l + l
 20   continue
      if ( j .le. ir ) then
         if ( j .lt. ir ) then
            if ( ra(j) .lt. ra(j+1)) j = j+1
         endif
         if ( rra .lt. ra(j)) then
            ra(i) = ra(j)
            i = j
            j = j + j
         else
            j = ir + 1
         endif
         go to 20
      endif
      ra(i) = rra
      goto 10
c
 30   return
      end
      
