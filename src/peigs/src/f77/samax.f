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
      real function samax(n,dx,incx)
c
c     returns the max of the absolute values.

      real dx(*),dtemp
      integer i,incx,n,nincx
c
      samax = 0.0
      dtemp = 0.0
      if( n.le.0 .or. incx.le.0 )return

      if(incx.ne.1) then
c
c        code for increment not equal to 1
c
        nincx = n*incx
        do 10 i = 1,nincx,incx
          dtemp = max( dtemp, abs(dx(i)) )
   10   continue
      else

c        code for increment equal to 1
c
        do 20 i = 1, n
          dtemp = max( dtemp, abs(dx(i)) )
   20   continue

      endif

      samax = dtemp

      return
      end
