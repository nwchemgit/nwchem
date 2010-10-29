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
      program machine
c
c
      double precision dlamch
      real slamch
c
      write(*,*)
      write(*,*) ' Double Precision results'
      write(*,1000) ' depsilon ', dlamch('e')
      write(*,1000) ' dbase ', dlamch('b')
      write(*,1000) ' dsafeulp ', dlamch('s')
      write(*,1000) ' dlamch(u)  ', dlamch('u')
c
      write(*,*)
      write(*,*) ' Single Precision results'
      write(*,1000) ' depsilon ', slamch('e')
      write(*,1000) ' dbase ', slamch('b')
      write(*,1000) ' dsafeulp ', slamch('s')
      write(*,1000) ' slamch(u)  ', slamch('u')
      write(*,*)
 1000 FORMAT( A12, 1X, 1P, E26.16 )
      stop
      end



