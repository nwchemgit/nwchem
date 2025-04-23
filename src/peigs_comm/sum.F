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
C	These are special purpose integer functions used in an mxcomb call
C	The routines must be integer functions, they must return 0, and
C	the result of a hopefully commutative and associative operation goes 
C	into the 1st parameter.  The 1st 2 input parameters can be any data 
C	type, and the # of bytes these parameters occupy is denoted in the 
C	3rd input parameter.  Usually the 3rd parameter is not used.
********************************************************************
      integer function sumi(a,b,size)
      integer a,b,size
      a = a + b
      sumi = 0
      return
      end
********************************************************************
      integer function sumd(a,b,size)
      double precision a,b
      integer size
      a = a + b
      sumd = 0
      return
      end
********************************************************************
      integer function sumdv(a,b,size,nelm)
c
      double precision a(nelm),b(nelm)
      integer size
      do 10 i = 1,nelm
         a(i) = a(i) + b(i)
 10   continue
      sumdv = 0
      return
      end
c     
      integer function sumiv(a,b,size,nelm)
      integer nelm
      integer a(nelm),b(nelm)
      integer size
c
      do 10 i = 1,nelm
         a(i) = a(i) + b(i)
 10   continue
      sumiv = 0
      return
      end
c      
C     Do this silliness to keep the transputers happy
********************************************************************
C	This function is defined to keep mxcomb consistent with mxcomb
C	by using a double precision complex function
c
      integer function sumdc(a,b,size)
      complex*16 a, b
      integer size
      a = a + b
      sumdc = 0
      return
      end
c

