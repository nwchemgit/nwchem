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
*  -- PFFT routine (version 1.0) --
*     Pacific Northwest Laboratory
*     April 5, 1995
*
*======================================================================
      INTEGER FUNCTION MAPCHK( MAP, N )
*
      INTEGER           N, MAP(N)

*    Make sure all entries in MAP are 0 < MAP(1:N) < MXNPRC(),

*    Returns MAPCHK =  0 if ok,
*                     -1 if bad.

      INTEGER           NAPROC, I

      INTEGER           MXNPRC
      EXTERNAL          MXNPRC

*-----------------------------------------------------------------------

      NAPROC = MXNPRC()

      MAPCHK = 0
*
      IF( N .LE. 0 )
     $  RETURN
*
      DO 10 I = 1, N
        IF( MAP(I) .LT. 0  .OR.  MAP(I) .GE. NAPROC )
     $    MAPCHK = -1
   10 CONTINUE
*
      RETURN
      END
