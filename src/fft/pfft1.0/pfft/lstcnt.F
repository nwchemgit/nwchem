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
      INTEGER FUNCTION LSTCNT( IFIND, LIST, N )
*
*  DESCRIPTION
*  ===========
*
*  Routine to determine how many times the INTEGER IFIND appears
*  in the INTEGER array LIST(1:N).
*
*  ARGUMENTS
*  =========

*   IFIND....... (Input) INTEGER
*                An integer.
*
*   LIST ....... (Input) INTEGER, DIMENSION ( N )
*                A list of integers
*                
*   N .......... (Input) INTEGER
*                Number of entries in LIST.
*                N > 0.
*
*   LSTCNT ..... (Output) INTEGER
*                The number of times that IFIND appears in LIST(1:N).
*                LSTCNT >= 0.
*
* ----------------------------------------------------------------------
*
* ARGUMENTS
* =========
      INTEGER            IFIND, N
      INTEGER            LIST( N )
*
* LOCAL VARIABLES
* ===============
*
      INTEGER            KOUNT, I
*
*
* EXECUTABLE CODE
* ======================================================================
*
      KOUNT = 0
      DO 10 I = 1, N
        IF( LIST(I) .EQ. IFIND )
     $    KOUNT = KOUNT + 1
   10 CONTINUE

      LSTCNT = KOUNT

      RETURN
      END
