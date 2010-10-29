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
      SUBROUTINE MAPFIL( ME, N, MAP, NME, MYMAP  )
*
      INTEGER           ME, N, NME
      INTEGER           MAP(N), MYMAP(*)
*    
*    Count the number of instance, NME of ME
*    in a MAP   with N elements.
*
*    Store the values of I such that MAP(I) = ME in
*          MYMAP(j), j = 1 to NME
*  
      INTEGER           I
*-----------------------------
*  
      NME = 0
*
      IF( N .LE. 0 )
     $  RETURN
*
      DO 10 I = 1, N
        IF( MAP(I) .EQ. ME ) THEN
          NME          = NME + 1
          MYMAP( NME ) = I
        ENDIF
   10 CONTINUE
*
      RETURN
      END
