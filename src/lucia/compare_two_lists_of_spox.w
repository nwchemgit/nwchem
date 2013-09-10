      PROGRAM COMPARE_TWO_LISTS_OF_SPOX
*
* Compare two lists of spin orbital excitations, read in from 
* input in a form comparable to standard output
*
* Jeppe Olsen, Sittingin a Hotel in Washington, wondering
* about the difference between two setups..
*
      INCLUDE 'implicit.inc'
      PARAMETER(NWORD=1000000)
*
      INTEGER IWORK(NWORD)
*
      READ(5,*) NGAS, NTP1, NTP2
      KCAAB1 = 1
      KCAAB2 = KCAAB + 4*NGAS*NTP1
      IFREE  = KCAAB2 + 4*NGAS*NTP2
*. Read in the lists is written in standard output
      CALL READ_LIST_OF_SPOX(5,NGAS,NTP1,IWORK(KCAAB1)
      CALL READ_LIST_OF_SPOX(5,NGAS,NTP2,IWORK(KCAAB2)
*. And compare the two lists
*. Those in 1 that is not in 2
      WRITE(6,*) ' Determing CAAB in 1 but not in 2'
      CALL COMPARE_TWO_LISTS_OF_SPOX_INNER
     &     (NGAS,IWORK(KCAAB1),NTP1,
     &           IWORK(KCAAB2),NTP2)
*. TO BE COMPLETED..
      SUBROUTINE 
      DO ICAAB = 1, NTPa
*. Skip an empty line and a line with text (Included...)
       READ(5,*)

      SUBROUTINE READ_LIST_OF_SPOX(IUNIT,NGAS,NSPOBEX_TP,ICAAB)
*
* Read a list of CAAB operators from some outputfile
*
*. Jeppe Olsen, Aug. 2009
*
      INCLUDE 'implicit.inc'
      INTEGER ICAAB(NGAS,4,NSPOBEX_TP)
      CHARACTER SOMETEXT(24)
*
      DO I = 1, NSPOBEX_TP
*. Skip two lines: one with empty and one with the text 
       READ(UNIT,*)
       READ(UNIT,*)
       DO ICOMP = 1, 4
         READ(UNIT,*) SOMETEXT, (ICAAB(IGAS,ICOMP,I),IGAS = 1, NGAS)
       END DO
      END DO
*
      RETURN
      END
