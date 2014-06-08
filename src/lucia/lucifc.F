*----------------------------------------------------------------------*
      SUBROUTINE WRTLUCIFC(ENERGY,CWFLABEL,IRELAX,
     &                     LUCMO,LSVCMO,
     &                     LU_FGEN,LSVFGEN,
     &                     LU_GAM1,LSVGAM1,
     &                     LU_RSPC,LSVRSPC,
     &                     LU_LAM2,FLAM2,LSVLAM2,
     &                     NSMOB,NSPIN,NTOOBS,NAOS,
     &                     NFRZ,NACT,NDEL)
*----------------------------------------------------------------------*
*     write LUCIA interface file (SIRIUS style)
*
*     Andreas, started 11-11-04 (Hellau!)
*----------------------------------------------------------------------*
      INCLUDE 'wrkspc.inc'

      INTEGER, PARAMETER ::
     &     NTEST = 00
      CHARACTER, PARAMETER ::  ! 'magic' label
     &     LABMAG(3)*8 = (/'********','********','********'/)

      LOGICAL, INTENT(IN) ::
     &     LSVCMO, LSVFGEN, LSVGAM1, LSVRSPC, LSVLAM2
      INTEGER, INTENT(IN) ::
     &     IRELAX, NSMOB, NSPIN,
     &     NTOOBS(NSMOB), NAOS(NSMOB), NFRZ(NSMOB),
     &     NACT(NSMOB),   NDEL(NSMOB),
     &     LUCMO, LU_FGEN, LU_GAM1, LU_RSPC, LU_LAM2
      REAL(8), INTENT(IN) ::
     &     ENERGY
      CHARACTER, INTENT(IN) ::
     &     CWFLABEL*80, FLAM2*8

      INTEGER ::
     &     LUCIFC, KLBUF, NBUF, NBUF_F, NBUF_G, NBUFCMO,
     &     II, IDUM, NAO, NTOOB

      INTEGER, EXTERNAL ::
     &     IOPEN_NUS
*----------------------------------------------------------------------*

      NAO = 0
      NTOOB = 0
      DO II = 1, NSMOB
        NAO = NAO + NAOS(II)
        NTOOB = NTOOB + NTOOBS(II)
      END DO

*----------------------------------------------------------------------*
*     Get buffer
*----------------------------------------------------------------------*
      IF (LSVCMO.OR.LSVFGEN.OR.LSVGAM1) THEN
        NBUFCMO = 0
        NBUF_F = 0
        NBUF_G = 0
        DO II = 1, NSMOB
          NBUFCMO = NBUFCMO + NAOS(II)*NTOOBS(II)
          NBUF_F  = NBUF_F  + NTOOBS(II)*NTOOBS(II)
        END DO
        NBUF_G = NBUF_F
        NBUF = 0
        IF (LSVCMO) NBUF = NBUFCMO
        IF (LSVFGEN) NBUF = MAX(NBUF,NBUF_F)
        IF (LSVGAM1) NBUF = MAX(NBUF,NBUF_G)
        IDUM = 0
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LUCIFC')
        CALL MEMMAN(KLBUF,NBUF,'ADDL  ',2,'BUFFER')
      END IF

*----------------------------------------------------------------------*
*     open file
*----------------------------------------------------------------------*
      LUCIFC = IOPEN_NUS('LUCIFC')

*----------------------------------------------------------------------*
*     Energy and general info
*----------------------------------------------------------------------*
      WRITE(LUCIFC) LABMAG,'GENERAL '
      WRITE(LUCIFC) CWFLABEL
      WRITE(LUCIFC) ENERGY
      WRITE(LUCIFC) IRELAX
      WRITE(LUCIFC) NTOOB,NAO,NSMOB,NSPIN,NAOS(1:NSMOB),NTOOBS(1:NSMOB),
     &              NFRZ(1:NSMOB),NACT(1:NSMOB),NDEL(1:NSMOB)

*----------------------------------------------------------------------*
*     Info on CMO-coefficients
*----------------------------------------------------------------------*
      IF(LSVCMO) THEN
        WRITE(LUCIFC) LABMAG,'TRANCMO '
        WRITE(LUCIFC) NTOOB,NAO,NSMOB,NSPIN,
     &       NAOS(1:NSMOB),NTOOBS(1:NSMOB)
        CALL VEC_FROM_DISC(WORK(KLBUF),NBUFCMO,1,-1,LUCMO)
        WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUFCMO-1)
        IF (NSPIN.EQ.2) THEN
          CALL VEC_FROM_DISC(WORK(KLBUF),NBUFCMO,0,-1,LUCMO)
          WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUFCMO-1)
        END IF
      END IF
        
*----------------------------------------------------------------------*
*     Info on general Fock-matrix
*----------------------------------------------------------------------*
      IF(LSVFGEN) THEN
        WRITE(LUCIFC) LABMAG,'INT FGEN'
        WRITE(LUCIFC) IRELAX,NTOOB,NSMOB,NSPIN,NTOOBS(1:NSMOB),
     &       NFRZ(1:NSMOB),NACT(1:NSMOB),NDEL(1:NSMOB)
        CALL VEC_FROM_DISC(WORK(KLBUF),NBUF_F,1,-1,LU_FGEN)
        WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUF_F-1)
        IF (NSPIN.EQ.2) THEN
          CALL VEC_FROM_DISC(WORK(KLBUF),NBUF_F,0,-1,LU_FGEN)
          WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUF_F-1)
        END IF
      END IF

*----------------------------------------------------------------------*
*     Info on one-particle density
*----------------------------------------------------------------------*
      IF (LSVGAM1) THEN
        WRITE(LUCIFC) LABMAG,'INT GAM1'
        WRITE(LUCIFC) IRELAX,NTOOB,NSMOB,NSPIN,NTOOBS(1:NSMOB),
     &       NFRZ(1:NSMOB),NACT(1:NSMOB),NDEL(1:NSMOB)
        CALL VEC_FROM_DISC(WORK(KLBUF),NBUF_G,1,-1,LU_GAM1)
        WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUF_G-1)
        IF (NSPIN.EQ.2) THEN
          CALL VEC_FROM_DISC(WORK(KLBUF),NBUF_G,0,-1,LU_GAM1)
          WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUF_G-1)
        END IF
      END IF

*----------------------------------------------------------------------*
*     Info on orbital response contribution to one-particle density
*----------------------------------------------------------------------*
      IF (LSVRSPC) THEN
        WRITE(LUCIFC) LABMAG,'INT KAP1'
        WRITE(LUCIFC) IRELAX,NTOOB,NSMOB,NSPIN,NTOOBS(1:NSMOB),
     &       NFRZ(1:NSMOB),NACT(1:NSMOB),NDEL(1:NSMOB)
        CALL VEC_FROM_DISC(WORK(KLBUF),NBUF_G,1,-1,LU_RSPC)
        WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUF_G-1)
        IF (NSPIN.EQ.2) THEN
          CALL VEC_FROM_DISC(WORK(KLBUF),NBUF_G,0,-1,LU_RSPC)
          WRITE(LUCIFC) WORK(KLBUF:KLBUF+NBUF_G-1)
        END IF
      END IF

*----------------------------------------------------------------------*
*     Info on irreducible two-particle density
*----------------------------------------------------------------------*
      IF (LSVLAM2) THEN
        WRITE(LUCIFC) LABMAG,'INT LAM2'
        WRITE(LUCIFC) IRELAX,NTOOB,NSMOB,NSPIN,NTOOBS(1:NSMOB),
     &       NFRZ(1:NSMOB),NACT(1:NSMOB),NDEL(1:NSMOB)
        WRITE(LUCIFC) FLAM2(1:8)
      END IF

      CALL RELUNIT(LUCIFC,'keep')

      IF (LSVCMO.OR.LSVFGEN.OR.LSVGAM1) THEN
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'LUCIFC')
      END IF

      END
c $Id$
