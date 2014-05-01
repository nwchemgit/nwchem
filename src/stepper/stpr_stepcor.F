C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_stepcor(A,INEED,G,C,ETOT,NATD,
     &    CONVGE,CONVGG,CONVGGM,
     &    cnvgrd, rtdb, step_number)
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "util.fh"
#include "cstprfiles.fh"
      integer rtdb
      integer step_number
      logical dbug
      LOGICAL TRACK,RPATH,BCKSTP,CNVGRD
      LOGICAL SECO, SWSEC, SWFIR
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &                  ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / OLDES / ELAST,E2PRED
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / THRESH / RCONVG,RSTD
      COMMON / SWITCH / SECO, SWSEC, SWFIR
      COMMON / MASS / TOTM, NUMAS
      DIMENSION G(NATD*3)
      DIMENSION C(3,NATD)
      DIMENSION A(INEED)
      DIMENSION IPT(25)
C
      dbug=.false.
C
C     Open input/output files.
C
C     Tape 41 is used to store step parameters PSTEPP.
C
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP41)
C
C     Tape 44 is used to store previous tracked mode in order to
C     find current hessian eigenvector with maximum overlap.
C
      OPEN (UNIT=ITAP44, FORM='UNFORMATTED', FILE=FILE44,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP44)
C
C     Tape 45 is a temp file used to store HESS, GRAD, COORD
C     and ENERGY during step determination.
C
      OPEN (UNIT=ITAP45, FORM='UNFORMATTED', FILE=FILE45,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP45)
#if defined(OPENCLOSE)
      CLOSE (UNIT=ITAP41, STATUS='KEEP')
      CLOSE (UNIT=ITAP44, STATUS='KEEP')
      CLOSE (UNIT=ITAP45, STATUS='KEEP')
#endif
C
C     Calculate memory offsets.
C
      IPT( 1) = 1
      IPT( 2) = IPT( 1) + NAT3SQ
      IPT( 3) = IPT( 2) + NAT3
      IPT( 4) = IPT( 3) + NAT3
      IPT( 5) = IPT( 4) + NAT3SQ
      IPT( 6) = IPT( 5) + NAT3
      IPT( 7) = IPT( 6) + NAT3
      IPT( 8) = IPT( 7) + NAT3
      IPT( 9) = IPT( 8) + NAT3
      IPT(10) = IPT( 9) + NATOM
      IPT(11) = IPT(10) + NAT3
      IPT(12) = IPT(11) + NAT3TR
      IPT(13) = IPT(12) + NAT3SQ
      IPT(14) = IPT(13) + NAT3SQ
      IPT(15) = IPT(14) + NATOM
      IPT(16) = IPT(15) + NUMAS
      IPT(17) = IPT(16) + NAT3TR
      IPT(18) = IPT(17) + NAT3TR
      IPT(19) = IPT(18) + 3
      IPT(20) = IPT(19) + 9
C
C     Set up atomic masses.
C
      CALL stpr_setmass(A(IPT(16)))
C
C     Get user input: TRUST,NMODE,TRACK,
C     SELECT, RPATH, STD, etc, and equivalence data coming
C     in from NWCHEM to data locations required by STEPPER.
C
      if(dbug) then
         grad_max = dabsmax(natd*3,g)
         write(6,*)' grad max in stepcor ',grad_max
      endif
      CALL stpr_face
     &       ( A(IPT( 1)), A(IPT( 2)), A(IPT( 5)), A(IPT( 9)),
     &       A(IPT(11)), A(IPT(15)), A(IPT(16)),
     &       G, C, ETOT, TRUST, ENERGY, CONVGE, CONVGG, CONVGGM,
     &       NATD, NMODE, TRACK, BCKSTP, RPATH , rtdb, step_number)
      if(dbug) then
         grad_max = dabsmax(natd*3,g)
         write(6,*)' grad max in stepcor ',grad_max
      endif
C
C     If analytical hessian is unavailable form update hessian.
C
      IF(ICALC.EQ.1.AND..NOT.BCKSTP)THEN
        CALL stpr_FORMUP (
     &         A(IPT(1)),  A(IPT(2)),  A(IPT(8)),  A(IPT(5)),
     &         A(IPT(3)), A(IPT(13)) , rtdb)
      ENDIF
C
C     Obtain step.
C
      CALL stpr_driver (
     &       A(IPT( 1)), A(IPT( 2)), A(IPT( 3)), A(IPT( 4)),
     &       A(IPT( 5)), A(IPT( 6)), A(IPT( 7)), A(IPT( 8)),
     &       A(IPT( 9)), A(IPT(10)), A(IPT(11)), A(IPT(12)),
     &       A(IPT(13)), A(IPT(14)), A(IPT(16)), A(IPT(17)),
     &       A(IPT(18)), A(IPT(19)), ENERGY, E2NEW, TRUST,
     &       CONVGE, CONVGG, CONVGGM,
     &       NMODE, TRACK, RPATH, BCKSTP, CNVGRD, rtdb )
C
C     Switch between first and second derivatives.
C
c      IF (SECO.OR.SWFIR.OR.SWSEC) THEN
c        CALL FIRSEC
c      ENDIF
C
C     Put back coordinates into SUPERMOLECULE array:
C
      CALL stpr_putbck  ( A(IPT(5)), C )
C
C     Write out new interatom distances and angles.
C
c      CALL DISTAN
c      CALL ANGLE
C
C
#if defined(OPENCLOSE)
C     Open input/output files.
C
C     Tape 41 is used to store step parameters PSTEPP.
C
      OPEN (UNIT=ITAP41, FORM='UNFORMATTED', FILE=FILE41,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP41)
C
C     Tape 44 is used to store previous tracked mode in order to
C     find current hessian eigenvector with maximum overlap.
C
      OPEN (UNIT=ITAP44, FORM='UNFORMATTED', FILE=FILE44,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP44)
C
C     Tape 45 is a temp file used to store HESS, GRAD, COORD
C     and ENERGY during step determination.
C
      OPEN (UNIT=ITAP45, FORM='UNFORMATTED', FILE=FILE45,
     &      ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
      REWIND (UNIT=ITAP45)
#endif
C
C     Close open files and return.
C
      CLOSE (UNIT=ITAP41, STATUS='KEEP')
      CLOSE (UNIT=ITAP44, STATUS='KEEP') 
      CLOSE (UNIT=ITAP45, STATUS='DELETE')
      ELAST = ENERGY
      E2PRED = E2NEW
      RETURN
      END
C> @}
