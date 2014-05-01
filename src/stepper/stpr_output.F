C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_output(STEP,COORD,BCKSTP, mxgrad,
     &                       convggm, convgg, convge )
c     $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "util.fh"
#include "mafdecls.fh"
      LOGICAL BCKSTP, oprint_info
      LOGICAL SECO, SWSEC, SWFIR
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / SWITCH / SECO, SWSEC, SWFIR
      double precision mxgrad, stepmax, stepnorm
      double precision convggm, convgg, convge
      DIMENSION STEP(NAT3),
     &     COORD(3,NATOM)
      double precision grad_max, grad_norm 
      external grad_max, grad_norm
c
      oprint_info = util_print('information',print_default)
c
      CALL stpr_getp(GMXMUM,GLNGTH,SLNGTH,EOLD,ENERGY,
     &     EDIFF,E2OLD,ERATIO,TRUSTO,TRUST,TRUSTN,
     &     E2NEW,NWALK,NNEG)
c
      IF ((ICALC.EQ.2).AND.(SECO.AND.(IWCTR.EQ.1))) THEN
         WRITE(6,*)'  Stepping switched to 1-st derivatives'
      ENDIF
      IF ((ICALC.EQ.2).AND.SWFIR) THEN
         WRITE(6,*)'  Using Analytical 1st derivatives. '
      ENDIF
      IF ((ICALC.EQ.1).AND.SWSEC) THEN
         WRITE(6,*)'  Using Analytical 1st and 2nd derivatives. '
      ENDIF

C     ECCE output
      stepmax  = grad_max  ( step, NATOM )
      stepnorm = grad_norm ( step, NATOM )
      if ( IWCTR .eq. 1 ) then
C       initial geometry
        call ecce_print2 ('geometry', mt_dbl, coord, 3, 3, NATOM)
      endif
      call ecce_print1 ('step max',  mt_dbl, stepmax,  1 )
      call ecce_print1 ('step norm', mt_dbl, stepnorm, 1 )
      call ecce_print1 ('delta energy', mt_dbl, EDIFF, 1 )

      if (oprint_info) then
*... does not reproduce correct output for reset and walk         WRITE(6,1000)NWALK
         WRITE(6,1000)iwctr
         WRITE(6,1006)EOLD
         WRITE(6,1007)ENERGY
         WRITE(6,1008)EDIFF,CONVGE 
         WRITE(6,1009)E2OLD
         WRITE(6,1010)ERATIO
         IF(BCKSTP)THEN
            WRITE(6,1012)TRUST
         ELSE
            WRITE(6,1011)TRUSTO
         ENDIF
         WRITE(6,1001)NNEG
         write(6,1032)mxgrad,convggm
         WRITE(6,1002)GLNGTH,convgg
      endif
      GLNGTH = GLNGTH/SQRT(DBLE(NATOM))
      if (oprint_info) then
         WRITE(6,1022)GLNGTH
         WRITE(6,1005)SLNGTH
         WRITE(6,1014)E2NEW
      endif
      if (util_print('old coords',print_high)) then
C     
C     Write old coordinates.
C     
         WRITE(6,1015)
         WRITE(6,1016)
         DO 5 I = 1,NATOM
            WRITE(6,1017)I,(COORD(J,I),J=1,3)
 5       CONTINUE
      endif
      if (util_print('step', print_high)) then
C     
C     Write step.
C     
         WRITE(6,1018)
         WRITE(6,1019)
         DO 6 I = 1,NATOM
            JB = (I-1)*3 + 1
            JE = I*3
            WRITE(6,1017)I,(STEP(J),J=JB,JE)
 6       CONTINUE
      endif
      RETURN
 1000 FORMAT('         Walk statistics for step: ',I3)
 1001 FORMAT('   Number of negative eigenvalues: ',I3)
 1002 FORMAT('                    Gradient norm: ',D20.12,
     &       ' threshold: ',d20.6)
 1022 FORMAT('                RMS Gradient norm: ',D20.12)
 1032 FORMAT('   Max absolute value of Gradient: ',D20.12,
     &       ' threshold: ',d20.6)
 1003 FORMAT('                  Value of lambda: ',D20.12)
 1004 FORMAT('          Value of scale factor A: ',D20.12)
 1005 FORMAT('  Largest component of step v-tor: ',D20.12)
 1006 FORMAT('     Previous actual total energy: ',D20.12)
 1007 FORMAT('      Current actual total energy: ',D20.12)
 1008 FORMAT('         Actual energy difference: ',D20.12,
     &       ' threshold: ',d20.6)
 1009 FORMAT('    Second order energy predicted: ',D20.12)
 1010 FORMAT('             Ratio of differences: ',D20.12)
 1011 FORMAT('               Fixed trust radius: ',D20.12)
 1012 FORMAT('   Modified trust radius for step: ',D20.12)
 1014 FORMAT('  New 2nd order energy prediction: ',D20.12)
 1015 FORMAT('  ',18('-'),' Input coordinates ',18('-'))
 1016 FORMAT('  ATOM:',12x,'X',14x,'Y',14x,'Z')
 1017 FORMAT('  ',I3,6X,3(1X,F14.10))
 1018 FORMAT('  ',20('-'),' Current step ',21('-'))
 1019 FORMAT('  STEP:',12x,'X',14x,'Y',14x,'Z')
      END
C> @}
