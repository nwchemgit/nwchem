*
* $Id: nwpw_message.F,v 1.2 2002-02-27 20:21:15 bylaska Exp $
*

      SUBROUTINE nwpw_MESSAGE(N)
      implicit none
#include "stdio.fh"
*     ================================================
*     THIS ROUTINE WRITES MESSAGE, TIME, AND DATE ON
*     I/O UNIT=10.  THIS ROUTINE IS MACHINE-DEPENDENT.
*     ================================================

      CHARACTER*17 MSG(6)
      integer N
      character*26 dd
      DATA MSG(1) /'JOB STARTED'/
      DATA MSG(2) /'ITERATION STARTED'/
      DATA MSG(3) /'ITERATION ENDED'/
      DATA MSG(4) /'JOB COMPLETED'/
      DATA MSG(5) /'JOB TERMINATED'/
      DATA MSG(6) /'ITERATION STARTED'/
 
      IF(N.GT.6.OR.N.LE.0) RETURN
      call util_date(dd)
      WRITE(luout,1000) MSG(N), dd
 1000 FORMAT(5X,'>>>  ',A,' AT ',A24,'  <<<')
      IF(N.EQ.2)  WRITE(luout,1001)
 1001 FORMAT(4x,'iter.',11x,'Energy ',7x,' DeltaE',7x,'DeltaRho '/
     >     4x,'------------------------------------------------------')
      IF(N.EQ.6)  WRITE(luout,1002)
 1002 FORMAT(4x,'iter.',9x,'KE+Energy',13x,'Energy',8x,'KE_psi',
     >     8x,'KE_ion',/
     >     4x, 
     > '-------------------------------------------------------------',
     > '---------')
      RETURN
      END
