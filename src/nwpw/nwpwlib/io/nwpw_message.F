*
* $Id: nwpw_message.F,v 1.10 2008-11-03 06:06:10 bylaska Exp $
*

      SUBROUTINE nwpw_MESSAGE(N)
      implicit none
#include "stdio.fh"
*     ================================================
*     THIS ROUTINE WRITES MESSAGE, TIME, AND DATE ON
*     I/O UNIT=10.  THIS ROUTINE IS MACHINE-DEPENDENT.
*     ================================================

      CHARACTER*17 MSG(10)
      integer N
      character*26 dd
      DATA MSG(1) /'JOB STARTED'/
      DATA MSG(2) /'ITERATION STARTED'/
      DATA MSG(3) /'ITERATION ENDED'/
      DATA MSG(4) /'JOB COMPLETED'/
      DATA MSG(5) /'JOB TERMINATED'/
      DATA MSG(6) /'ITERATION STARTED'/
      DATA MSG(7) /'ITERATION STARTED'/
      DATA MSG(8) /'ITERATION STARTED'/
      DATA MSG(9) /'ITERATION STARTED'/
      DATA MSG(10) /'ITERATION STARTED'/
 
      IF(N.GT.10.OR.N.LE.0) RETURN
      call util_date(dd)
      WRITE(luout,1000) MSG(N), dd
      IF(N.EQ.2)  WRITE(luout,1001)
      IF(N.EQ.6)  WRITE(luout,1002)
      IF(N.EQ.7)  WRITE(luout,1003)
      IF(N.EQ.8)  WRITE(luout,1004)
      IF(N.EQ.9)  WRITE(luout,1005)
      IF(N.EQ.10) WRITE(luout,1006)
      RETURN

 1000 FORMAT(5X,'>>>  ',A,' AT ',A24,'  <<<')
 1001 FORMAT(4x,'iter.',11x,'Energy ',7x,' DeltaE',7x,'DeltaRho '/
     >     4x,'------------------------------------------------------')
 1002 FORMAT(4x,'iter.',9x,'KE+Energy',13x,'Energy',8x,'KE_psi',
     >     8x,'KE_ion', 3x,'Temperature'/
     >     4x, 
     > '-------------------------------------------------------------',
     > '-----------------------')
 1003 FORMAT(4x,'iter.',11x,'Energy ',7x,' DeltaE',7x,'DeltaRho ',
     >                                             7x,'DeltaV  '/
     >     4x,'------------------------------------------------------',
     > '---------------')
 1004 FORMAT(4x,'iter.',11x,'Energy ',7x,' DeltaE',7x,'DeltaVhxc',
     >                                             7x,'alpha   '/
     >     4x,'------------------------------------------------------',
     > '---------------')
 1005 FORMAT(4x,'iter.',9x,'KE+Energy',13x,'Energy',8x,'KE_psi',
     >     8x,'KE_ion', 3x,'Temperature',3x,'Pressure'/
     >     4x, 
     > '-------------------------------------------------------------',
     > '----------------------------------')
 1006 FORMAT(4x,'iter.',9x,'KE+Energy',13x,'Energy',
     >     8x,'KE_ion', 3x,'Temperature'/
     >     4x, 
     > '-------------------------------------------------------------',
     > '-----------------------')


      END
