      SUBROUTINE print_output(istep,nprnt,ekin,ecoul,eshrt,ebond,eshel)

      implicit none

      include 'p_const.inc'
      include 'p_input.inc'
      include 'cm_temp.inc'

      integer  istep,nprnt

      real*8  ekin,ecoul,eshrt,ebond,eshel
      real*8  stptmp,stpke,stppe,stpte

      logical lnew

      data lnew/.true./

      save lnew

      stptmp=2.0*ekin/boltzmann/degfree
      stpke=ekin/convfct2
      stppe=ecoul+eshrt+ebond+eshel
      stpte=stppe+stpke
 
      if(lnew)write(output,'(/,a12,8a12)')'    STEP    ',
     $  '  TOTAL E.  ',' KINETIC E. ','POTENTIAL E.','  COUL. E.  ',
     $  '   VDW E.   ','  BOND E.   ','  SHELL E.  ',' TEMPERATURE'
      if(lnew)lnew=.false.

      if(mod(istep,nprnt).eq.0)write(output,'(i12,8f12.4)')
     $  istep,stpte,stpke,stppe,ecoul,eshrt,ebond,eshel,stptmp

      return

      END
c $Id$
