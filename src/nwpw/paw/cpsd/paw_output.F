*
* $Id: paw_output.F,v 1.3 2003-03-11 18:09:35 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *         psi_marat_write		*
*     *             			*
*     ***********************************

      subroutine paw_energy_output(un,nion,nel,E)
      implicit none 
      integer un
      integer nion
      integer nel
      double precision E(*)

*     *** total energy ***
      write(un,1430) E(1),E(1)/nion

c     write(un,1453)

*     *** kinetic energy ***
      write(un,1431) E(2)+E(3),(E(2)+E(3))/nel
*     *** plane wave  ***
      write(un,1432) E(2),E(2)/nel
*     *** local basis ***
      write(un,1433) E(3),E(3)/nel
      
*     *** coulomb energy ***
      write(un,1434)
*     *** plane wave  ***
      write(un,1435) E(4),E(4)/nel
*     *** local basis ***
      write(un,1436) E(5),E(5)/nel
*     *** multipole ***
      write(un,1437) E(6),E(6)/nel

*     *** exch-corr energy ***
      write(un,1438)
*     *** plane wave  ***
      write(un,1439) E(7),E(7)/nel
*     *** local basis ***
      write(un,1440) E(8),E(8)/nel

*     *** local pseudopot energy ***
      write(un,1441)
*     *** plane wave  ***
      write(un,1442) E(9),E(9)/nel
*     *** local basis ***
      write(un,1443) E(10),E(10)/nel

 1430 FORMAT(//' Total Energy:              ',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/'  kinetic energy:           ',E19.10,' (',E15.5,'/ion)')
 1432 FORMAT('      plane wave:           ',E19.10,' (',E15.5,'/ion)')
 1433 FORMAT('     local basis:           ',E19.10,' (',E15.5,'/ion)')
 1434 FORMAT(/'  coulomb energy:')
 1435 FORMAT('      plane wave:           ',E19.10,' (',E15.5,'/ion)')
 1436 FORMAT('     local basis:           ',E19.10,' (',E15.5,'/ion)')
 1437 FORMAT('       multipole:           ',E19.10,' (',E15.5,'/ion)')
 1438 FORMAT(/' exch-corr energy:')
 1439 FORMAT('       plane wave:          ',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT('      local basis:          ',E19.10,' (',E15.5,'/ion)')
 1441 FORMAT(/' pseudopot energy:')
 1442 FORMAT('       plane wave:          ',E19.10,' (',E15.5,'/ion)')
 1443 FORMAT('      local basis:          ',E19.10,' (',E15.5,'/ion)')
 1453 FORMAT('--------------------------------------------------------')
      end


