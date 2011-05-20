*
* $Id$
*

*     ***********************************
*     *                                 *
*     *         psi_energy_write        *
*     *                                 *
*     ***********************************

      subroutine paw_energy_output(un,nion,nel,E)
      implicit none 
      integer un
      integer nion
      integer nel
      double precision E(*),etmp
      
*     *** total energy ***
      write(un,1430) E(1),E(1)/nion
      write(un,1440) E(11),E(11)/nel

      etmp = E(4)+E(5)+E(6)
      write(un,1450) etmp,etmp/nel

      etmp = E(7)+E(8)
      write(un,1460) etmp,etmp/nel


      write(un,1530) E(2),E(2)/nel
      write(un,1531) E(4),E(4)/nel
      write(un,1532) E(7),E(7)/nel
      write(un,1533) E(9),E(9)/nel
  
      write(un,1540) E(3),E(3)/nel
      write(un,1541) E(5),E(5)/nel
      write(un,1542) E(8),E(8)/nel
      write(un,1543) E(10),E(10)/nel

      write(un,1550) E(6),E(6)/nion

c*     *** kinetic energy ***
c      write(un,1431) E(2)+E(3),(E(2)+E(3))/nel
c*     *** plane wave  ***
c      write(un,1432) E(2),E(2)/nel
c*     *** local basis ***
c      write(un,1433) E(3),E(3)/nel
c      
c*     *** coulomb energy ***
c      write(un,1434)
c*     *** plane wave  ***
c      write(un,1435) E(4),E(4)/nel
c*     *** local basis ***
c      write(un,1436) E(5),E(5)/nel
c*     *** multipole ***
c      write(un,1437) E(6),E(6)/nel
c
c*     *** exch-corr energy ***
c      write(un,1438)
c*     *** plane wave  ***
c      write(un,1439) E(7),E(7)/nel
c*     *** local basis ***
c      write(un,1440) E(8),E(8)/nel
c
c*     *** local pseudopot energy ***
c      write(un,1441)
c*     *** plane wave  ***
c      write(un,1442) E(9),E(9)/nel
c*     *** local basis ***
c      write(un,1443) E(10),E(10)/nel

 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')

 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1491 FORMAT( ' K.S. Vl+Vqm/mm      :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')

 1530 FORMAT(/' kinetic (planewave) :',E19.10,' (',E15.5,'/electron)')
 1531 FORMAT( ' coulomb (planewave) :',E19.10,' (',E15.5,'/electron)')
 1532 FORMAT( ' exc-cor (planewave) :',E19.10,' (',E15.5,'/electron)')
 1533 FORMAT( ' pseudo  (planewave) :',E19.10,' (',E15.5,'/electron)')

 1540 FORMAT(/' kinetic (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1541 FORMAT( ' coulomb (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1542 FORMAT( ' exc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1543 FORMAT( ' pseudo (loc. basis) :',E19.10,' (',E15.5,'/electron)')

 1550 FORMAT(/' coulomb (multipole) :',E19.10,' (',E15.5,'/ion)')

      end


