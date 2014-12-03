c     $Id$                       


*     ***************************
*     *                         *
*     *    paw_cgsd_noit_energy *
*     *                         *
*     ***************************

      real*8 function paw_cgsd_noit_energy()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ***
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,oprint
      real*8  E(30),en(2),en1(2),en2(2),EV
      integer rtdb,i,NN,neall

*     **** external functions ****
      logical  control_print
      integer  ion_nion,paw_psi_ne,paw_psi_ispin,control_rtdb
      real*8   paw_psi_eigenvalue
      external control_print
      external ion_nion,paw_psi_ne,paw_psi_ispin,control_rtdb
      external paw_psi_eigenvalue

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))
   
*     **** generate phaze factors and local psp and core density ****
      call paw_electron_gen_vl_potential()


      call paw_nominimize(E)


*     **** diagonalize hamiltonian and rotate psi  ****
*     **** but skip if doing an SIC calculations   ****
c     BUG!!!!!!
      call paw_psi_spin_density(en1,en2)
      en(1) = en1(1)+en2(1)
      en(2) = en1(2)+en2(2)

      call paw_psi_1gen_hml()
      call paw_psi_diagonalize_hml()


*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (oprint) then
         neall = (paw_psi_ne(1)+paw_psi_ne(2))
         write(luout,1304)
         write(luout,1410)
         write(luout,*)
c         write(luout,1321) r_charge,' (real space)'
         write(luout,1320) en(1), en(paw_psi_ispin()),' (real space)'
         write(luout,1322) en1(1),en1(paw_psi_ispin()),' (real space)'
         write(luout,1323) en2(1),en2(paw_psi_ispin()),' (real space)'
         call paw_energy_output(luout,ion_nion(),neall,E)

*        **** write out KS eigenvalues ****
         write(luout,1500)
         NN=paw_psi_ne(1)-paw_psi_ne(2)
         EV=27.2116d0
         do I=1,NN
           write(luout,1510) paw_psi_eigenvalue(1,I),
     >                       paw_psi_eigenvalue(1,I)*EV
         end do
         do I=1,paw_psi_ne(2)
           write(luout,1510)  paw_psi_eigenvalue(1,I+NN),
     >                        paw_psi_eigenvalue(1,I+NN)*EV,
     >                        paw_psi_eigenvalue(2,I),
     >                        paw_psi_eigenvalue(2,I)*EV
         end do

      end if



*     **** debug - energies written to rtdb for numerical stress ****
      rtdb = control_rtdb()
      value =           btdb_put(rtdb,'paw:E_ke',      mt_dbl,1,
     >                           E(2)+E(3))
      value = value.and.btdb_put(rtdb,'paw:E_hartree', mt_dbl,1,
     >                           E(4)+E(5))
      value = value.and.btdb_put(rtdb,'paw:E_mult',    mt_dbl,1,
     >                           E(6))
      value = value.and.btdb_put(rtdb,'paw:E_xc',      mt_dbl,1,
     >                           E(7)+E(8))
      value = value.and.btdb_put(rtdb,'paw:E_local',   mt_dbl,1,
     >                           E(9)+E(10))
      if (.not. value) call errquit(
     >   'paw_cgsd_noit_energy: numerical stress - error writing rtdb',
     >    0,RTDB_ERR)


      paw_cgsd_noit_energy = E(1)
      return
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'======================')
 1301 FORMAT(//'== Energy Calculation ==')
 1302 FORMAT(  '======================')
 1304 FORMAT(/) 
 1305 FORMAT(10X,'====== Grassmann conjugate gradient iteration ======')
 1306 FORMAT(10X,'============ Grassmann lmbfgs iteration ============')
 1307 FORMAT(10X,'============ Grassmann diis iteration ==============')
 1308 FORMAT(10X,'======= Stiefel conjugate gradient iteration =======')
 1309 FORMAT(10X,'============= Kohn-Sham scf iteration ==============')
 1310 FORMAT(I8,E20.10,3E15.5)
 1311 FORMAT(10X,'======= Projected Grassmann lmbfgs iteration =======')
 1312 FORMAT(10X,'============= Stiefel lmbfgs iteration =============')
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1322 FORMAT('     plane-wave part:         ',F11.5,'       ',F11.5,A)
 1323 FORMAT('      augmented part:         ',F11.5,'       ',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Diagonal Elements of Hamiltonian matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1341 FORMAT(I3,2I3,' H=',E16.6)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
c1410 FORMAT(10X,'=============  summary of results  =================')
 1410 FORMAT('==  Summary Of Results  ==')
 1420 FORMAT( ' final position of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1431 FORMAT(/' QM Energies')
 1432 FORMAT( '------------')
 1433 FORMAT( ' total  QM energy    :',E19.10,' (',E15.5,'/ion)')
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
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 
 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' QM/MM energy                     :',E19.10)
 1703 FORMAT(/' MM/MM energy                     :',E19.10)
 1704 FORMAT(/' MM Polarization energy           :',E19.10)
 1705 FORMAT( ' MM Vibration energy              :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy            :',E19.10)
 1707 FORMAT( ' - MM Charge Field/QM Electron    :',E19.10)
 1708 FORMAT( ' - MM Charge Field/QM Ion         :',E19.10)
 1709 FORMAT( ' - MM LJ/QM LJ                    :',E19.10)
 1710 FORMAT( ' - MM Charge Field/MM Charge Field:',E19.10)
 1711 FORMAT( ' - MM LJ/MM LJ                    :',E19.10)
 
 1800 FORMAT(/' Charge Field Energies')
 1801 FORMAT( ' ---------------------')
 1802 FORMAT( ' - Charge Field/Electron    :',E19.10)
 1803 FORMAT( ' - Charge Field/Ion         :',E19.10)
 1804 FORMAT( ' - Charge Field/Charge Field:',E19.10)
 1805 FORMAT( ' Charge Field Energy        :',E19.10)

      end

