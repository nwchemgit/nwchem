c     $Id: cgsd_noit_energy.F,v 1.11 2008-11-17 17:25:44 bylaska Exp $                       
*     ***************************
*     *				*
*     *	   cgsd_noit_energy	*
*     *				*
*     ***************************

      real*8 function cgsd_noit_energy()
      implicit none

#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"

      integer MASTER,taskid
      parameter (MASTER=0)

      logical stalled,value,oprint,psitmp
      integer minimizer
      integer i,j,ms,neall,NN
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(30)
      integer rtdb


*     **** external functions ****
      logical     psp_semicore,pspw_charge_found,pspw_qmmm_found
      logical     control_out_of_time,pspw_SIC,control_DOS,pspw_HFX
      logical     control_print,control_psi_tmp,control_fractional
      character*4 ion_aname
      integer     control_it_in, control_it_out,control_minimizer
      integer     ion_nion,ion_katm,psi_ne,psi_ispin,control_rtdb
      integer     psi_iptr_hml
      real*8      control_tole,control_tolc,psi_occupation
      real*8      ion_rion,psi_eigenvalue,psi_hml,ion_amass
      external psp_semicore,pspw_charge_found,pspw_qmmm_found
      external control_out_of_time,pspw_SIC,control_DOS,pspw_HFX
      external ion_aname
      external control_it_in, control_it_out,control_minimizer
      external ion_nion,ion_katm,psi_ne,psi_ispin,control_rtdb
      external psi_iptr_hml
      external control_print,control_psi_tmp,control_fractional
      external control_tole,control_tolc,psi_occupation
      external ion_rion,psi_eigenvalue,psi_hml,ion_amass
      character*4 psi_ab_irrep_name
      external    psi_ab_irrep_name
      integer     control_symmetry,control_version
      external    control_symmetry,control_version
      logical  psp_pawexist
      external psp_pawexist
      real*8   psi_smearfermi,psi_smearcorrection
      external psi_smearfermi,psi_smearcorrection



      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))


*     **** generate phaze factors and local psp and core density ****
      call dcopy(30,0.0d0,0,E,1)
      call phafac()
      if (control_version().eq.3) call ewald_phafac()
      call electron_gen_vl_potential()
      if (psp_semicore(0)) call semicore_density_update()


      call nominimize(E)

*     **** diagonalize hamiltonian ****
      call psi_spin_density(en)
      call psi_1gen_hml()
      call psi_diagonalize_hml()
      if (control_fractional()) then
         call psi_1define_occupation(.true.)
        E(28) = psi_smearcorrection()
        E(1)  = E(1) + E(28)
      end if

      if (pspw_SIC()) then
        call psi_1gen_hml_g()
      end if

*     **** geometrical center and center of mass of the cluster ****
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (control_symmetry().gt.0) call psi_ab_gen_irrep_names(.false.)
      if (oprint) then
         neall = (psi_ne(1)+psi_ne(2))
         write(luout,1304)
         write(luout,1410)

         write(luout,*)
         write(luout,1320) en(1),en(psi_ispin()),' (real space)'
         write(luout,1430) E(1),E(1)/ion_nion()

         if (pspw_qmmm_found()) then
            write(luout,1431)
            write(luout,1432)
            write(luout,1433)  E(1)-E(11),
     >                        (E(1)-E(11))/ion_nion()
         end if

         if (pspw_charge_found()) then
            write(luout,1431)
            write(luout,1432)
            write(luout,1433) (E(1)-E(19)-E(20)-E(21)),
     >         (E(1)-E(19)-E(20)-E(21))/ion_nion()
         end if

         write(luout,1440) E(2),E(2)/neall
         write(luout,1450) E(3),E(3)/neall
         write(luout,1460) E(4),E(4)/neall
         if (pspw_SIC()) then
           write(luout,1455) E(22),E(22)/neall
           write(luout,1456) E(23),E(23)/neall
         end if
         if (pspw_HFX()) then
           write(luout,1457) E(26),E(26)/neall
         end if
         write(luout,1470) E(5),E(5)/ion_nion()
         if (control_fractional()) write(luout,1471) E(28),E(28)/neall

         write(luout,1480) E(6),E(6)/neall
         write(luout,1490) E(7),E(7)/neall
         write(luout,1495) E(8),E(8)/neall
         write(luout,1496) E(9),E(9)/neall
         write(luout,1497) E(10),E(10)/neall
         if (pspw_SIC())  then
           write(luout,1499) E(24),E(24)/neall
           write(luout,1501) E(25),E(25)/neall
         end if
         if (pspw_HFX())  then
           write(luout,1502) E(27),E(27)/neall
         end if

         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial

c        **** paw energies ****
         if (psp_pawexist()) then
         write(luout,1680) E(14),E(14)/ion_nion()
         write(luout,1681) E(15),E(15)/ion_nion()
         end if

        if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
c
c            write(luout,1703) E(14)+E(15)
c            write(luout,1710) E(14)
c            write (luout,1711) E(15)
c
c            write(luout,1704) E(16)
c            write(luout,1705) E(17)
c            write(luout,1706) E(18)
        end if
        if (pspw_charge_found()) then
            write(luout,1800)
            write(luout,1801)
            write(luout,1805) E(19)+E(20)+E(21)
            write(luout,1802) E(19)
            write(luout,1803) E(20)
            write(luout,1804) E(21)
        end if

*        **** write out <psi|H|psi> matrix ****
         if (pspw_SIC()) then
           do ms=1,psi_ispin()
             if (ms.eq.1) write(luout,1331)
             if (ms.eq.2) write(luout,1332)
             !*** call Gainsville matrix output ***
             call output(dbl_mb(psi_iptr_hml(ms,1,1)),
     >                    1,psi_ne(ms),1,psi_ne(ms),
     >                    psi_ne(ms),psi_ne(ms),1)

c            do j=1,psi_ne(ms)
c            do i=1,psi_ne(ms)
c              write(luout,1341) ms,i,j,psi_hml(ms,i,j)
c            end do
c            end do

           end do

         end if


*        **** write out KS eigenvalues ****
           NN=psi_ne(1)-psi_ne(2)
           EV=27.2116d0
          if (control_fractional()) then
             if (psi_ispin().eq.1) then
               write(luout,1507) psi_smearfermi(1),psi_smearfermi(1)*EV
             else
             write(luout,1507) psi_smearfermi(1),psi_smearfermi(1)*EV,
     >                         psi_smearfermi(2),psi_smearfermi(2)*EV
             end if
           end if

           write(luout,1500)
           if (control_symmetry().gt.0) then
             if (control_fractional()) then
             do I=1,NN
               write(luout,1513) psi_eigenvalue(1,I),
     >                           psi_eigenvalue(1,I)*EV,
     >                           psi_ab_irrep_name(I),
     >                           psi_occupation(1,I)
             end do
             do I=1,psi_ne(2)
               write(luout,1513)  psi_eigenvalue(1,I+NN),
     >                            psi_eigenvalue(1,I+NN)*EV,
     >                            psi_ab_irrep_name(I),
     >                            psi_occupation(1,I+NN),
     >                            psi_eigenvalue(2,I),
     >                            psi_eigenvalue(2,I)*EV,
     >                            psi_ab_irrep_name(I+psi_ne(1)),
     >                            psi_occupation(2,I)
             end do
             else
             do I=1,NN
                 write(luout,1512) psi_eigenvalue(1,I),
     >                           psi_eigenvalue(1,I)*EV,
     >                           psi_ab_irrep_name(I)
             end do
             do I=1,psi_ne(2)
               write(luout,1512)  psi_eigenvalue(1,I+NN),
     >                            psi_eigenvalue(1,I+NN)*EV,
     >                            psi_ab_irrep_name(I),
     >                            psi_eigenvalue(2,I),
     >                            psi_eigenvalue(2,I)*EV,
     >                            psi_ab_irrep_name(I+psi_ne(1))
             end do
             end if


           else
             if (control_fractional()) then
             do I=1,NN
               write(luout,1511) psi_eigenvalue(1,I),
     >                           psi_eigenvalue(1,I)*EV,
     >                           psi_occupation(1,I)
             end do
             do I=1,psi_ne(2)
               write(luout,1511)  psi_eigenvalue(1,I+NN),
     >                            psi_eigenvalue(1,I+NN)*EV,
     >                            psi_occupation(1,I+NN),
     >                            psi_eigenvalue(2,I),
     >                            psi_eigenvalue(2,I)*EV,
     >                            psi_occupation(2,I)
             end do
             else
             do I=1,NN
               write(luout,1510) psi_eigenvalue(1,I),
     >                           psi_eigenvalue(1,I)*EV
             end do
             do I=1,psi_ne(2)
               write(luout,1510)  psi_eigenvalue(1,I+NN),
     >                            psi_eigenvalue(1,I+NN)*EV,
     >                            psi_eigenvalue(2,I),
     >                            psi_eigenvalue(2,I)*EV
             end do
             end if


           end if

      end if
      if (control_symmetry().gt.0) call psi_ab_kill_irrep_names()





*     **** debug - energies written to rtdb for numerical stress ****
      rtdb = control_rtdb()
      value =           rtdb_put(rtdb,'pspw:E_ke',      mt_dbl,1,E(6))
      value = value.and.rtdb_put(rtdb,'pspw:E_hartree', mt_dbl,1,E(3))
      value = value.and.rtdb_put(rtdb,'pspw:E_xc',      mt_dbl,1,E(4))
      value = value.and.rtdb_put(rtdb,'pspw:E_ewald',   mt_dbl,1,E(5))
      value = value.and.rtdb_put(rtdb,'pspw:E_local',   mt_dbl,1,E(7))
      value = value.and.rtdb_put(rtdb,'pspw:E_nonlocal',mt_dbl,1,E(8))
      value = value.and.rtdb_put(rtdb,'pspw:E_hsic',    mt_dbl,1,E(22))
      value = value.and.rtdb_put(rtdb,'pspw:E_xcsic',   mt_dbl,1,E(23))
      value = value.and.rtdb_put(rtdb,'pspw:E_hfx',     mt_dbl,1,E(26))
      if (.not. value) call errquit(
     >   'cgsd_noit_energy: numerical stress - error writing rtdb',0,
     >       RTDB_ERR)

*     **** debug - energies written to rtdb for numerical forces ****
      value = value.and.rtdb_put(rtdb,'pspw:E_LJ',   mt_dbl,1,E(11))
      value = value.and.rtdb_put(rtdb,'pspw:E_Q',    mt_dbl,1,E(12))
      if (.not. value) call errquit(
     >   'cgsd_noit_energy: numerical force - error writing rtdb',0,
     >       RTDB_ERR)


      cgsd_noit_energy = E(1)
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
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1331 FORMAT(/' Elements of Hamiltonian matrix (up/restricted)')
 1332 FORMAT(/' Elements of Hamiltonian matrix (down)')
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
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT( ' smearing  energy    :',E19.10,' (',E15.5,'/electron)')
 1480 FORMAT(/' kinetic (planewave) :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' V_local (planewave) :',E19.10,' (',E15.5,'/electron)')
 1491 FORMAT( ' K.S. Vl+Vqm/mm      :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' V_nl    (planewave) :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' V_coul  (planewave) :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' V_xc    (planewave) :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1499 FORMAT( ' K.S. SIC-hartree energy  :',E19.10,
     >        ' (',E15.5,'/electron)')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1507 FORMAT(/' Fermi energy =',2(E18.7,' (',F8.3,'eV)'))
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV)  occ=',F5.3))
 1512 FORMAT(2(E18.7,' (',F8.3,'eV)',A4))
 1513 FORMAT(2(E18.7,' (',F8.3,'eV)',A4,' occ=',F5.3))

 1680 FORMAT(/' kinetic (loc. basis):',E19.10,' (',E15.5,'/ion)')
 1681 FORMAT( ' kinetic       (core):',E19.10,' (',E15.5,'/ion)')
 1690 FORMAT( ' V_local (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1696 FORMAT( ' coulomb (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1697 FORMAT( ' exc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')



 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy                        :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy          :',E19.10)
 1704 FORMAT( ' MM Vibrational energy            :',E19.10)
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

