c
c     $Id$                       
c

*     ***********************************
*     *			 		*
*     *		ofpw_cgsd_energy	*
*     *					*
*     ***********************************

      real*8 function ofpw_cgsd_energy(newpsi)
      implicit none
      logical newpsi

#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"

      integer MASTER
      parameter (MASTER=0)

      logical stalled,value,oprint,psitmp
      integer taskid
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
      logical     control_print,control_psi_tmp
      character*4 ion_aname
      integer     control_it_in, control_it_out,control_minimizer
      integer     ion_nion,ion_katm,psi_ne,psi_ispin,control_rtdb
      integer     psi_iptr_hml
      real*8      control_tole,control_tolc
      real*8      ion_rion,psi_eigenvalue,psi_hml,ion_amass
      external psp_semicore,pspw_charge_found,pspw_qmmm_found
      external control_out_of_time,pspw_SIC,control_DOS,pspw_HFX
      external ion_aname
      external control_it_in, control_it_out,control_minimizer
      external ion_nion,ion_katm,psi_ne,psi_ispin,control_rtdb
      external psi_iptr_hml
      external control_print,control_psi_tmp
      external control_tole,control_tolc
      external ion_rion,psi_eigenvalue,psi_hml,ion_amass
      character*4 psi_ab_irrep_name
      external    psi_ab_irrep_name
      integer     control_symmetry
      external    control_symmetry


      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))

      psitmp = control_psi_tmp()


*     **** set the minimizer ****
      call dcopy(30,0.0d0,0,E,1)
      minimizer = control_minimizer()

*     **** generate phaze factors and local psp and core density ****
      call phafac()
      !call electron_gen_vl_potential()
      if (psp_semicore(0)) call semicore_density_update()


*     :::::::::::  begin electron iteration  :::::::::::::::::::::::
      if (oprint) then
         !write(luout,1300)
         write(luout,1301)
         !write(luout,1302)
         write(luout,1304)
         if (minimizer.eq.1) write(luout,1305)
         if (minimizer.eq.2) write(luout,1306)
         if (minimizer.eq.3) write(luout,1307)
         if (minimizer.eq.4) write(luout,1308)
         if (minimizer.eq.5) write(luout,1309)
         if (minimizer.eq.6) write(luout,1311)
         if (minimizer.eq.7) write(luout,1312)
         call util_flush(luout)
      end if

      stalled = .false.
      deltae  = -1.0d-03
      icount=0
      bfgscount=0
      it_in  = control_it_in()
      it_out = control_it_out()
      tole   = control_tole()
      tolc   = control_tolc()
      E(1)=0.0d0
      if (oprint) then
        if (minimizer.eq.5) then
          call nwpw_message(7)
        else
          call nwpw_message(2)
        end if
      end if

      if (oprint) call util_flush(luout)
      !if (minimizer.gt.1) call pspw_Grsm_list_start()
      if (minimizer.eq.5) it_out = 1
      if (newpsi) then
        !call sdminimize(15)
      end if
   2  continue
         icount = icount + 1
         if (stalled) then
           !call sdminimize(0)
           bfgscount = 0 
         end if

         deltae_old = deltae
         if (minimizer.eq.1) then
           bfgscount = bfgscount + 1
           call ofpw_cgminimize(E,deltae,deltac,bfgscount,it_in)
         else if (minimizer.eq.2) then
           bfgscount = bfgscount + 1
           !call ofpw_bfgsminimize(E,deltae,deltac,bfgscount,minimizer)
           write(*,*) "should not be here either"
         else 
           bfgscount = bfgscount + 1
           write(*,*) "should not be here"
         end if

         if ((dabs(deltae).gt.dabs(deltae_old)).or.
     >       (dabs(deltae).gt.1.0d-2)          .or.
     >       (deltae.gt.0.0d0)) then
            stalled = .true.
         else
            stalled = .false.
         end if

         if ((oprint).and.(minimizer.ne.5)) then 
           write(luout,1310) icount*it_in,E(1),deltae,deltac
           call util_flush(luout)
         end if
         if (deltae.gt.0.0d0) then
            if (oprint) 
     >       write(luout,*) 
     >       ' *** energy going up. iteration not terminated'
             stalled = .true.
            !go to 3
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.tole).and.
     >       (deltac.lt.tolc)) then
            if (oprint) 
     >       write(luout,*) ' *** tolerance ok. iteration terminated'
            go to 3
         end if
         if (control_out_of_time()) then
            if (oprint) 
     >       write(luout,*) ' *** out of time. iteration terminated'
           go to 3
         end if

      if (icount.lt.it_out) go to 2
      if (oprint) 
     >write(luout,*) '*** arrived at the Maximum iteration.  terminated'

*     :::::::::::  end of electron iteration loop  :::::::::::::::::::::

   3  continue
      !if (minimizer.gt.1) call pspw_Grsm_list_end()
      if (oprint) CALL NWPW_MESSAGE(3)


*     **** diagonalize hamiltonian and rotate psi  ****
*     **** but skip if doing an SIC calculations   ****
      !call psi_spin_density(en)
      en(1) = 0.0d0
      en(2) = 0.0d0
      

*     **** geometrical center and center of mass of the cluster ****
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      !if (control_symmetry().gt.0) call psi_ab_gen_irrep_names(.false.)
      if (oprint) then
         !neall = (psi_ne(1)+psi_ne(2))
         neall = 0
         write(luout,1304)
         write(luout,1410)

         write(luout,*)
         !write(luout,1320) en(1),en(psi_ispin()),' (real space)'
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
         write(luout,1470) E(5),E(5)/ion_nion()
         write(luout,1480) E(6),E(6)/neall
         write(luout,1490) E(7),E(7)/neall
         write(luout,1495) E(8),E(8)/neall
         write(luout,1496) E(9),E(9)/neall
         write(luout,1497) E(10),E(10)/neall
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial

        if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
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

*        **** write out KS eigenvalues ****

      end if
      !if (control_symmetry().gt.0) call psi_ab_kill_irrep_names()

      rtdb = control_rtdb()
      !if (control_DOS()) call psi_DOS(rtdb)

c*     **** debug - energies written to rtdb for numerical stress ****
c      rtdb = control_rtdb()
c      value =           btdb_put(rtdb,'pspw:E_ke',      mt_dbl,1,E(6))
c      value = value.and.btdb_put(rtdb,'pspw:E_hartree', mt_dbl,1,E(3))
c      value = value.and.btdb_put(rtdb,'pspw:E_xc',      mt_dbl,1,E(4))
c      value = value.and.btdb_put(rtdb,'pspw:E_ewald',   mt_dbl,1,E(5))
c      value = value.and.btdb_put(rtdb,'pspw:E_local',   mt_dbl,1,E(7))
c      value = value.and.btdb_put(rtdb,'pspw:E_nonlocal',mt_dbl,1,E(8))
c      if (.not. value) call errquit(
c     >   'cgsd_noit_energy: numerical stress - error writing rtdb',0)

      ofpw_cgsd_energy = E(1)
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
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV)  occ=',F5.3))
 1512 FORMAT(2(E18.7,' (',F8.3,'eV)',A4))
 1513 FORMAT(2(E18.7,' (',F8.3,'eV)',A4,' occ=',F5.3))

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

