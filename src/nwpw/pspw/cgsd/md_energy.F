c
c     $Id$                       
c

*     ***************************
*     *				*
*     *		md_energy	*
*     *				*
*     ***************************
      real*8 function md_energy(E)
      implicit none
      real*8 E(*)

#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      integer MASTER
      parameter (MASTER=0)

      logical value,oprint
      integer taskid
      integer i,j,ms
      integer icount
      real*8  EV,virial
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  e_lj,e_q,e_spring,eion
      integer rtdb

*     **** external functions ****
      logical  control_print
      integer  control_version
      integer  ion_nion,ion_katm,control_rtdb
      real*8   ion_rion,ion_amass,ewald_e,ion_ion_e
      external control_print
      external control_version
      external ion_nion,ion_katm,control_rtdb
      external ion_rion,ion_amass,ewald_e,ion_ion_e


*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_mmq_Q_E,pspw_qmmm_spring_E
      real*8   pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_mmq_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))

*     **** set the minimizer ****
      !call dcopy(30,0.0d0,0,E,1)

*     **** generate phaze factors ****
      !if (control_version().eq.3) call ewald_phafac()

*     :::::::::::  get energy  :::::::::::::::::::::::
      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()


      E(1)  = 0.0d0
      E(2)  = eion
      E(3)  = 0.0d0
      E(4)  = 0.0d0
      E(5)  = eion
      E(6)  = 0.0d0
      E(7)  = 0.0d0
      E(8)  = 0.0d0
      E(9)  = 0.0d0
      E(10) = 0.0d0
      if (pspw_qmmm_found()) then
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_mmq_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(2)  = E(2) + e_lj + e_q + e_spring
         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring
         E(14) = pspw_qmmm_LJ_Emix()
         E(14) = E(14) + pspw_qmmm_Q_Emix()
         !E(23) = E(23) + E(14)
         !E(24) = E(24) + E(14)*E(14)

      end if

      !E(25) = E(25) + E(1)
      !E(26) = E(26) + E(1)*E(1)



*     **** get pspw_charge energies ****
c      if (pspw_charge_found()) then
c         E(19)  = 0.0d0
c         E(20)  = pspw_charge_Energy_ion()
c         E(21)  = pspw_charge_Energy_charge()
c         E(1)   = E(1) + E(20) + E(21)
c      end if


c*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
c      if (oprint) then
c         write(luout,1304)
c         write(luout,1410)
c
c         write(luout,*)
c         write(luout,1430) E(1),E(1)/ion_nion()
c
c
c         if (pspw_charge_found()) then
c            write(luout,1431)
c            write(luout,1432)
c            write(luout,1433) (E(1)-E(19)-E(20)-E(21)),
c     >         (E(1)-E(19)-E(20)-E(21))/ion_nion()
c         end if
c
c         write(luout,1470) E(5),E(5)/ion_nion()
c         if (pspw_qmmm_found()) then
c            write(luout,1700)
c            write(luout,1701)
c            write(luout,1702) E(11)
c            write(luout,1703) E(12)
c            write(luout,1704) E(13)
c         end if
c         if (pspw_charge_found()) then
c            write(luout,1800)
c            write(luout,1801)
c            write(luout,1805) E(19)+E(20)+E(21)
c            write(luout,1802) E(19)
c            write(luout,1803) E(20)
c            write(luout,1804) E(21)
c         end if
c
c      end if
c      call ecce_print1 ('total energy', mt_dbl, E(1), 1)
c      call ecce_print1 ('nuclear repulsion energy', mt_dbl, E(5), 1)

      md_energy = E(2)
      return

 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'======================')
 1301 FORMAT(//'== Energy Calculation ==')
 1302 FORMAT(  '======================')
 1304 FORMAT(/)
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
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
 1491 FORMAT( ' Vl+Vqm/mm           :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' V_nl    (planewave) :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' V_Coul  (planewave) :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' V_xc.   (planewave) :',E19.10,' (',E15.5,'/electron)')
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


 1700 FORMAT(/' LJ/residualCoulomb/vib/CAV Energies')
 1701 FORMAT( ' -----------------------------------')
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



*     *******************************
*     *				    *
*     *	    md_energy_gradient      *
*     *				    *
*     *******************************
      subroutine md_energy_gradient(G1)
      implicit none
      real*8 G1(3,*)

#include "stdio.fh"
#include "util.fh"

      logical allow_translation,lprint,mprint
      integer MASTER
      parameter (MASTER=0)
      integer i,k,taskid,nion,nion1
      integer i1
      real*8  GG,fmax,fatom
      real*8  fmx,fmy,fmz
      real*8  fmx2,fmy2,fmz2

*     **** external functions ****
      logical     psp_semicore,pspw_charge_found,pspw_qmmm_found
      logical     control_allow_translation,ion_q_FixIon,control_print
      character*4 ion_aname,pspw_charge_aname
      integer     ion_katm,ion_nion,control_version
      integer     pspw_charge_nion
      real*8      ion_rion,pspw_charge_rion
      real*8      pspw_charge_charge
      external psp_semicore,pspw_charge_found,pspw_qmmm_found
      external control_allow_translation,ion_q_FixIon,control_print
      external ion_aname,pspw_charge_aname
      external ion_katm,ion_nion,control_version
      external pspw_charge_nion
      external ion_rion,pspw_charge_rion
      external pspw_charge_charge
      logical  pspw_bqext
      external pspw_bqext
      !integer  ion_nion_mm,ion_nion_qm
      !external ion_nion_mm,ion_nion_qm

      allow_translation = control_allow_translation()
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

*     **** generate phaze factors ****
      if (control_version().eq.3) call ewald_phafac()

      call dcopy(3*nion,0.0d0,0,G1,1)
      if (control_version().eq.3) call ewald_f(G1)
      if (control_version().eq.4) call ion_ion_f(G1)
      if (pspw_qmmm_found()) call pspw_qmmm_mmq_fion(G1)

c      if (pspw_charge_found()) then
c        if(pspw_bqext()) then
c           call pspw_charge_charge_Fion(G1)
c         else
c           nion1 = ion_nion()
c           call pspw_charge_Fion_Fcharge(G1,G1(1,nion1+1))
c           call pspw_charge_Fcharge(G1(1,nion1+1))
c         end if
c      end if

*     **** remove ion forces using ion_FixIon ****
      call ion_FixIon(G1)

      if (.not.allow_translation) then
        call center_F_mass(G1,fmx,fmy,fmz)
        do i=1,nion
         G1(1,i) = G1(1,i) - fmx
         G1(2,i) = G1(2,i) - fmy
         G1(3,i) = G1(3,i) - fmz
        end do
      end if
      call center_F_mass(G1,fmx2,fmy2,fmz2)

      GG = 0.0d0
      fmax = 0.0d0
      do i=1,nion
         GG = GG + G1(1,i)**2 + G1(2,i)**2 + G1(3,i)**2
         fatom = dsqrt(G1(1,i)**2+G1(2,i)**2 +G1(3,i)**2)
         if (fatom.gt.fmax)  fmax = fatom
      end do

      call Parallel_taskid(taskid)
      mprint = ((taskid.eq.MASTER).and.control_print(print_high))
      lprint = ((taskid.eq.MASTER).and.control_print(print_medium))

      if (taskid.eq.MASTER) then
        if (mprint) write(luout,1301)
        if (lprint) then
        write(luout,1304)
        if (.not.allow_translation) write(luout,1400) fmx,fmy,fmz
        write(luout,1304)
        write(luout,1410)
        end if

*      **** print out positions ***
        if (mprint) then
        write(luout,1420)
        do I=1,ion_nion()
          if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3)
          else
           write(6,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3)
          end if
        end do

c*       **** print out charge positions ***
c        if (pspw_charge_found().and.(.not.pspw_bqext())) then
c          do i=1,pspw_charge_nion()
c            i1 = ion_nion() + i
c            if (ion_q_FixIon(i1)) then
c            write(luout,1193) i1,pspw_charge_aname(i),
c     >                    (pspw_charge_rion(K,i),K=1,3),
c     >                    pspw_charge_charge(i)
c            else
c            write(luout,1192) i1,pspw_charge_aname(i),
c     >                    (pspw_charge_rion(K,i),K=1,3),
c     >                    pspw_charge_charge(i)
c            end if
c          end do
c        end if

        end if

c*       **** print out forces ***
        if (lprint) then
        write(luout,1421)
        write(luout,1190)(i,ion_aname(I),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())

c*       **** print out charge forces ***
c        if (pspw_charge_found().and.(.not.pspw_bqext())) then
c          do i=1,pspw_charge_nion()
c            i1 = ion_nion() + i
c            write(luout,1190) i1,pspw_charge_aname(i),
c     >                    (G1(K,i1),K=1,3)
c          end do
c        end if
c
c        write(luout,1210) fmx2,fmy2,fmz2  
c        write(luout,1425)
c        write(luout,1426) dsqrt(GG),
c     >                    dsqrt(GG)/dble(nion),
c     >                    fmax,fmax*(27.2116d0/0.529177d0)

        end if
      end if

c     call dscal(3*nion,(-1.0d0),G1,1)

      return
 1190 FORMAT(5X, I4, A5,  ' (',3F11.5,' )')
 1191 FORMAT(5X, I4, A5,  ' (',3F11.5,' ) - fixed')
 1192 FORMAT(5X, I4, A5,  ' (',3F11.5,' ) q=',F8.3)
 1193 FORMAT(5X, I4, A5,  ' (',3F11.5,' ) q=',F8.3,' - fixed')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'========================')
 1301 FORMAT(//'== Gradient Calculation ==')
 1302 FORMAT(  '========================')
 1304 FORMAT(/)
 1400 FORMAT('Translation force removed: (',3F11.5,')')
 1410 FORMAT(10X,'=============  Ion Gradients =================')
 1425 FORMAT(10X,'===============================================')
 1426 FORMAT(10X,'|F|       =',E15.6,
     >      /10x,'|F|/nion  =',E15.6,
     >      /10x,'max|Fatom|=',E15.6,1x,'(',F8.3,'eV/Angstrom)'//)
 1420 FORMAT( ' Ion Positions:')
 1421 FORMAT( ' Ion Forces:')
      end

