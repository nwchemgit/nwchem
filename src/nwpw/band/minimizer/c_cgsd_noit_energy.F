c     $Id$                       

*     ***************************
*     *				*
*     *	   c_cgsd_noit_energy	*
*     *				*
*     ***************************

      real*8 function c_cgsd_noit_energy()
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"
#include "errquit.fh"

      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,mprint
      real*8  E(50),virial,EV
      real*8  f0,f1,f2,f3,f4,f5,f6
      integer rtdb,neall
      integer i,NN,nb

*     **** external functions ****
       integer  control_rtdb
       external control_rtdb
      logical     band_HFX
      logical     cpsp_semicore,control_out_of_time,cpsi_spin_orbit
      logical     control_fractional,control_print,brillioun_print
      character*4 ion_aname
      integer     control_it_in,control_it_out,control_minimizer
      integer     ion_nion,ion_katm,cpsi_ne,cpsi_ispin
      real*8      control_tole,control_tolc,cpsi_occupation_brdcst
      real*8      ion_rion,cpsi_eigenvalue_brdcst,cpsi_svector_brdcst
      real*8      ion_amass
      external    band_HFX
      external    cpsp_semicore,control_out_of_time,cpsi_spin_orbit
      external    control_fractional,control_print,brillioun_print
      external ion_aname
      external control_it_in,control_it_out,control_minimizer
      external ion_nion,ion_katm,cpsi_ne,cpsi_ispin
      external control_tole,control_tolc,cpsi_occupation_brdcst
      external ion_rion,cpsi_eigenvalue_brdcst,cpsi_svector_brdcst
      external ion_amass
      integer  brillioun_nbrillioun
      real*8   brillioun_k_brdcst,brillioun_ks_brdcst
      real*8   brillioun_weight_brdcst
      external brillioun_nbrillioun
      external brillioun_k_brdcst,brillioun_ks_brdcst
      external brillioun_weight_brdcst
      real*8   cpsi_smearfermi,cpsi_smearcorrection
      external cpsi_smearfermi,cpsi_smearcorrection
      logical  ion_disp_on
      external ion_disp_on

      real*8   cpsi_1ke
      external cpsi_1ke
      real*8   ttss


      call Parallel_taskid(taskid)

      mprint = (taskid.eq.MASTER).and.control_print(print_medium)


*     **** generate phaze factors and local psp and core density ****
c      call cphafac()
c      call cphafac_k()
      call c_electron_gen_vl_potential()
      if (cpsp_semicore(0)) call c_semicore_density_update()


      call c_nominimize(E)

*     **** diagonalize hamiltonian ****
      call cpsi_1gen_hml()
      call cpsi_diagonalize_hml()
      if (control_fractional()) then
         call cpsi_1define_occupation(0.0d0,.false.)
         E(12) = cpsi_smearcorrection()
         E(1)  = E(1) + E(12)
      end if


*     **** generate svectors ****
      if (cpsi_spin_orbit()) then
        call cpsi_1rotate2()
        call cpsi_2to1()
        call cpsi_1gen_svector()
      end if


**:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (mprint) then
         neall = (cpsi_ne(1)+cpsi_ne(2))
         write(luout,1304)
         write(luout,1410)

         write(luout,*)
c         write(luout,1320) en(1),en(cpsi_ispin()),' (real space)'
         write(luout,1430) E(1),E(1)/ion_nion()
         write(luout,1440) E(2),E(2)/neall
         write(luout,1450) E(3),E(3)/neall
         write(luout,1460) E(4),E(4)/neall
         if (band_HFX()) then
           write(luout,1457) E(26),E(26)/neall
         end if
         write(luout,1470) E(5),E(5)/ion_nion()
         if (control_fractional()) write(luout,1471) E(12),E(12)/neall

         write(luout,1480) E(6),E(6)/neall
         write(luout,1490) E(7),E(7)/neall
         write(luout,1495) E(8),E(8)/neall
         write(luout,1496) E(9),E(9)/neall
         write(luout,1497) E(10),E(10)/neall
         write(luout,1499) E(11),E(11)/neall
         virial = (E(11)+E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial
         if (band_HFX())  then
           write(luout,1502) E(27),E(27)/neall
         end if
        if (ion_disp_on()) then
            write(luout,1720) E(33)
        end if
      end if

      NN=cpsi_ne(1)-cpsi_ne(2)
      EV=27.2116d0
      if (mprint) then
        if (control_fractional()) then
          if (cpsi_ispin().eq.1) then
            write(luout,1507) cpsi_smearfermi(1),cpsi_smearfermi(1)*EV
          else
            write(luout,1507) cpsi_smearfermi(1),cpsi_smearfermi(1)*EV,
     >                        cpsi_smearfermi(2),cpsi_smearfermi(2)*EV
          end if
        end if
      end if

      if (brillioun_print()) then
      do nb=1,brillioun_nbrillioun()
         f0 = brillioun_weight_brdcst(nb) 
         f1 = brillioun_ks_brdcst(1,nb) 
         f2 = brillioun_ks_brdcst(2,nb) 
         f3 = brillioun_ks_brdcst(3,nb) 
         f4 = brillioun_k_brdcst(1,nb) 
         f5 = brillioun_k_brdcst(2,nb) 
         f6 = brillioun_k_brdcst(3,nb) 
         if (mprint) then
            write(luout,1508) nb,f0,f1,f2,f3,f4,f5,f6
            write(luout,1500)
         end if
         if (cpsi_spin_orbit()) then
           if (mprint) write(luout,1511)
           do i=0,cpsi_ne(1)-1
             f1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i) 
             f2 = cpsi_svector_brdcst(nb,cpsi_ne(1)-i,3) 
             f3 = cpsi_svector_brdcst(nb,cpsi_ne(1)-i,2) 
             f4 = cpsi_svector_brdcst(nb,cpsi_ne(1)-i,1) 
             f5 = cpsi_occupation_brdcst(nb,cpsi_ne(1)-i,1)
             if (mprint)
     >          write(luout,1512) f1,f1*EV,f2,f3,f4,f5
           end do
         else
           do i=0,NN-1
             f1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i) 
             f2 = cpsi_occupation_brdcst(nb,1,cpsi_ne(1)-i)
             if (mprint) write(luout,1510) f1,f1*EV,f2
           end do
           do i=0,cpsi_ne(2)-1
             f1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i-NN) 
             f2 = cpsi_occupation_brdcst(nb,1,cpsi_ne(1)-i-NN) 
             f3 = cpsi_eigenvalue_brdcst(nb,2,cpsi_ne(2)-i) 
             f4 = cpsi_occupation_brdcst(nb,2,cpsi_ne(2)-i)
             if (mprint) 
     >          write(luout,1510) f1,f1*EV,f2,f3,f3*EV,f4
           end do
         end if
      end do

      else
         if (mprint) write(luout,1501)
      endif




*     **** diagonalize hamiltonian and rotate psi ****
c     call psi_spin_density(en)
c     call psi_1gen_hml()
c     call psi_diagonalize_hml()
c
c     call psi_1rotate2()
c     call psi_2to1()


*     **** debug - energies written to rtdb for numerical stress ****
      rtdb = control_rtdb()
      value =           btdb_put(rtdb,'band:E_ke',      mt_dbl,1,E(6))
      value = value.and.btdb_put(rtdb,'band:E_hartree', mt_dbl,1,E(3))
      value = value.and.btdb_put(rtdb,'band:E_xc',      mt_dbl,1,E(4))
      value = value.and.btdb_put(rtdb,'band:E_ewald',   mt_dbl,1,E(5))
      value = value.and.btdb_put(rtdb,'band:E_local',   mt_dbl,1,E(7))
      value = value.and.btdb_put(rtdb,'band:E_nonlocal',mt_dbl,1,E(8))
      value = value.and.btdb_put(rtdb,'band:E_disp',    mt_dbl,1,E(33))
      if (.not. value) call errquit(
     >   'c_cgsd_noit_energy: numerical stress - error writing rtdb',0,
     &       RTDB_ERR)


      c_cgsd_noit_energy = E(1)
      return
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'======================')
 1301 FORMAT(  '= energy calculation =')
 1302 FORMAT(  '======================')
 1304 FORMAT(/)
 1305 FORMAT(10X,'== bundled Grassmann conjugate gradient iteration ==')
 1306 FORMAT(10X,'======== bundled Grassmann lmbfgs iteration ========')
 1307 FORMAT(10X,'============ Grassmann diis iteration ==============')
 1308 FORMAT(10X,'=== bundled Stiefel conjugate gradient iteration ===')
 1309 FORMAT(10X,'=========== band by band diis iteration ============')
 1310 FORMAT(I8,E20.10,3E15.5)
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1330 FORMAT(/' comparison between hamiltonian and lambda matrix')
 1340 FORMAT(I3,2I3,' H=',E16.7,', L=',E16.7,', H-L=',E16.7)
 1350 FORMAT(/' orthonormality')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT(I3)
 1380 FORMAT(' ''',a,'''',I4)
 1390 FORMAT(I3)
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT( ' final position of ions:')
 1430 FORMAT(/' total     energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT( ' smearing  energy    :',E19.10,' (',E15.5,'/electron)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1491 FORMAT( ' K.S. Vl+Vqm/mm      :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1499 FORMAT( ' K.S. V_nlso energy  :',E19.10,' (',E15.5,'/electron)')
 1500 FORMAT(/' orbital energies:')
 1501 FORMAT(/' orbital energies not printed',
     >        ' - number of k-points is very large')
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1507 FORMAT(/' Fermi energy =',2(E18.7,' (',F8.3,'eV)'))
 1508 FORMAT(/' Brillouin zone point: ',i6,
     >       /'    weight=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1511 FORMAT(33x,"Spin(Sz,Sy,Sz)")
 1512 FORMAT(E18.7,' (',F8.3,' eV) (',F8.3,' zhat + ',F8.3,' yhat +',
     > F8.3,' xhat ) occ=',F5.3)

 1720 FORMAT(/' Dispersion energy   :',E19.10)

      end

