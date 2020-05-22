c     $Id$                       


*     ***************************
*     *                         *
*     *    paw_cgsd_energy      *
*     *                         *
*     ***************************

      real*8 function paw_cgsd_energy(newpsi)
      implicit none

#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"
#include "frac_occ.fh"
c#include "global.fh"
c#include "btdb.fh"
#include "errquit.fh"

      logical newpsi
      integer MASTER
      parameter (MASTER=0)

      logical stalled,value,oprint,psitmp
      integer taskid
      integer minimizer
      integer i,j,ms,neall,NN,ispin
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2),en1(2),en2(2),r_charge
      real*8  E(30),eorb,ecoul,exc
      integer rtdb

*     **** external functions ****
      logical     control_out_of_time
      logical     control_print,control_psi_tmp
      character*4 ion_aname
      integer     control_it_in, control_it_out,control_minimizer
      integer     ion_nion,ion_katm,paw_psi_ne,paw_psi_ispin
      integer     control_rtdb
      real*8      control_tole,control_tolc
      real*8      ion_rion,paw_psi_eigenvalue,ion_amass
      external control_out_of_time
      external ion_aname
      external control_it_in, control_it_out,control_minimizer
      external ion_nion,ion_katm,paw_psi_ne,paw_psi_ispin,control_rtdb
      external control_print,control_psi_tmp
      external control_tole,control_tolc
      external ion_rion,paw_psi_eigenvalue,ion_amass


      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))

      psitmp = control_psi_tmp()

*     **** set the minimizer ****
      minimizer = control_minimizer()

*     **** generate phaze factors and local psp and core density ****
      call paw_electron_gen_vl_potential()


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
      if (minimizer.gt.1) call paw_Grsm_list_start()
      if (minimizer.eq.5) it_out = 1
      if (newpsi) then
        call paw_sdminimize(15)
        !call bybminimize(E,deltae,deltac,1,.true.,2)
      end if
   2  continue
         icount = icount + 1
         if (stalled) then
           call paw_sdminimize(0)
           bfgscount = 0 
         end if

         deltae_old = deltae
         if (minimizer.eq.1) then
           bfgscount = bfgscount + 1
           call paw_cgminimize(E,deltae,deltac,bfgscount,it_in)
         else if (minimizer.eq.2) then
           bfgscount = bfgscount + 1
           call paw_bfgsminimize(E,deltae,deltac,bfgscount,minimizer)
         else if (minimizer.eq.3) then
           bfgscount = bfgscount + 1
           !call bfgsminimize(E,deltae,deltac,bfgscount,minimizer)
         else if (minimizer.eq.4) then
           bfgscount = bfgscount + 1
           !call cgminimize2(E,deltae,deltac,bfgscount)
         else if (minimizer.eq.5) then
           bfgscount = bfgscount + 1
           !call bybminimize(E,deltae,deltac,bfgscount,.false.,1)
         else if (minimizer.eq.6) then
           bfgscount = bfgscount + 1
           !call bfgsminimize3(E,deltae,deltac,bfgscount,minimizer)
         else if (minimizer.eq.7) then
           bfgscount = bfgscount + 1
           !call bfgsminimize2(E,deltae,deltac,bfgscount,minimizer)
         end if

         if ((dabs(deltae).gt.dabs(deltae_old)).or.
     >       (dabs(deltae).gt.1.0d-2)          .or.
     >       (deltae.gt.0.0d0)) then
            stalled = .true.
         else
            stalled = .false.
         end if
        
         if (psitmp) then
            call psi_tmp_write()
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
      if (minimizer.gt.1) call paw_Grsm_list_end()
      if (oprint) CALL NWPW_MESSAGE(3)


*     **** diagonalize hamiltonian and rotate psi  ****
*     **** but skip if doing an SIC calculations   ****
c     BUG!!!!!!
      call paw_psi_spin_density(en1,en2)
      en(1) = en1(1)+en2(1)
      en(2) = en1(2)+en2(2)

      call paw_psi_1gen_hml()

      call paw_psi_diagonalize_hml()
      call paw_psi_1rotate2()
      call paw_psi_2to1()

 

*     **** geometrical center and center of mass of the cluster ****
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)



*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (oprint) then
         neall = (paw_psi_ne(1)+paw_psi_ne(2))
         write(6,1304)
         write(6,1410)
         write(6,*)
c         write(6,1321) r_charge,' (real space)'
         write(6,1320) en(1), en(paw_psi_ispin()),' (real space)'
         write(6,1322) en1(1),en1(paw_psi_ispin()),' (real space)'
         write(6,1323) en2(1),en2(paw_psi_ispin()),' (real space)'

         call paw_energy_output(luout,ion_nion(),neall,E)

*        **** write out KS eigenvalues ****
         write(6,1500)
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


      paw_cgsd_energy = E(1)
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

*     ***********************************
*     *                                 *
*     *     paw_cgsd_energy_gradient    *
*     *                                 *
*     ***********************************

      subroutine paw_cgsd_energy_gradient(G1)
      implicit none
      real*8 G1(3,*)

#include "stdio.fh"
#include "util.fh"

      logical allow_translation,lprint,mprint
      integer MASTER
      parameter (MASTER=0)
      integer i,k,taskid,nion,nion1
      integer i1
      real*8  GG
      real*8  fmx,fmy,fmz
      real*8  fmx2,fmy2,fmz2

*     **** external functions ****
      logical     control_allow_translation,ion_q_FixIon,control_print
      character*4 ion_aname
      integer     ion_katm,ion_nion,control_version
      real*8      ion_rion
      external control_allow_translation,ion_q_FixIon,control_print
      external ion_aname
      external ion_katm,ion_nion,control_version
      external ion_rion


      nion = ion_nion()
      call dcopy(3*nion,0.0d0,0,G1,1)
      call paw_set_mult_energy_coeff()
      call paw_psi_1force(G1)

*     **** remove ion forces using ion_FixIon ****
      allow_translation = control_allow_translation()
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
      do i=1,nion
         GG = GG + G1(1,i)**2 + G1(2,i)**2 + G1(3,i)**2
      end do

      call Parallel_taskid(taskid)
      mprint = ((taskid.eq.MASTER).and.control_print(print_medium))
      lprint = ((taskid.eq.MASTER).and.control_print(print_low))

      if (taskid.eq.MASTER) then
        if (mprint) then
        !write(luout,1300)
        write(luout,1301)
        !write(luout,1302)
        end if

        if (lprint) then
        write(luout,1304)
        if (.not.allow_translation) write(luout,1400) fmx,fmy,fmz
        write(luout,1304)
        write(luout,1410)
        end if

        if (mprint) then
        write(luout,1420)
        do I=1,ion_nion()
          if (ion_q_FixIon(I)) then
           write(6,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3)
          else
           write(6,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3)
          end if
        end do
        end if

        if (lprint) then
        write(luout,1421)
        write(luout,1190)(i,ion_aname(I),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())


        write(luout,1210) fmx2,fmy2,fmz2  
        write(luout,1425)
        write(luout,1426) dsqrt(GG)
        end if
      end if

c     call dscal(3*nion,(-1.0d0),G1,1)

      return
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1191 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - fixed')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'========================')
 1301 FORMAT(//'== Gradient Calculation ==')
 1302 FORMAT(  '========================')
 1304 FORMAT(/)
 1400 FORMAT('Translation force removed: (',3F11.5,')')
 1410 FORMAT(10X,'=============  Ion Gradients =================')
 1425 FORMAT(10X,'===============================================')
 1426 FORMAT(10X,'|F|=',E21.10//)
 1420 FORMAT( ' Ion Positions:')
 1421 FORMAT( ' Ion Forces:')
      end


