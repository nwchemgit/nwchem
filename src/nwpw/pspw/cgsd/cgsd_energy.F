c     $Id: cgsd_energy.F,v 1.7 2001-12-21 16:50:29 bylaska Exp $                       
*     ***************************
*     *							*
*     *		cgsd_energy			*
*     *							*
*     ***************************

      real*8 function cgsd_energy(newpsi)
      implicit none
      logical newpsi
      integer MASTER
      parameter (MASTER=0)

      logical stalled
      integer taskid
      integer minimizer
      integer i,k,neall,NN
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(20)

*     **** external functions ****
      logical     psp_semicore,pspw_qmmm_found
      character*2 ion_aname
      integer     control_it_in, control_it_out,control_minimizer
      integer     ion_nion,ion_katm,psi_ne,psi_ispin
      real*8      control_tole,control_tolc
      real*8      ion_rion,psi_eigenvalue,ion_amass
      external psp_semicore,pspw_qmmm_found
      external ion_aname
      external control_it_in, control_it_out,control_minimizer
      external ion_nion,ion_katm,psi_ne,psi_ispin
      external control_tole,control_tolc
      external ion_rion,psi_eigenvalue,ion_amass


      call Parallel_taskid(taskid)

*     **** set the minimizer ****
      minimizer = control_minimizer()

*     **** generate phaze factors and local psp and core density ****
      call phafac()
      call electron_gen_vl_potential()
      if (psp_semicore(0)) call semicore_density_update()


*     :::::::::::  begin electron iteration  :::::::::::::::::::::::
      if (taskid.eq.MASTER) then
         write(6,1300)
         write(6,1301)
         write(6,1302)
         write(6,1304)
         if (minimizer.eq.1) write(6,1305)
         if (minimizer.eq.2) write(6,1306)
         call util_flush(6)
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
      if (taskid.eq.MASTER) call message(2)
      if (taskid.eq.MASTER) call util_flush(6)
      if (newpsi) call sdminimize(5)
   2  continue
         icount = icount + 1
         if (stalled) then
           call sdminimize(0)
           bfgscount = 0 
         end if

         deltae_old = deltae
         if (minimizer.eq.1) then
           call cgminimize(E,deltae,deltac)
         else if (minimizer.eq.2) then
           bfgscount = bfgscount + 1
           call bfgsminimize(E,deltae,deltac,bfgscount)
         end if

         if ((dabs(deltae).gt.dabs(deltae_old)).or.
     >       (dabs(deltae).gt.1.0d-3)) then
            stalled = .true.
         else
            stalled = .false.
         end if

         if (taskid.eq.MASTER) then 
           write(6,1310) icount*it_in,E(1),deltae,deltac
           call util_flush(6)
         end if
         if (deltae.gt.0.0d0) then
            if (taskid.eq.MASTER) 
     >       write(6,*) ' *** Energy going up.  iteration terminated.'
            go to 3
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.tole).and.
     >       (deltac.lt.tolc)) then
            if (taskid.eq.MASTER) 
     >       write(6,*) ' *** tolerance ok.     iteration terminated.'
            go to 3
         end if
      if (icount.lt.it_out) go to 2
      if (taskid.eq.MASTER) 
     > write(6,*) '*** arived at the Maximum iiteration.   terminated.'

*     :::::::::::  end of electron iteration loop  :::::::::::::::::::::

   3  continue
      if (taskid.eq.MASTER) CALL MESSAGE(3)


*     **** diagonalize hamiltonian and rotate psi ****
      call psi_spin_density(en)
      call psi_1gen_hml()
      call psi_diagonalize_hml()

      call psi_1rotate2()
      call psi_2to1()

*     **** geometrical center and center of mass of the cluster ****
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

*:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (taskid.eq.MASTER) then
         neall = (psi_ne(1)+psi_ne(2))
         write(6,1304)
         write(6,1410)
         write(6,1420)
         write(6,1190)(i,ion_aname(I),
     >                 (ion_rion(K,I),K=1,3),
     >                 ion_amass(i)/1822.89d0,
     >                I=1,ion_nion())
         write(6,1200) CX,CY,CZ
         write(6,1210) GX,GY,GZ

         call pspw_qmmm_Print(6)

         write(6,*)
         write(6,1320) en(1),en(psi_ispin()),' (real space)'
         write(6,1430) E(1),E(1)/ion_nion()
         write(6,1440) E(2),E(2)/neall
         write(6,1450) E(3),E(3)/neall
         write(6,1460) E(4),E(4)/neall
         write(6,1470) E(5),E(5)/ion_nion()

         write(6,1480) E(6),E(6)/neall
         if (pspw_qmmm_found()) then
            write(6,1491) E(7),E(7)/neall
         else
            write(6,1490) E(7),E(7)/neall
         end if
         write(6,1495) E(8),E(8)/neall
         write(6,1496) E(9),E(9)/neall
         write(6,1497) E(10),E(10)/neall
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(6,1498) virial

        if (pspw_qmmm_found()) then
            write(6,1700)
            write(6,1701)
            write(6,1702) E(11)
            write(6,1703) E(12)
            write(6,1704) E(13)
            write(6,1705) E(14)
            write(6,1706) E(15)
         end if


         write(6,1500)
         NN=psi_ne(1)-psi_ne(2)
         EV=27.2116
         DO I=1,NN
           write(6,1510) psi_eigenvalue(1,I),
     >                   psi_eigenvalue(1,I)*EV
         end do
         DO I=1,psi_ne(2)
           write(6,1510)  psi_eigenvalue(1,I+NN),
     >                    psi_eigenvalue(1,I+NN)*EV,
     >                    psi_eigenvalue(2,I),
     >                    psi_eigenvalue(2,I)*EV
         end do

      end if

      cgsd_energy = E(1)
      return



 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'======================')
 1301 FORMAT(  '= energy calculation =')
 1302 FORMAT(  '======================')
 1304 FORMAT(/)
 1305 FORMAT(10X,'========== conjugate gradient iteration ============')
 1306 FORMAT(10X,'============== lmbfgs iteration ====================')
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
 1460 FORMAT( ' exc-corr  energy    :',E19.10,' (',E15.5,'/electron)')
 1470 FORMAT( ' ion-ion   energy    :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT(/' K.S. kinetic energy :',E19.10,' (',E15.5,'/electron)')
 1490 FORMAT( ' K.S. V_l  energy    :',E19.10,' (',E15.5,'/electron)')
 1491 FORMAT( ' K.S. Vl+Vqm/mm      :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' QM/MM energy           :',E19.10)
 1703 FORMAT( ' MM/MM energy           :',E19.10)
 1704 FORMAT( ' MM Polarization energy :',E19.10)
 1705 FORMAT( ' MM Vibration energy    :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy  :',E19.10)
      end

*     *******************************
*     *								*
*     *	    cgsd_energy_gradient	*
*     *								*
*     *******************************

      subroutine cgsd_energy_gradient(G1)
      implicit none
      real*8 G1(3,*)

      logical allow_translation
      integer MASTER
      parameter (MASTER=0)
      integer i,k,taskid,nion
      integer i1
      real*8  GG
      real*8  fmx,fmy,fmz
      real*8  fmx2,fmy2,fmz2

*     **** external functions ****
      logical     psp_semicore,pspw_qmmm_found
      logical     control_allow_translation
      character*2 ion_aname,pspw_qmmm_aname
      integer     ion_katm,ion_nion,control_version
      integer     pspw_qmmm_nion
      real*8      ion_rion,pspw_qmmm_rion
      external psp_semicore,pspw_qmmm_found
      external control_allow_translation
      external ion_aname,pspw_qmmm_aname
      external ion_katm,ion_nion,control_version
      external pspw_qmmm_nion
      external ion_rion,pspw_qmmm_rion


      allow_translation = control_allow_translation()
      nion = ion_nion()
      if (pspw_qmmm_found()) nion = nion + pspw_qmmm_nion()

      call dcopy(3*nion,0.0d0,0,G1,1)

      call psi_1force(G1)
      if (psp_semicore(0)) call electron_semicoreforce(G1)

      if (control_version().eq.3) call ewald_f(G1)
      if (control_version().eq.4) call ion_ion_f(G1)
      if (pspw_qmmm_found()) then
         call pspw_qmmm_Fion(G1)
         call rho_1Fmm(G1(1,ion_nion()+1))
      end if

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
      if (taskid.eq.MASTER) then
        write(6,1300)
        write(6,1301)
        write(6,1302)
        write(6,1304)
        if (.not.allow_translation) write(6,1400) fmx,fmy,fmz
        write(6,1304)
        write(6,1410)
        write(6,1420)
        write(6,1190)(i,ion_aname(I),
     >                  (ion_rion(K,I),K=1,3),I=1,ion_nion())
        if (pspw_qmmm_found()) then
          do i=1,pspw_qmmm_nion()
            i1 = ion_nion() + i
            write(6,1190) i1,pspw_qmmm_aname(i),
     >                    (pspw_qmmm_rion(K,i),K=1,3)
          end do
        end if

        write(6,1421)
        write(6,1190)(i,ion_aname(I),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())

        if (pspw_qmmm_found()) then
          do i=1,pspw_qmmm_nion()
            i1 = ion_nion() + i
            write(6,1190) i1,pspw_qmmm_aname(i),
     >                    (G1(K,i1),K=1,3)
          end do
        end if

        write(6,1210) fmx2,fmy2,fmz2  
        write(6,1425)
        write(6,1426) dsqrt(GG)
      end if

c     call dscal(3*nion,(-1.0d0),G1,1)

      return
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'========================')
 1301 FORMAT(  '= Gradient calculation =')
 1302 FORMAT(  '========================')
 1304 FORMAT(/)
 1400 FORMAT('Translation force removed: (',3F11.5,')')
 1410 FORMAT(10X,'=============  Ion Gradients =================')
 1425 FORMAT(10X,'===============================================')
 1426 FORMAT(10X,'|F|=',E21.10//)
 1420 FORMAT( ' Ion Positions:')
 1421 FORMAT( ' Ion Forces:')
      end


*     ***************************
*     *							*
*     *	    cgsd_energy_stress	*
*     *							*
*     ***************************

      subroutine cgsd_energy_stress(stress)
      implicit none
      real*8 stress(3,3)

      integer taskid,MASTER
      parameter (MASTER=0)

*     **** local variables ****
      integer u,v,s
      real*8  tstress(3,3),ht(3,3),scal

*     **** external functions ****
      real*8   lattice_omega,lattice_unita
      external lattice_omega,lattice_unita


      call dcopy(9,0.0d0,0,stress,1)

*     **** Kinetic energy component : dE_kin/dhuv ****
      call psi_1ke_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      write(*,*) "KE "
      call write_stress(tstress)


*     **** Coulomb energy component : dE_Coul/dhuv ****
      call psi_1coulomb_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      write(*,*) "+coulomb"
      call write_stress(tstress)

*     **** Core-correction Coulomb energy component : dE_core/dhuv ****
*     call electron_semicore_stress(tstress)
*     call daxpy(9,1.0d0,tstress,1,stress,1)

*     **** Local pseudo energy component : dE_local/dhuv ****
      call dng_1vlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      write(*,*) "+local"
      call write_stress(tstress)

*     **** Nonlocal pseudo energy component : dE_nolocal/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      call psi_1vnonlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      write(*,*) "nonlocal"
      call write_stress(tstress)

*     **** xc energy component : dE_xc/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      write(*,*) " xc 0"
      call write_stress(tstress)
      call rho_1exc_stress(tstress)
      write(*,*) " xc a"
      call write_stress(tstress)

      call daxpy(9,1.0d0,tstress,1,stress,1)
      write(*,*) " xc"
      call write_stress(tstress)

*     **** Ewald energy component : dE_ewald/dhuv ****
      call ewald_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      write(*,*) "+ewald"
      call write_stress(tstress)

**     **** define ht ****
*      do v=1,3
*      do u=1,3
*         ht(u,v)=lattice_unita(v,u)
*      end do
*      end do
*
**     **** define stress tensor ****
*      call dcopy(9,0.0d0,0,tstress,1)
*      scal = -1.0d0/lattice_omega()
*      do v=1,3
*      do u=1,3
*         do s=1,3
*            tstress(u,v) = tstress(u,v) + scal*stress(u,s)*ht(s,v)
*         end do
*      end do
*      end do
*      call dcopy(9,tstress,1,stress,1)

      call dscal(9,-1.0d0,stress,1)
      call write_stress(stress)

      return
      end

      subroutine write_stress(stress)
      implicit none
      real*8 stress(3,3)

      integer MASTER,taskid
      parameter (MASTER=0)

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         write(*,1900)
         write(*,1901)
         write(*,1902)
         write(*,1904)
         write(*,1910)
         write(*,1911) stress(1,1),stress(1,2),stress(1,3)
         write(*,1912) stress(2,1),stress(2,2),stress(2,3)
         write(*,1912) stress(3,1),stress(3,2),stress(3,3)
         write(*,1915)
         write(*,1916) dsqrt(
     >                  stress(1,1)**2+stress(1,2)**2+stress(1,3)**2
     >                 +stress(2,1)**2+stress(2,2)**2+stress(2,3)**2
     >                 +stress(3,1)**2+stress(3,2)**2+stress(3,3)**2)

 1900 FORMAT(//'======================')
 1901 FORMAT(  '= Stress calculation =')
 1902 FORMAT(  '======================')
 1904 FORMAT(/)
 1910 FORMAT(10X,'=============  Stress Tensor =================')
 1911 FORMAT(5X,' S = ',' (',3F11.5,' )')
 1912 FORMAT(5X,'     ',' (',3F11.5,' )')
 1915 FORMAT(10X,'==============================================')
 1916 FORMAT(10X,'|S| = ',E11.5)
      end if

      return
      end

