c     $Id: cgsd_energy.F,v 1.1 2001-08-30 00:25:11 edo Exp $                       
*     ***************************
*     *							*
*     *		cgsd_energy			*
*     *							*
*     ***************************

      real*8 function cgsd_energy()
      implicit none
      integer MASTER
      parameter (MASTER=0)

      integer taskid
      integer i,k,neall,NN
      integer it_in,it_out,icount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(20)

*     **** external functions ****
      logical     psp_semicore,Waterpsp_found
      character*2 psp_atom
      integer     control_it_in, control_it_out
      integer     ion_nion,ion_katm,psi_ne,psi_ispin
      real*8      control_tole,control_tolc
      real*8      ion_rion,psi_eigenvalue
      external psp_semicore,Waterpsp_found
      external psp_atom
      external control_it_in, control_it_out
      external ion_nion,ion_katm,psi_ne,psi_ispin
      external control_tole,control_tolc
      external ion_rion,psi_eigenvalue


      call Parallel_taskid(taskid)

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
         write(6,1305)
         call util_flush(6)
      end if

      icount=0
      it_in  = control_it_in()
      it_out = control_it_out()
      tole   = control_tole()
      tolc   = control_tolc()
      E(1)=0.0d0
      if (taskid.eq.MASTER) call message(2)
      if (taskid.eq.MASTER) call util_flush(6)
   2  continue
         icount = icount + 1
         call cgminimize(E,deltae,deltac)

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
         write(6,1190)(i,psp_atom(ion_katm(I)),
     >                  (ion_rion(K,I),K=1,3),I=1,ion_nion())
         write(6,1200) CX,CY,CZ
         write(6,1210) GX,GY,GZ

         call Waterpsp_Print(6)

         write(6,*)
         write(6,1320) en(1),en(psi_ispin()),' (real space)'
         write(6,1430) E(1),E(1)/ion_nion()
         write(6,1440) E(2),E(2)/neall
         write(6,1450) E(3),E(3)/neall
         write(6,1460) E(4),E(4)/neall
         write(6,1470) E(5),E(5)/ion_nion()

         write(6,1480) E(6),E(6)/neall
         if (Waterpsp_found()) then
            write(6,1491) E(7),E(7)/neall
         else
            write(6,1490) E(7),E(7)/neall
         end if
         write(6,1495) E(8),E(8)/neall
         write(6,1496) E(9),E(9)/neall
         write(6,1497) E(10),E(10)/neall
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(6,1498) virial

         if (Waterpsp_found()) then
            write(6,1701) E(12),E(12)/ion_nion()
            write(6,1702) E(13)
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



 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
 1200 FORMAT(5X,'  G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'======================')
 1301 FORMAT(  '= energy calculation =')
 1302 FORMAT(  '======================')
 1304 FORMAT(/)
 1305 FORMAT(10X,'============ electron iteration ====================')
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
 1491 FORMAT( ' K.S. Vl+Vwater      :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1500 FORMAT(/' orbital energies:')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1700 FORMAT(/' water-elc energy    :',E19.10,' (',E15.5,'/electron)')
 1701 FORMAT(/' water-ion energy    :',E19.10,' (',E15.5,'/ion)')
 1702 FORMAT( ' water-water energy  :',E19.10)
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
      integer i1,i2,i3
      real*8  GG
      real*8  fmx,fmy,fmz
      real*8  fmx2,fmy2,fmz2

*     **** external functions ****
      logical     psp_semicore,Waterpsp_found,control_allow_translation
      character*2 psp_atom
      integer     ion_katm,ion_nion,control_version
      integer     Waterpsp_nwater
      real*8      ion_rion,Waterpsp_rwater
      external psp_semicore,Waterpsp_found,control_allow_translation
      external psp_atom
      external ion_katm,ion_nion,control_version
      external Waterpsp_nwater
      external ion_rion,Waterpsp_rwater


      allow_translation = control_allow_translation()
      nion = ion_nion()
      if (Waterpsp_found()) nion = nion + 3*Waterpsp_nwater()

      call dcopy(3*nion,0.0d0,0,G1,1)

      call psi_1force(G1)
      if (psp_semicore(0)) call electron_semicoreforce(G1)

      if (control_version().eq.3) call ewald_f(G1)
      if (control_version().eq.4) call ion_ion_f(G1)
      if (Waterpsp_found()) then
         call Waterpsp_Fion(G1)
         call rho_1Fwater(G1(1,ion_nion()+1))
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
        write(6,1190)(i,psp_atom(ion_katm(I)),
     >                  (ion_rion(K,I),K=1,3),I=1,ion_nion())
        if (Waterpsp_found()) then
          do i=1,Waterpsp_nwater()
            i1 = ion_nion() + 3*(i-1) + 1
            i2 = i1+1
            i3 = i1+2
            write(6,1190) i1,'O^',(Waterpsp_rwater(K,1,i),K=1,3)
            write(6,1190) i2,'H^',(Waterpsp_rwater(K,2,i),K=1,3)
            write(6,1190) i3,'H^',(Waterpsp_rwater(K,3,i),K=1,3)
          end do
        end if

        write(6,1421)
        write(6,1190)(i,psp_atom(ion_katm(I)),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())

        if (Waterpsp_found()) then
          do i=1,Waterpsp_nwater()
            i1 = ion_nion() + 3*(i-1) + 1
            i2 = i1+1
            i3 = i1+2
            write(6,1190) i1,'O^',(G1(K,i1),K=1,3)
            write(6,1190) i2,'H^',(G1(K,i2),K=1,3)
            write(6,1190) i3,'H^',(G1(K,i3),K=1,3)
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
c     write(*,*) "KE "
c     call write_stress(stress)


*     **** Coulomb energy component : dE_Coul/dhuv ****
      call psi_1coulomb_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c     write(*,*) "+coulomb"
c     call write_stress(stress)

*     **** Core-correction Coulomb energy component : dE_core/dhuv ****
*     call electron_semicore_stress(tstress)
*     call daxpy(9,1.0d0,tstress,1,stress,1)

*     **** Local pseudo energy component : dE_local/dhuv ****
      call dng_1vlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c     write(*,*) "+local"
c     call write_stress(stress)

*     **** Nonlocal pseudo energy component : dE_nolocal/dhuv ****
      call psi_1vnonlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c     write(*,*) "nonlocal"
c     call write_stress(tstress)

*     **** xc energy component : dE_xc/dhuv ****
      call rho_1exc_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c     write(*,*) "+xc"
c     call write_stress(stress)

*     **** Ewald energy component : dE_ewald/dhuv ****
      call ewald_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c     write(*,*) "+ewald"
c     call write_stress(stress)

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

