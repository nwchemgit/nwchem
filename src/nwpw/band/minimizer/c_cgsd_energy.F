*
*     $Id: c_cgsd_energy.F,v 1.1 2001-12-17 23:31:13 bylaska Exp $                       
*

*     ***************************
*     *							*
*     *		c_cgsd_energy			*
*     *							*
*     ***************************

      real*8 function c_cgsd_energy(newpsi)
      implicit none
      logical newpsi
      integer MASTER
      parameter (MASTER=0)

      logical stalled
      integer taskid
      integer i,k,neall,NN,nb
      integer it_in,it_out,icount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(20)

*     **** external functions ****
      logical     cpsp_semicore
      character*2 ion_aname
      integer     control_it_in, control_it_out
      integer     ion_nion,ion_katm,cpsi_ne,cpsi_ispin
      real*8      control_tole,control_tolc
      real*8      ion_rion,cpsi_eigenvalue,ion_amass
      external cpsp_semicore
      external ion_aname
      external control_it_in, control_it_out
      external ion_nion,ion_katm,cpsi_ne,cpsi_ispin
      external control_tole,control_tolc
      external ion_rion,cpsi_eigenvalue,ion_amass
      integer  brillioun_nbrillioun
      real*8   brillioun_k,brillioun_ks,brillioun_weight
      external brillioun_nbrillioun
      external brillioun_k,brillioun_ks,brillioun_weight


      call Parallel_taskid(taskid)


*     **** generate phaze factors and local psp and core density ****
      call cphafac()
      call c_electron_gen_vl_potential()
      if (cpsp_semicore(0)) call c_semicore_density_update()


*     :::::::::::  begin electron iteration  :::::::::::::::::::::::
      if (taskid.eq.MASTER) then
         write(6,1300)
         write(6,1301)
         write(6,1302)
         write(6,1304)
         write(6,1305)
         call util_flush(6)
      end if

      stalled = .false.
      deltae  = -1.0d-03
      icount=0
      it_in  = control_it_in()
      it_out = control_it_out()
      tole   = control_tole()
      tolc   = control_tolc()
      E(1)=0.0d0
      if (taskid.eq.MASTER) call message(2)
      if (taskid.eq.MASTER) call util_flush(6)
      if (newpsi) call c_sdminimize(5)
   2  continue
         icount = icount + 1
         if (stalled) call c_sdminimize(0)

         deltae_old = deltae
         call c_cgminimize(E,deltae,deltac)

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
      call cpsi_spin_density(en)
      call cpsi_1gen_hml()
      call cpsi_diagonalize_hml()

      call cpsi_1rotate2()
      call cpsi_2to1()

*     **** geometrical center and center of mass of the cluster ****
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

**:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (taskid.eq.MASTER) then
         neall = (cpsi_ne(1)+cpsi_ne(2))
         write(6,1304)
         write(6,1410)
         write(6,1420)
         write(6,1190)(i,ion_aname(I),
     >                 (ion_rion(K,I),K=1,3),
     >                 ion_amass(i)/1822.89d0,
     >                I=1,ion_nion())
         write(6,1200) CX,CY,CZ
         write(6,1210) GX,GY,GZ


         write(6,*)
         write(6,1320) en(1),en(cpsi_ispin()),' (real space)'
         write(6,1430) E(1),E(1)/ion_nion()
         write(6,1440) E(2),E(2)/neall
         write(6,1450) E(3),E(3)/neall
         write(6,1460) E(4),E(4)/neall
         write(6,1470) E(5),E(5)/ion_nion()

         write(6,1480) E(6),E(6)/neall
         write(6,1490) E(7),E(7)/neall
         write(6,1495) E(8),E(8)/neall
         write(6,1496) E(9),E(9)/neall
         write(6,1497) E(10),E(10)/neall
         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(6,1498) virial

         NN=cpsi_ne(1)-cpsi_ne(2)
         EV=27.2116
         do nb=1,brillioun_nbrillioun()
           write(6,1508) nb,
     >                   brillioun_weight(nb),
     >                   brillioun_ks(1,nb),
     >                   brillioun_ks(2,nb),
     >                   brillioun_ks(3,nb),
     >                   brillioun_k(1,nb),
     >                   brillioun_k(2,nb),
     >                   brillioun_k(3,nb)
           write(6,1500)
           do i=0,NN-1
             write(6,1510) cpsi_eigenvalue(nb,1,cpsi_ne(1)-i),
     >                     cpsi_eigenvalue(nb,1,cpsi_ne(1)-i)*EV
           end do
           do i=0,cpsi_ne(2)-1
             write(6,1510)  cpsi_eigenvalue(nb,1,cpsi_ne(1)-i-NN),
     >                      cpsi_eigenvalue(nb,1,cpsi_ne(1)-i-NN)*EV,
     >                      cpsi_eigenvalue(nb,2,cpsi_ne(2)-i),
     >                      cpsi_eigenvalue(nb,2,cpsi_ne(2)-i)*EV
           end do
         end do

      end if

      c_cgsd_energy = E(1)
      return



 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
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
 1491 FORMAT( ' K.S. Vl+Vqm/mm      :',E19.10,' (',E15.5,'/electron)')
 1495 FORMAT( ' K.S. V_nl energy    :',E19.10,' (',E15.5,'/electron)')
 1496 FORMAT( ' K.S. V_Hart energy  :',E19.10,' (',E15.5,'/electron)')
 1497 FORMAT( ' K.S. V_xc energy    :',E19.10,' (',E15.5,'/electron)')
 1498 FORMAT( ' Virial Coefficient  :',E19.10)
 1500 FORMAT(/' orbital energies:')
 1508 FORMAT(/' Brillion zone point: ',i2,
     >       /'    weight=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))

      end


