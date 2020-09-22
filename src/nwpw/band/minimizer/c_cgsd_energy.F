
*     $Id$                       
*

*     ***************************
*     *				*
*     *	    c_cgsd_energy	*
*     *				*
*     ***************************

      real*8 function c_cgsd_energy(newpsi)
      implicit none
      logical newpsi

#include "stdio.fh"
#include "util.fh"

      integer MASTER
      parameter (MASTER=0)

      logical stalled,mprint
      integer taskid,minimizer,rtdb
      integer i,k,neall,NN,nb
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2)
      real*8  E(50),f0,f1,f2,f3,f4,f5,f6

*     **** external functions ****
      logical     cpsp_semicore,control_out_of_time
      character*4 ion_aname
      integer     control_it_in,control_it_out,control_minimizer
      integer     ion_nion,ion_katm,cpsi_ne,cpsi_ispin
      real*8      control_tole,control_tolc
      real*8      ion_rion,cpsi_eigenvalue_brdcst,ion_amass
      real*8      cpsi_svector_brdcst
      real*8      cpsi_occupation_brdcst
      external cpsp_semicore,control_out_of_time
      external ion_aname
      external control_it_in,control_it_out,control_minimizer
      external ion_nion,ion_katm,cpsi_ne,cpsi_ispin
      external control_tole,control_tolc
      external ion_rion,cpsi_eigenvalue_brdcst,ion_amass
      external cpsi_svector_brdcst
      external cpsi_occupation_brdcst
      logical  brillioun_print,control_print
      integer  brillioun_nbrillioun
      real*8   brillioun_k_brdcst,brillioun_ks_brdcst
      real*8   brillioun_weight_brdcst
      external brillioun_print,control_print
      external brillioun_nbrillioun
      external brillioun_k_brdcst,brillioun_ks_brdcst
      external brillioun_weight_brdcst
      logical  cpsi_spin_orbit,control_fractional,band_hfx
      external cpsi_spin_orbit,control_fractional,band_hfx
      real*8   cpsi_smearfermi
      external cpsi_smearfermi
      integer  control_rtdb
      external control_rtdb
      logical  control_DOS,ion_disp_on
      external control_DOS,ion_disp_on

      real*8   cpsi_1ke,cpsi_1energy
      external cpsi_1ke,cpsi_1energy

      call Parallel_taskid(taskid)

      mprint = (taskid.eq.MASTER).and.control_print(print_medium)

*     **** set the minimizer ****
      minimizer = control_minimizer()


*     **** generate phaze factors and local psp and core density ****
      call c_electron_gen_vl_potential()
      if (cpsp_semicore(0)) call c_semicore_density_update()


*     :::::::::::  begin electron iteration  :::::::::::::::::::::::
      if (mprint) then
         write(luout,1301)
         write(luout,1304)
         if (minimizer.eq.1) write(luout,1305)
         if (minimizer.eq.2) write(luout,1306)
         if (minimizer.eq.5) write(luout,1309)
         if (minimizer.eq.8) write(luout,1311)
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
      if (mprint) then
        if ((minimizer.eq.5).or.(minimizer.eq.8)) then
          call nwpw_message(7)
        else
          call nwpw_message(2)
        end if
      end if
      if (mprint) call util_flush(luout)
      if (minimizer.eq.2) call band_Grsm_list_start()
      if ((minimizer.eq.5).or.(minimizer.eq.8)) it_out = 1

      if (newpsi) then
         call c_sdminimize(10)
         call c_bybminimize0()
         if (control_fractional()) then
            call cpsi_1toelectron()
            call cpsi_1gen_hml()
            call cpsi_diagonalize_hml()
            call cpsi_1rotate2()
            call cpsi_2to1()
            !call cpsi_1define_occupation(-1.0d0,.false.)
            call cpsi_1define_occupation(1.0d0,.false.)
         end if
      end if

   2  continue
         icount = icount + 1
         if (stalled) then
           call c_sdminimize(0)
           bfgscount = 0
         end if

         deltae_old = deltae
         
         if (minimizer.eq.1) then
           bfgscount = bfgscount + 1
           call c_cgminimize(E,deltae,deltac,bfgscount,it_in)

         else if (minimizer.eq.2) then
           bfgscount = bfgscount + 1
           call c_bfgsminimize(E,deltae,deltac,bfgscount,minimizer)

         else if (minimizer.eq.3) then
           bfgscount = bfgscount + 1
           write(*,*) "Bundled DIIS minimizer not implemented"
           stop
         else if (minimizer.eq.4) then
           bfgscount = bfgscount + 1
           write(*,*) "Bundled Stiefel minimizer not implemented"
           stop
         else if (minimizer.eq.5) then
           bfgscount = bfgscount + 1
           call c_bybminimize(E,deltae,deltac,bfgscount,.false.,1)
         else if (minimizer.eq.6) then
           bfgscount = bfgscount + 1
           write(*,*) "minimizer 6 not implemented"
           stop
         else if (minimizer.eq.7) then
           bfgscount = bfgscount + 1
           write(*,*) "minimizer 7 not implemented"
           stop
         else if (minimizer.eq.8) then
           bfgscount = bfgscount + 1
           call c_bybminimize2(E,deltae,deltac,bfgscount,.false.,1)

         end if

         if ((dabs(deltae).gt.dabs(deltae_old)).or.
     >       (dabs(deltae).gt.1.0d-2)          .or.
     >       (deltae.gt.0.0d0)) then
            stalled = .true.
         else
            stalled = .false.
         end if

         if ((mprint).and.(minimizer.ne.5)
     >                         .and.(minimizer.ne.8)) then 
           write(luout,1310) icount*it_in,E(1),deltae,deltac
           call util_flush(luout)
         end if
         if (deltae.gt.0.0d0) then
            if (mprint) 
     >       write(luout,*) 
     >       ' *** energy going up. iteration terminated.'
            go to 3
         end if
         deltae = dabs(deltae)
         if ((deltae.lt.tole).and.
     >       (deltac.lt.tolc)) then
            if (mprint) 
     >       write(luout,*) 
     >       ' *** tolerance ok.    iteration terminated.'
            go to 3
         end if
         if (control_out_of_time()) then
            if (mprint) 
     >       write(luout,*) 
     >       ' *** out of time.     iteration terminated.'
            go to 3
         end if
      if (icount.lt.it_out) go to 2
      if (mprint) 
     > write(luout,*) 
     >  '*** arrived at the Maximum iteration.   terminated.'

*     :::::::::::  end of electron iteration loop  :::::::::::::::::::::

   3  continue
      if (minimizer.eq.2) call band_Grsm_list_end()
      if (mprint) CALL nwpw_MESSAGE(3)

      call cpsi_spin_density(en)

*     **** diagonalize hamiltonian and rotate psi ****
      call cpsi_1gen_hml()
      call cpsi_diagonalize_hml()
      if (.not.control_fractional()) then
         call cpsi_1rotate2()
         call cpsi_2to1()
      end if
      !if (control_fractional()) call cpsi_1define_occupation(.false.,1)

*     *** generate svectors ***
      if (cpsi_spin_orbit()) call cpsi_1gen_svector()

*     **** geometrical center and center of mass of the cluster ****
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)

**:::::::::::::::::   report summary of results  :::::::::::::::::::::::
      if (mprint) then
         neall = (cpsi_ne(1)+cpsi_ne(2))
         write(luout,1304)
         write(luout,1410)
        
c         write(luout,1420)                 ! do not write out ion in epilogue
c         write(luout,1190)(i,ion_aname(I),
c     >                 (ion_rion(K,I),K=1,3),
c     >                 ion_amass(i)/1822.89d0,
c     >                I=1,ion_nion())
c         write(luout,1200) CX,CY,CZ
c         write(luout,1210) GX,GY,GZ

         write(luout,*)
         write(luout,1320) en(1),en(cpsi_ispin()),' (real space)'
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
         write(luout,1498) E(11),E(11)/neall
         if (band_HFX())  then
           write(luout,1502) E(27),E(27)/neall
         end if

         virial = (E(11)+E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1499) virial

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
c          if (mprint) write(luout,1511)
          do i=0,cpsi_ne(1)-1
            f1 = cpsi_eigenvalue_brdcst(nb,1,cpsi_ne(1)-i) 
            f2 = cpsi_svector_brdcst(nb,cpsi_ne(1)-i,1)
            f3 = cpsi_svector_brdcst(nb,cpsi_ne(1)-i,2)
            f4 = cpsi_svector_brdcst(nb,cpsi_ne(1)-i,3)
            f0 = dsqrt(f2*f2 + f3*f3 + f4*f4)
            f5 = cpsi_occupation_brdcst(nb,cpsi_ne(1)-i,1)       
            if (mprint) write(luout,1512) f1,f1*EV,f0,f2,f3,f4,f5
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
            if (mprint) write(luout,1510) f1,f1*EV,f2,f3,f3*EV,f4
          end do
        end if
      end do
      else
         if (mprint) write(luout,1501) 
      endif

      rtdb = control_rtdb()
      if (control_DOS()) call cpsi_DOS(rtdb)

      c_cgsd_energy = E(1)
      return



 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F6.3,' ')
 1200 FORMAT(5X,'   G.C.  ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'======================')
 1301 FORMAT(  '= energy calculation =')
 1302 FORMAT(  '======================')
 1304 FORMAT(/)
 1305 FORMAT(10X,'== bundled Grassmann conjugate gradient iteration ==')
 1306 FORMAT(10X,'======== bundled Grassmann lmbfgs iteration ========')
 1307 FORMAT(10X,'============ Grassmann diis iteration ==============')
 1308 FORMAT(10X,'=== bundled Stiefel conjugate gradient iteration ===')
 1309 FORMAT(10X,'=========== scf minimization (potential) ===========')
 1310 FORMAT(I8,E20.10,3E15.5)
 1311 FORMAT(10X,'=========== scf minimization (density)  ============')
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
 1498 FORMAT( ' K.S. V_so energy    :',E19.10,' (',E15.5,'/electron)')
 1499 FORMAT( ' Virial Coefficient  :',E19.10)
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
 1510 FORMAT(4(E18.7,' (',F8.3,'eV) occ=',F5.3))
 1511 FORMAT(33x,"Spin(Sz,Sy,Sz)")
c 1512 FORMAT(E18.7,' (',F8.3,' eV) (',F8.3,' zhat + ',F8.3,' yhat +',
c     > F8.3,' xhat ) occ=',F5.3)
 1512 FORMAT(E18.7,' (',F8.3,' eV) (|s| =',F6.3,
     >       ', s = <',F7.3,',',F7.3,',',F7.3,'> ) occ=',F5.3)

 1720 FORMAT(/' Dispersion energy   :',E19.10)

      end


