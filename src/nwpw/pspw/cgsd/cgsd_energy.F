c
c     $Id$                       
c

*     ***************************
*     *                         *
*     *        cgsd_energy      *
*     *                         *
*     ***************************

      real*8 function cgsd_energy(newpsi)
      implicit none

#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"
#include "errquit.fh"
c#include "global.fh"
#include "btdb.fh"

      logical newpsi
      integer MASTER,taskid
      parameter (MASTER=0)

      logical stalled,value,oprint,psitmp,failed,field_exist
      integer minimizer
      integer i,j,ms,neall,NN,nx,ny,nz,n2ft3d,tmp1(2)
      integer it_in,it_out,icount,bfgscount
      real*8  EV,virial,sumall,scal
      real*8  tole,tolc,deltae,deltac,deltae_old
      real*8  cx,cy,cz
      real*8  gx,gy,gz
      real*8  en(2),en1(2),en2(2)
      real*8  E(60),Egas
      integer rtdb

*     **** external functions ****
      logical     psp_semicore,pspw_charge_found,pspw_qmmm_found
      logical     psp_pawexist,psp_U_psputerm,meta_found
      logical     control_out_of_time,pspw_SIC,control_DOS,pspw_HFX
      logical     control_print,control_psi_tmp,control_fractional
      character*4 ion_aname
      integer     control_it_in, control_it_out,control_minimizer
      integer     ion_nion,ion_katm,psi_ne,psi_ispin,control_rtdb
      integer     psi_iptr_hml
      real*8      control_tole,control_tolc,psi_occupation
      real*8      ion_rion,psi_eigenvalue,psi_hml,ion_amass
      real*8      psi_smearfermi,lattice_omega
      external psp_semicore,pspw_charge_found,pspw_qmmm_found
      external psp_pawexist,psp_U_psputerm,meta_found
      external control_out_of_time,pspw_SIC,control_DOS,pspw_HFX
      external ion_aname
      external control_it_in, control_it_out,control_minimizer
      external ion_nion,ion_katm,psi_ne,psi_ispin,control_rtdb
      external psi_iptr_hml
      external control_print,control_psi_tmp,control_fractional
      external control_tole,control_tolc,psi_occupation
      external ion_rion,psi_eigenvalue,psi_hml,ion_amass
      external psi_smearfermi,lattice_omega
      character*4 psi_ab_irrep_name
      external    psi_ab_irrep_name
      integer     control_symmetry,control_version
      external    control_symmetry,control_version
      logical  ion_disp_on,nwpw_cosmo_on,nwpw_cosmo_firsttime
      external ion_disp_on,nwpw_cosmo_on,nwpw_cosmo_firsttime
      logical  nwpw_born_on
      external nwpw_born_on
      real*8   nwpw_cosmo_screen,nwpw_born_screen
      external nwpw_cosmo_screen,nwpw_born_screen

      logical  control_hess_model,pspw_Efield_found
      external control_hess_model,pspw_Efield_found
      real*8   ion_hess_energy
      external ion_hess_energy
      real*8   psi_1energy,control_gas_energy,psi_2energy
      external psi_1energy,control_gas_energy,psi_2energy
      logical  pspw_V_APC_on
      external pspw_V_APC_on
      logical  control_use_fractional_rho,psi_try_read_density
      external control_use_fractional_rho,psi_try_read_density
      integer  pspw_Efield_type
      external pspw_Efield_type
      real*8   ewald_e,psi_hml_value
      external ewald_e,psi_hml_value

      !**** hack for implementing hessian energy model ****
      if (control_hess_model()) then
         cgsd_energy = ion_hess_energy()
         return
      end if

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))

      psitmp = control_psi_tmp()

      field_exist = pspw_charge_found().or.pspw_Efield_found()


*     **** set the minimizer ****
      call dcopy(40,0.0d0,0,E,1)
      minimizer = control_minimizer()

*     **** generate phaze factors and local psp and core density ****
      call phafac()
      if (control_version().eq.3) call ewald_phafac()
      call electron_gen_vl_potential()
      if (psp_semicore(0)) call semicore_density_update()
      if (psp_pawexist())
     >   call psp_dE_ncmp_vloc_Qlm(psi_ispin(),.false.,en)


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
         if (minimizer.eq.8) write(luout,1313)
         if (minimizer.eq.9) write(luout,1315)
         if (minimizer.eq.10) write(luout,1316)
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
        if ((minimizer.eq.5).or.(minimizer.eq.8)) then
          call nwpw_message(7)
        else
          call nwpw_message(2)
        end if
      end if

      if (oprint) call util_flush(luout)
      if (minimizer.gt.1) call pspw_Grsm_list_start()
      if (minimizer.eq.5) it_out = 1
      if (minimizer.eq.8) it_out = 1
      if ((newpsi).or.(nwpw_cosmo_firsttime())) then
         call pspw_Lin_HFX_reset()
         if (minimizer.lt.4) call sdminimize(15)
         call bybminimize0()
         if (control_fractional()) then
            call psiocc_1to2()
            call psi_1toelectron()
            call psi_1gen_hml()
            call psi_diagonalize_hml_assending()
            call psi_1rotate2()
            call psi_2to1()
            call psi_1define_occupation(-1.0d0,.false.)
         end if
      end if

   2  continue
         icount = icount + 1
         if (stalled) then
           !call pspw_hfx_localize_set_findcolumns(.true.)
           call pspw_Lin_HFX_reset()
           call sdminimize(0)
           bfgscount = 0 
         end if


         deltae_old = deltae
         if (minimizer.eq.1) then
           bfgscount = bfgscount + 1
           call cgminimize(E,deltae,deltac,bfgscount,it_in)
         else if (minimizer.eq.2) then
           bfgscount = bfgscount + 1
           call bfgsminimize(E,deltae,deltac,bfgscount,minimizer)
         else if (minimizer.eq.3) then
           bfgscount = bfgscount + 1
           call bfgsminimize(E,deltae,deltac,bfgscount,minimizer)
         else if (minimizer.eq.4) then
           bfgscount = bfgscount + 1
           call cgminimize2(E,deltae,deltac,bfgscount)
           if (control_fractional()) then
              call psi_1gen_hml()
              call psi_diagonalize_hml_assending()
              call psi_1define_occupation(-1.0d0,.false.)
              !call psi_1rotate2()
              !call psi_2to1()
           end if
         else if (minimizer.eq.5) then
           bfgscount = bfgscount + 1
           call bybminimize(E,deltae,deltac,bfgscount,.false.,1,failed)

c          **** reset minimizer if failed ****
           if (failed) then
              minimizer = 1
              it_out    = control_it_out()
              call pspw_Grsm_list_end()
              if (oprint) then
                 write(luout,1314)
                 write(luout,*)
                 write(luout,1305)
                 call nwpw_message(2)
              end if
           end if

         else if (minimizer.eq.6) then
           bfgscount = bfgscount + 1
           call bfgsminimize3(E,deltae,deltac,bfgscount,minimizer)
         else if (minimizer.eq.7) then
           bfgscount = bfgscount + 1
           call bfgsminimize2(E,deltae,deltac,bfgscount,minimizer)

         else if (minimizer.eq.8) then
           bfgscount = bfgscount + 1
           value = .false.
           if (control_fractional()) then
              if (control_use_fractional_rho()) then
                 value = psi_try_read_density(2)
              end if
           end if
!$OMP PARALLEL 
           call bybminimize2(E,deltae,deltac,
     >                       bfgscount,.false.,1,value,failed)
!$OMP END PARALLEL 
           if (control_fractional()) then
              if (control_use_fractional_rho()) then
                 call psi_write_density(1)
              end if
           end if

c          **** reset minimizer if failed ****
           if (failed) then
              minimizer = 1
              it_out    = control_it_out()
              call pspw_Grsm_list_end()
              if (oprint) then
                 write(luout,1314)
                 write(luout,*)
                 write(luout,1305)
                 call nwpw_message(2)
              end if
           end if

         else if (minimizer.eq.9) then
           bfgscount = bfgscount + 1
           call cgminimizeABC(E,deltae,deltac,bfgscount,it_in)

         else if (minimizer.eq.10) then
           bfgscount = bfgscount + 1
           call bfgsminimizeABC(E,deltae,deltac,bfgscount,minimizer)
         end if

         if ((dabs(deltae).gt.dabs(deltae_old)).or.
     >       (dabs(deltae).gt.1.0d-2)          .or.
     >       (deltae.gt.0.0d0)) then
            stalled = .true.
         else
            stalled = .false.
         end if
         !stalled = .false. !debug
        
         if (psitmp) then
            call psi_tmp_write()
         end if

         if ((oprint).and.(minimizer.ne.5).and.(minimizer.ne.8)) then 
           write(luout,1310) icount*it_in,E(1),deltae,deltac
           call util_flush(luout)

           call ecce_print1('iteration counter', mt_int,icount*it_in,1)
           call ecce_print1('iterative total energy difference',
     >                      mt_dbl, deltae, 1)
           call ecce_print1('iterative total density difference',
     >                      mt_dbl,deltac, 1)
         end if
         if (deltae.gt.0.0d0) then
            if (oprint) 
     >       write(luout,*) 
     >       ' *** energy going up. iteration not terminated'
             stalled = .true.
             if (icount.lt.it_out) go to 2
         end if
         deltae = dabs(deltae)
         !if (deltae.lt.1.0d-3) then
         !   call pspw_hfx_localize_set_findcolumns(.false.)
         !end if
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
      if (minimizer.gt.1) call pspw_Grsm_list_end()
      if (oprint) CALL NWPW_MESSAGE(3)


*     **** diagonalize hamiltonian and rotate psi  ****
*     **** but skip if doing an SIC calculations   ****
      call psi_spin_density(en1)
      if (psp_pawexist()) then
         call D3dB_n2ft3d(1,n2ft3d)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal = lattice_omega()/dble(nx*ny*nz)
         !en2(1) = 0.0
         !en2(2) = 0.0
         if (.not.BA_push_get(mt_dbl,n2ft3d,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit(
     >        'cgsd_energy: out of stack memory',0,MA_ERR)
         call psi_1qlm_atom()
         do ms=1,psi_ispin()
           call nwpw_compcharge_gen_dn_cmp_smooth_ms(ms,dbl_mb(tmp1(1)))
           call Pack_c_unpack(0,dbl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
           call D3dB_r_Zero_Ends(1,dbl_mb(tmp1(1)))
           call D3dB_r_dsum(1,dbl_mb(tmp1(1)),sumall)
           en2(ms) = sumall*scal
         end do
         if (.not.BA_pop_stack(tmp1(2)))
     >   call errquit(
     >        'cgsd_energy: popping stack memory',0,MA_ERR)
      else
         en2(1) = 0.0d0
         en2(2) = 0.0d0
      end if
      en(1) = en1(1) + en2(1)
      en(2) = en1(2) + en2(2)
      call psi_1gen_hml()
      call psi_diagonalize_hml()

*     !*** make sure occupation is in descending order ***
      if (control_fractional())  then
         !call psi_1reverse_occupation()
         call psi_1desending_occupation()
      end if


      if (pspw_SIC()) then
         call psi_1gen_hml_g()
      else if (.not.control_fractional()) then
         call psi_1rotate2()
         call psi_2to1()
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
         if (psp_pawexist()) then
            write(luout,1322) en1(1),en1(psi_ispin()),' (real space)'
            write(luout,1323) en2(1),en2(psi_ispin()),' (real space)'
         end if

*       **** write APC potential and charges ***
        if (pspw_V_APC_on()) call pspw_shortprint_APC(luout)

         write(luout,1430) E(1),E(1)/ion_nion()


         if (pspw_qmmm_found()) then
            write(luout,1431)
            write(luout,1432)
            write(luout,1433)  E(1)-E(11),
     >                        (E(1)-E(11))/ion_nion()
         end if

         if (field_exist) then
            write(luout,1431)
            write(luout,1432)
            if (pspw_Efield_type().eq.0) then
               write(luout,1433) (E(1)-E(49)-E(51)),
     >         (E(1)-E(49)-E(51))/ion_nion()
            else
               write(luout,1433) (E(1)-E(49)-E(50)-E(51)),
     >         (E(1)-E(49)-E(50)-E(51))/ion_nion()
            end if
         end if

         !*** write-out total paw energy - includes core energies ***
         if (psp_pawexist()) 
     >      write(luout,1434) (E(1)+E(36)+E(45)),
     >                        (E(1)+E(36)+E(45))/ion_nion()

         write(luout,1440) E(2),E(2)/neall
         write(luout,1450) E(3),E(3)/neall
         write(luout,1460) E(4),E(4)/neall
         if (pspw_SIC()) then
           write(luout,1455) E(16),E(16)/neall
           write(luout,1456) E(17),E(17)/neall
         end if
         if (pspw_HFX()) then
           write(luout,1457) E(26),E(26)/neall
         end if
         if (psp_U_psputerm()) then
           write(luout,1458) E(29),E(29)/neall
         end if
         if (meta_found()) then
           write(luout,1459) E(31),E(31)/ion_nion()
         end if
         if (pspw_V_APC_on()) then
           write(luout,1505) E(22),E(22)/ion_nion()
         end if
         write(luout,1470) E(5),E(5)/ion_nion()
         if (control_fractional()) write(luout,1471) E(28),E(28)/neall



         write(luout,1480) E(6),E(6)/neall

c        if (pspw_qmmm_found()) then
c           write(luout,1491) E(7),E(7)/neall
c        else
c           write(luout,1490) E(7),E(7)/neall
c        end if
         write(luout,1490) E(7),E(7)/neall

         write(luout,1495) E(8),E(8)/neall
         write(luout,1496) E(9),E(9)/neall
         write(luout,1497) E(10),E(10)/neall
         if (pspw_SIC())  then
           write(luout,1499) E(18),E(18)/neall
           write(luout,1501) E(19),E(19)/neall
         end if
         if (pspw_HFX())  then
           write(luout,1502) E(27),E(27)/neall
         end if
         if (psp_U_psputerm()) then
           write(luout,1503) E(30),E(30)/neall
         end if
         if (meta_found()) then
           write(luout,1504) E(32),E(32)/neall
         end if
         if (pspw_V_APC_on()) then
           write(luout,1506) E(23),E(23)/neall
         end if

         virial = (E(10)+E(9)+E(8)+E(7))/E(6)
         write(luout,1498) virial

c        **** paw energies ****
         if (psp_pawexist()) then
            write(luout,1680) E(35),E(35)/neall
            write(luout,1681) E(36),E(36)/ion_nion()
            write(luout,1682) E(37),E(37)/neall

            write(luout,1683) E(38),E(38)/neall
            write(luout,1684) E(39),E(39)/neall

            write(luout,1685) E(40),E(40)/neall
            write(luout,1686) E(41),E(41)/neall
            write(luout,1687) E(42),E(42)/neall
            write(luout,1688) E(43),E(43)/neall
            write(luout,1689) E(44),E(44)/ion_nion()
            write(luout,1690) E(45),E(45)/ion_nion()
         end if

         if (pspw_qmmm_found()) then
            write(luout,1700)
            write(luout,1701)
            write(luout,1702) E(11)
            write(luout,1703) E(12)
            write(luout,1704) E(13)
            write(luout,1705) E(14)
c
c            write(luout,1703) E(14)+E(15)
c            write(luout,1710) E(14)
c            write (luout,1711) E(15)
c
c            write(luout,1704) E(16)
c            write(luout,1705) E(17)
c            write(luout,1706) E(18)
         end if
         if (ion_disp_on()) then
            write(luout,1720) E(33)
         end if

         if (field_exist) then
            write(luout,1800)
            write(luout,1801)
           if (pspw_Efield_type().eq.0) then
              write(luout,1805) E(49)+E(51)
              write(luout,1806) E(49)
              write(luout,1807) E(50)
              !write(luout,1804) E(51)
              call dipole_Efield_print(luout)
            else
               write(luout,1805) E(49)+E(50)+E(51)
               write(luout,1802) E(49)
               write(luout,1803) E(50)
               write(luout,1804) E(51)
            end if
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
        
*     **** write out COSMO energies ****
      if (nwpw_cosmo_on()) then
         Egas = control_gas_energy()
         if (oprint) then
            write(luout,1730)
            write(luout,1731) nwpw_cosmo_screen()
            write(luout,1732) E(46)
            write(luout,1733) E(47)
            write(luout,1734) E(48)
            write(luout,1735) E(46)+E(47)+E(48),
     >                        (E(46)+E(47)+E(48))*27.2116d0*23.06d0
            if (dabs(Egas).gt.1.0d-6)
     >         write(luout,1736)  E(1)-Egas,
     >                           (E(1)-Egas)*27.2116d0*23.06d0
            call nwpw_cosmo_print(luout,Egas,E(1))
         end if
      end if

*     **** write out extended Born solvation energies ****
      if (nwpw_born_on()) then
         if (pspw_V_APC_on()) then
            Egas = control_gas_energy()
            if (oprint) then
               write(luout,1740)
               write(luout,1741) nwpw_born_screen()
               write(luout,1745) E(22),E(22)*27.2116d0*23.06d0
               if (dabs(Egas).gt.1.0d-6)
     >            write(luout,1746) E(1)-Egas,
     >                             (E(1)-Egas)*27.2116d0*23.06d0
               call nwpw_born_print(luout,Egas,E(1))
            end if
         else
            call rho_1dngen_APC()
            call pspw_print_APC(luout)
         end if
      end if

      if (control_symmetry().gt.0) call psi_ab_kill_irrep_names()
      if (control_fractional()) call psi_1assending_occupation()

      rtdb = control_rtdb()
      if (control_DOS()) call psi_DOS(rtdb)

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

*     **** debug - energies written to rtdb for numerical stress ****
      rtdb = control_rtdb()
      value =           btdb_put(rtdb,'pspw:E_ke',      mt_dbl,1,E(6))
      value = value.and.btdb_put(rtdb,'pspw:E_hartree', mt_dbl,1,E(3))
      value = value.and.btdb_put(rtdb,'pspw:E_xc',      mt_dbl,1,E(4))
      value = value.and.btdb_put(rtdb,'pspw:E_ewald',   mt_dbl,1,E(5))
      value = value.and.btdb_put(rtdb,'pspw:E_local',   mt_dbl,1,E(7))
      value = value.and.btdb_put(rtdb,'pspw:E_nonlocal',mt_dbl,1,E(8))
      value = value.and.btdb_put(rtdb,'pspw:E_hsic',    mt_dbl,1,E(16))
      value = value.and.btdb_put(rtdb,'pspw:E_xcsic',   mt_dbl,1,E(17))
      value = value.and.btdb_put(rtdb,'pspw:E_hfx',     mt_dbl,1,E(26))
      value = value.and.btdb_put(rtdb,'pspw:E_dftu',    mt_dbl,1,E(29))
      value = value.and.btdb_put(rtdb,'pspw:E_disp',    mt_dbl,1,E(33))
      if (.not. value) call errquit(
     >   'cgsd_energy: numerical stress - error writing rtdb',0,
     >       RTDB_ERR)

*     **** debug - energies written to rtdb for numerical forces ****
      value = value.and.btdb_put(rtdb,'pspw:E_LJ',   mt_dbl,1,E(11))
      value = value.and.btdb_put(rtdb,'pspw:E_Q',    mt_dbl,1,E(12))
      value = value.and.btdb_put(rtdb,'pspw:E_spring', mt_dbl,1,E(13))
      if (.not. value) call errquit(
     >   'cgsd_noit_energy: numerical force - error writing rtdb',0,
     >       RTDB_ERR)

*     **** debug paw - energies written to rtdb for numerical forces ****
      value = value.and.btdb_put(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E(38))
      value = value.and.
     >        btdb_put(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E(40))
      value = value.and.btdb_put(rtdb,'pspw:E_cmp_cmp',  mt_dbl,1,E(41))
      value = value.and.btdb_put(rtdb,'pspw:E_cmp_pw',   mt_dbl,1,E(42))
      value = value.and.btdb_put(rtdb,'pspw:E_exc_atom', mt_dbl,1,E(43))
      value = value.and.btdb_put(rtdb,'pspw:E_pxc_atom', mt_dbl,1,E(44))
      value = value.and.btdb_put(rtdb,'pspw:kin_loc',  mt_dbl,1,E(35))
      value = value.and.btdb_put(rtdb,'pspw:kin_core', mt_dbl,1,E(36))
      value = value.and.btdb_put(rtdb,'pspw:val_core', mt_dbl,1,E(37))
      value = value.and.btdb_put(rtdb,'pspw:ion_core', mt_dbl,1,E(45))
      value = value.and.btdb_put(rtdb,'pspw:vloc_atom',mt_dbl,1,E(39))
      value = value.and.btdb_put(rtdb,'pspw:E_all',       mt_dbl,1,E(1))
      if (.not. value) call errquit(
     >   'cgsd_energy: numerical force - error writing rtdb',0,
     >       RTDB_ERR)

      call ecce_print1 ('total energy', mt_dbl, E(1), 1)
      call ecce_print1 ('nuclear repulsion energy', mt_dbl, E(5), 1)
      call ecce_print1 ('coulomb energy', mt_dbl, E(3), 1)
      call ecce_print1 ('exchange energy', mt_dbl, E(4), 1)
      call ecce_print1 ('correlation energy', mt_dbl, E(4), 1)


      cgsd_energy = E(1)
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
 1309 FORMAT(10X,'======= Kohn-Sham scf iteration (potential) ========')
 1310 FORMAT(I8,E20.10,3E15.5)
 1311 FORMAT(10X,'======= Projected Grassmann lmbfgs iteration =======')
 1312 FORMAT(10X,'============= Stiefel lmbfgs iteration =============')
 1313 FORMAT(10X,'======== Kohn-Sham scf iteration (density) =========')
 1314 FORMAT(10X,'======== Kohn-Sham scf iteration failed    =========')
 1315 FORMAT(10X,'==== Grassmann cg (Stich linesearch) iteration  ====')
 1316 FORMAT(10X,'===== Grassmann lmbfgs (Stich l.s.) iteration  =====')
 1320 FORMAT(' number of electrons: spin up=',F11.5,'  down=',F11.5,A)
 1321 FORMAT(' total charge of system:',F11.5,A)
 1322 FORMAT('     plane-wave part:         ',F11.5,'       ',F11.5,A)
 1323 FORMAT('      augmented part:         ',F11.5,'       ',F11.5,A)
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
 1434 FORMAT(/' total paw energy    :',E19.10,' (',E15.5,'/ion)')
 1440 FORMAT( ' total orbital energy:',E19.10,' (',E15.5,'/electron)')
 1450 FORMAT( ' hartree   energy    :',E19.10,' (',E15.5,'/electron)')
 1455 FORMAT( ' SIC-hartree energy  :',E19.10,' (',E15.5,'/electron)')
 1456 FORMAT( ' SIC-exc-corr energy :',E19.10,' (',E15.5,'/electron)')
 1457 FORMAT( ' HF exchange energy  :',E19.10,' (',E15.5,'/electron)')
 1458 FORMAT( ' DFT+U     energy    :',E19.10,' (',E15.5,'/electron)')
 1459 FORMAT( ' Metadynamics energy :',E19.10,' (',E15.5,'/electron)')
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
 1500 FORMAT(/' orbital energies:')
 1501 FORMAT( ' K.S. SIC-exc-corr energy :',E19.10,
     >        ' (',E15.5,'/electron)')
 1502 FORMAT( ' K.S. HFX energy     :',E19.10,
     >        ' (',E15.5,'/electron)')
 1503 FORMAT( ' K.S. DFT+U energy   :',E19.10,
     >        ' (',E15.5,'/electron)')
 1504 FORMAT( ' K.S. Metadynamics energy   :',E19.10,
     >        ' (',E15.5,'/ion)')
 1505 FORMAT( ' APC energy          :',E19.10,' (',E15.5,'/ion)')
 1506 FORMAT( ' K.S. V_APC energy   :',E19.10,' (',E15.5,'/ion)')
 1507 FORMAT(/' Fermi energy =',2(E18.7,' (',F8.3,'eV)'))
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))
 1511 FORMAT(2(E18.7,' (',F8.3,'eV)  occ=',F5.3))
 1512 FORMAT(2(E18.7,' (',F8.3,'eV)',A4))
 1513 FORMAT(2(E18.7,' (',F8.3,'eV)',A4,' occ=',F5.3))

 1680 FORMAT(/' kinetic (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1681 FORMAT( ' kinetic       (core):',E19.10,' (',E15.5,'/ion)')
 1682 FORMAT( ' valence-core        :',E19.10,' (',E15.5,'/electron)')
 1683 FORMAT( ' ncmp*V_local        :',E19.10,' (',E15.5,'/electron)')
 1684 FORMAT( ' V_lpaw  (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1685 FORMAT( ' coulomb (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1686 FORMAT( ' coulomb  (ncmp-ncmp):',E19.10,' (',E15.5,'/electron)')
 1687 FORMAT( ' coulomb    (ncmp-pw):',E19.10,' (',E15.5,'/electron)')
 1688 FORMAT( ' exc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1689 FORMAT( ' pxc-cor (loc. basis):',E19.10,' (',E15.5,'/electron)')
 1690 FORMAT( ' ion-core            :',E19.10,' (',E15.5,'/ion)')
 1691 FORMAT( ' V_local PAW residual:',E19.10,' (',E15.5,'/electron)')

 1700 FORMAT(/' QM/MM-pol-vib/CAV Energies')
 1701 FORMAT( ' --------------------------')
 1702 FORMAT( ' LJ energy                        :',E19.10)
 1703 FORMAT( ' Residual Coulomb energy          :',E19.10)
 1704 FORMAT( ' MM Vibrational energy            :',E19.10)
 1705 FORMAT( ' QM/MM coupling energy            :',E19.10)
 1706 FORMAT( ' (QM+MM)/Cavity energy            :',E19.10)
 1707 FORMAT( ' - MM Charge Field/QM Electron    :',E19.10)
 1708 FORMAT( ' - MM Charge Field/QM Ion         :',E19.10)
 1709 FORMAT( ' - MM LJ/QM LJ                    :',E19.10)
 1710 FORMAT( ' - MM Charge Field/MM Charge Field:',E19.10)
 1711 FORMAT( ' - MM LJ/MM LJ                    :',E19.10)

 1720 FORMAT(/' Dispersion energy   :',E19.10)

 1730 FORMAT(/' COSMO energies:')
 1731 FORMAT(5x,' screen=(epsilon-1)/(epsilon+1/2):',F11.6)
 1732 FORMAT(5x,' screen*Qelc*B*q energy  :',E19.10)
 1733 FORMAT(5x,' screen*Qion*B*q energy  :',E19.10)
 1734 FORMAT(5x,' screen*0.5*q*A*q energy :',E19.10)
 1735 FORMAT(5x,' solvation energy (w/o QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')
 1736 FORMAT(5x,' solvation energy (w/  QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')
c 1735 FORMAT(5x,' solvation energy (w/o QM polarization) :',E19.10,
c     >   /45x,' (',F8.3,' kcal/mol)')
c 1736 FORMAT(5x,' solvation energy (w/  QM polarization) :',E19.10,
c     >   /45x,' (',F8.3,' kcal/mol)')

 1740 FORMAT(/' extended Born solvation energies:')
 1741 FORMAT(5x,' screen=(epsilon-1)/(epsilon):',F11.6)
 1745 FORMAT(5x,' solvation energy (w/o QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')
 1746 FORMAT(5x,' solvation energy (w/  QM polarization) :',E19.10,
     >   ' (',F8.3,' kcal/mol)')



 1800 FORMAT(/' Charge+Electric Field Energies')
 1801 FORMAT( ' ------------------------------')
 1802 FORMAT( ' - Charge+Electric Field/Electron    :',E19.10)
 1803 FORMAT( ' - Charge+Electric Field/Ion         :',E19.10)
 1804 FORMAT( ' - Charge+Electric Field/Charge Field:',E19.10)
 1805 FORMAT( ' Charge+Electric Field Energy        :',E19.10)
 1806 FORMAT( ' - Electric Field/Resta_Dipole       :',E19.10)
 1807 FORMAT( ' - K.S. Electric Field/Resta_Dipole  :',E19.10)
   
      end


      subroutine psi_printeigs_debug()
      implicit none

#include "stdio.fh"

      integer MASTER,taskid
      parameter (MASTER=0)

      integer i,NN
      real*8 EV

      integer  psi_ne
      external psi_ne
      real*8   psi_eigenvalue,psi_occupation
      external psi_eigenvalue,psi_occupation

      NN=psi_ne(1)-psi_ne(2)
      EV=27.2116d0

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         do i=1,NN
            write(luout,1511) psi_eigenvalue(1,i),
     >                        psi_eigenvalue(1,i)*EV,
     >                        psi_occupation(1,i)
         end do
         do i=1,psi_ne(2)
            write(luout,1511)  psi_eigenvalue(1,i+NN),
     >                         psi_eigenvalue(1,i+NN)*EV,
     >                         psi_occupation(1,i+NN),
     >                         psi_eigenvalue(2,i),
     >                         psi_eigenvalue(2,i)*EV,
     >                         psi_occupation(2,i)
         end do
      end if

      return
 1511 FORMAT(2(E18.7,' (',F8.3,'eV)  occ=',e18.12))
      end
         

*     *******************************
*     *				    *
*     *	    cgsd_energy_gradient    *
*     *				    *
*     *******************************

      subroutine cgsd_energy_gradient(G1)
      implicit none
      real*8 G1(3,*)

#include "stdio.fh"
#include "util.fh"

      logical allow_translation,lprint,mprint,dprint,debug
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
      logical  pspw_bqext,ion_disp_on
      external pspw_bqext,ion_disp_on
      logical  control_hess_model,pspw_Efield_found
      external control_hess_model,pspw_Efield_found

      !**** hack for implementing hessian model ****
      if (control_hess_model()) then
         call ion_hess_force(G1)
         return
      end if

      allow_translation = control_allow_translation()
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      call Parallel_taskid(taskid)
      debug  = control_print(print_debug)
      dprint = ((taskid.eq.MASTER).and.control_print(print_debug))
      mprint = ((taskid.eq.MASTER).and.control_print(print_high))
      lprint = ((taskid.eq.MASTER).and.control_print(print_medium))

c*     **** debug ****
c       if (debug) then
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          call psi_1force_local(G1)
c          if (dprint) call write_force(nion,G1,"local")
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          call psi_1force_nonlocal(G1)
c          if (dprint) call write_force(nion,G1,"nonlocal")
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          if (psp_semicore(0)) call electron_semicoreforce(G1)
c          if (dprint) call write_force(nion,G1,"xc")
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          if (control_version().eq.3) call ewald_f(G1)
c          if (control_version().eq.4) call ion_ion_f(G1)
c          if (dprint) call write_force(nion,G1,"ewald ")
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          if (pspw_qmmm_found()) call pspw_qmmm_LJ_fion(G1)
c          if (dprint) call write_force(nion,G1,"LJ ")
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          if (pspw_qmmm_found()) call pspw_qmmm_Q_fion(G1)
c          if (dprint) call write_force(nion,G1,"Q ")
c          call dcopy(3*nion,0.0d0,0,G1,1)
c          call psi_1force_psp_U_v_nonlocal(G1)
c          if (dprint) call write_force(nion,G1,"psp_U_v_nonlocal")
c       end if
c*     **** debug ****
      
      call dcopy(3*nion,0.0d0,0,G1,1)

!$OMP PARALLEL
      call psi_1force(G1)
      if (psp_semicore(0)) call electron_semicoreforce(G1)

      if (control_version().eq.3) call ewald_f(G1)
      if (control_version().eq.4) call ion_ion_f(G1)
!$OMP END PARALLEL
      if (pspw_qmmm_found()) call pspw_qmmm_fion(G1)
      if (ion_disp_on()) call ion_disp_force(G1)


*     **** get forces from external charges ****
      if (pspw_charge_found()) then
        if(pspw_bqext()) then
           call pspw_charge_charge_Fion(G1)
         else
           nion1 = ion_nion()
           call pspw_charge_Fion_Fcharge(G1,G1(1,nion1+1))
           call pspw_charge_Fcharge(G1(1,nion1+1))
           call rho_1Fcharge(G1(1,nion1+1))
         end if
      end if

*     **** get forces from external Efield ****
      if (pspw_Efield_found()) then
         call pspw_Efield_Fion(G1)
      end if


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
    

*       **** print out charge positions ***
        if (pspw_charge_found().and.(.not.pspw_bqext())) then
          do i=1,pspw_charge_nion()
            i1 = ion_nion() + i
            if (ion_q_FixIon(i1)) then
            write(luout,1193) i1,pspw_charge_aname(i),
     >                    (pspw_charge_rion(K,i),K=1,3),
     >                    pspw_charge_charge(i)
            else
            write(luout,1192) i1,pspw_charge_aname(i),
     >                    (pspw_charge_rion(K,i),K=1,3),
     >                    pspw_charge_charge(i)
            end if
          end do
        end if
        end if


        if (lprint) then
        write(luout,1421)
        write(luout,1190)(i,ion_aname(I),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())

*       **** print out charge forces ***
        if (pspw_charge_found().and.(.not.pspw_bqext())) then
          do i=1,pspw_charge_nion()
            i1 = ion_nion() + i
            write(luout,1190) i1,pspw_charge_aname(i),
     >                    (G1(K,i1),K=1,3)
          end do
        end if

        write(luout,1210) fmx2,fmy2,fmz2  
        write(luout,1425)
        write(luout,1426) dsqrt(GG),
     >                    dsqrt(GG)/dble(nion),
     >                    fmax,fmax*(27.2116d0/0.529177d0)
        end if
      end if

c     call dscal(3*nion,(-1.0d0),G1,1)

      return
 1190 FORMAT(5X, I4, A5,  ' (',3F12.6,' )')
 1191 FORMAT(5X, I4, A5,  ' (',3F12.6,' ) - fixed')
 1192 FORMAT(5X, I4, A5,  ' (',3F12.6,' ) q=',F8.3)
 1193 FORMAT(5X, I4, A5,  ' (',3F12.6,' ) q=',F8.3,' - fixed')
 1210 FORMAT(5X,'   C.O.M.',' (',3F12.6,' )')
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


*     ***************************
*     *				*
*     *	    cgsd_energy_stress	*
*     *				*
*     ***************************

      subroutine cgsd_energy_stress(stress,lstress)
      implicit none
      real*8 stress(3,3)
      real*8 lstress(6)

#include "stdio.fh"
#include "util.fh"

      integer taskid,MASTER
      parameter (MASTER=0)

*     **** local variables ****
      logical oprint,hprint
      integer u,v,s
      real*8  tstress(3,3),ht(3,3),scal
      real*8  sigma(3,3),xcstress(3,3)

*     **** external functions ****
      logical  psp_semicore,pspw_SIC,pspw_HFX,control_print,ion_disp_on
      integer  control_optimize_cell_strategy
      integer  control_optimize_lattice_vectors
      integer  control_optimize_lattice
      real*8   lattice_omega,lattice_unita
      external psp_semicore,pspw_SIC,pspw_HFX,control_print,ion_disp_on
      external control_optimize_cell_strategy
      external control_optimize_lattice_vectors
      external control_optimize_lattice
      external lattice_omega,lattice_unita


      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))
      hprint = ((taskid.eq.MASTER).and.control_print(print_high))

      call dcopy(9,0.0d0,0,stress,1)

*     **** Kinetic energy component : dE_kin/dhuv ****
      call psi_1ke_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"ke")

*     **** Coulomb energy component : dE_Coul/dhuv ****
      call psi_1coulomb_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"coulomb")

*     **** Local pseudo energy component : dE_local/dhuv ****
      call dng_1vlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"local")

*     **** Nonlocal pseudo energy component : dE_nolocal/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      call psi_1vnonlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"nonlocal")

*     **** xc energy component : dE_xc/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      call rho_1exc_stress(tstress)
      call dcopy(9,tstress,1,xcstress,1)
      if (hprint) call write_stress(tstress,"xc1")

*     **** Core-correction Coulomb energy component : dE_core/dhuv ****
      if (psp_semicore(0)) then
        call rho_1semicore_stress(tstress)
        if (hprint) call write_stress(tstress,"xc2")
        call daxpy(9,1.0d0,tstress,1,xcstress,1)
      end if

      call daxpy(9,1.0d0,xcstress,1,stress,1)
      if (hprint) call write_stress(xcstress,"xc")

*     **** Ewald energy component : dE_ewald/dhuv ****
      call ewald_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"ewald")

*     **** SIC energy component : dE_SIC/dhuv ****
      if (pspw_SIC())  then
        call electron_SIC_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
        if (hprint) call write_stress(tstress,"SIC")
      end if

*     **** HFX energy component : dE_HFX/dhuv ****
      if (pspw_hfx())  then
        call electron_HFX_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
        if (hprint) call write_stress(tstress,"HFX")
      end if

*     **** DFT+U energy component : dE_dftu/dhuv ****


*     **** Dispersion energy component : dE_disp/dhuv ****
      if (ion_disp_on()) then
        call ion_disp_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
        if (hprint) call write_stress(tstress,"Dispersion")
      end if


**     **** define ht ****
       do v=1,3
       do u=1,3
          ht(u,v)=lattice_unita(v,u)
       end do
       end do

      call cell_lattice_gradient(stress,lstress)

*      **** apply lattice_vector constraints ****
       call zero_cell_constraint(stress,lstress)
 
**     **** define stress tensor ****
       call dcopy(9,0.0d0,0,sigma,1)
       scal = -1.0d0/lattice_omega()
       do v=1,3
       do u=1,3
          do s=1,3
             sigma(u,v) = sigma(u,v) + scal*stress(u,s)*ht(s,v)
          end do
       end do
       end do

c     *** call dscal(9,-1.0d0,stress,1)
      if (oprint) call write_stress(stress,"total gradient")

      if (oprint) then
        write(luout,*)
        write(luout,1811) "dE/da     = ",lstress(1)
        write(luout,1811) "dE/db     = ",lstress(2)
        write(luout,1811) "dE/dc     = ",lstress(3)
        write(luout,1811) "dE/dalpha = ",lstress(4)
        write(luout,1811) "dE/dbeta  = ",lstress(5)
        write(luout,1811) "dE/dgamma = ",lstress(6)
        write(luout,*)
      end if

      if (hprint) call write_stress(sigma,"Internal Stress Tensor")

      return
 1800 FORMAT('Lattice Constraint: gradient(',i1,',',i1,')-->0')
 1801 FORMAT('Lattice Constraint: gradient(',i1,')-->0')
 1811 FORMAT(5X,A14,F11.5)
      end

      subroutine write_stress(stress,name)
      implicit none
      real*8 stress(3,3)
      character*(*) name

#include "stdio.fh"

      integer MASTER,taskid
      parameter (MASTER=0)
      real*8 pressure,autoMbar,autoGPa,autoatm
      parameter (autoMbar=294.214239071d0)
      parameter (autoGPa=autoMbar*100.0d0)
      parameter (autoatm =290.360032539d6)

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         write(luout,1900)
         write(luout,1901)
         write(luout,1902)
         write(luout,1904)
         write(luout,1910) name
         write(luout,1911) stress(1,1),stress(1,2),stress(1,3)
         write(luout,1912) stress(2,1),stress(2,2),stress(2,3)
         write(luout,1912) stress(3,1),stress(3,2),stress(3,3)
         write(luout,1915)
         write(luout,1916) dsqrt(
     >                  stress(1,1)**2+stress(1,2)**2+stress(1,3)**2
     >                 +stress(2,1)**2+stress(2,2)**2+stress(2,3)**2
     >                 +stress(3,1)**2+stress(3,2)**2+stress(3,3)**2)
         pressure = (stress(1,1) + stress(2,2) + stress(3,3))/3.0d0
         write(luout,1917) pressure,
     >                 pressure*autoMbar,
     >                 pressure*autoGPa,
     >                 pressure*autoatm

 1900 FORMAT(//'======================')
 1901 FORMAT(  '= Stress calculation =')
 1902 FORMAT(  '======================')
 1904 FORMAT(/)
 1910 FORMAT(10X,'============= ',A,' ==============')
 1911 FORMAT(5X,' S = ',' (',3F11.5,' )')
 1912 FORMAT(5X,'     ',' (',3F11.5,' )')
 1915 FORMAT(10X,'===================================================')
 1916 FORMAT(10X,'|S|      = ',E11.5)
 1917 FORMAT(10X,'pressure = ',E9.3,' au',
     >      /10X,'         = ',E9.3,' Mbar'
     >      /10X,'         = ',E9.3,' GPa',
     >      /10X,'         = ',E9.3,' atm'/)
      end if

      return
      end



      subroutine write_lstress(lstress,name)
      implicit none
      real*8 lstress(6)
      character*(*) name

#include "stdio.fh"

      integer MASTER,taskid
      parameter (MASTER=0)
      real*8 pressure,autoMbar,autoGPa,autoatm
      parameter (autoMbar=294.214239071d0)
      parameter (autoGPa=autoMbar*100.0d0)
      parameter (autoatm =290.360032539d6)

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         write(luout,1900)
         write(luout,1901)
         write(luout,1902)
         write(luout,1904)
         write(luout,1910) name
         write(luout,1811) "dE/da     = ",lstress(1)
         write(luout,1811) "dE/db     = ",lstress(2)
         write(luout,1811) "dE/dc     = ",lstress(3)
         write(luout,1811) "dE/dalpha = ",lstress(4)
         write(luout,1811) "dE/dbeta  = ",lstress(5)
         write(luout,1811) "dE/dgamma = ",lstress(6)
         write(luout,1915)
         write(luout,*)

 1811 FORMAT(5X,A14,F11.5)
 1900 FORMAT(//'=======================')
 1901 FORMAT(  '= LStress calculation =')
 1902 FORMAT(  '=======================')
 1904 FORMAT(/)
 1910 FORMAT(10X,'============= ',A,' ==============')
 1911 FORMAT(5X,' S = ',' (',3F11.5,' )')
 1912 FORMAT(5X,'     ',' (',3F11.5,' )')
 1915 FORMAT(10X,'===================================================')
 1916 FORMAT(10X,'|S|      = ',E11.5)
 1917 FORMAT(10X,'pressure = ',E9.3,' au',
     >      /10X,'         = ',E9.3,' Mbar'
     >      /10X,'         = ',E9.3,' GPa',
     >      /10X,'         = ',E9.3,' atm'/)
      end if

      return
      end




*     ***********************************************
*     *                                             *
*     *            cgsd_pressure_stress             *
*     *                                             *
*     ***********************************************

      subroutine cgsd_pressure_stress(ispin,ne,psi,rhoall,dng,xcp,
     >                                evnl,exc,pxc,
     >                                pressure,p1,p2,stress)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi(*),rhoall(*),dng(*)
      real*8     xcp(*)
      real*8     evnl,exc,pxc
      real*8     pressure,p1,p2 
      real*8     stress(3,3)

      integer taskid,MASTER
      parameter (MASTER=0)

*     **** local variables ****
      integer u,v,s,ii,gga,nfft3d
      real*8  tstress(3,3),ht(3,3),scal,scal2
      real*8  sigma(3,3),xcstress(3,3),hm(3,3),pi

*     **** external functions ****
      logical  psp_semicore,pspw_SIC,pspw_HFX,ion_disp_on
      integer  ion_nion,control_gga
      real*8   lattice_omega,lattice_unita,ion_vion,ion_amass
      real*8   lattice_unitg
      external psp_semicore,pspw_SIC,pspw_HFX
      external ion_nion,control_gga,ion_disp_on
      external lattice_omega,lattice_unita,ion_vion,ion_amass
      external lattice_unitg


      call dcopy(9,0.0d0,0,stress,1)

*     **** Kinetic energy component : dE_kin/dhuv ****
      call ke_euv(ispin,ne,psi,tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c      p1 = (tstress(1,1)+tstress(2,2)+tstress(3,3))/3.0d0
c      write(*,*) "ke p=",p1

*     **** Coulomb energy component : dE_Coul/dhuv ****
      call coulomb_euv(dng,tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c      p1 = (tstress(1,1)+tstress(2,2)+tstress(3,3))/3.0d0
c      write(*,*) "coulomb p=",p1

*     **** Local pseudo energy component : dE_local/dhuv ****
      call v_local_euv(dng,tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c      p1 = (tstress(1,1)+tstress(2,2)+tstress(3,3))/3.0d0
c      write(*,*) "local p=",p1

*     **** Nonlocal pseudo energy component : dE_nolocal/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      call v_nonlocal_euv_2(ispin,ne,psi,tstress) 

*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do
      do v=1,3
      do u=1,3
         tstress(u,v) = tstress(u,v) - evnl*hm(u,v)
      end do
      end do
      call daxpy(9,1.0d0,tstress,1,stress,1)
c      p1 = (tstress(1,1)+tstress(2,2)+tstress(3,3))/3.0d0
c      write(*,*) "nolocal p=",p1


*     **** xc energy component : dE_xc/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)

*     **** LDA part ****
      do v=1,3
      do u=1,3
         tstress(u,v) = (exc-pxc)*hm(u,v)
      end do
      end do

*     **** PBE96 GGA part ****
      gga = control_gga()
      call D3dB_nfft3d(1,nfft3d)
      if ((gga.ge.10).and.(gga.lt.100)) then
       call v_bwexc_euv(gga,2*nfft3d,ispin,rhoall,
     >                  1.0d0,1.0d0,xcstress)
       do v=1,3
       do u=1,3
          tstress(u,v) = tstress(u,v) + xcstress(u,v)
       end do
       end do
      end if

      if (gga.eq.110) then
       call v_bwexc_euv(10,2*nfft3d,ispin,rhoall,
     >                  0.75d0,1.0d0,xcstress)
       do v=1,3
       do u=1,3
          tstress(u,v) = tstress(u,v) + xcstress(u,v)
       end do
       end do
      end if

      if (gga.eq.112) then
       call v_bwexc_euv(12,2*nfft3d,ispin,rhoall,
     >                  0.75d0,1.0d0,xcstress)
       do v=1,3
       do u=1,3
          tstress(u,v) = tstress(u,v) + xcstress(u,v)
       end do
       end do
      end if

      if (gga.eq.114) then
       call v_bwexc_euv(14,2*nfft3d,ispin,rhoall,
     >                  1.00d0,1.0d0,xcstress)
       do v=1,3
       do u=1,3
          tstress(u,v) = tstress(u,v) + xcstress(u,v)
       end do
       end do
      end if


      if (gga.eq.115) then
       call v_bwexc_euv(15,2*nfft3d,ispin,rhoall,
     >                  0.80d0,1.0d0,xcstress)
       do v=1,3
       do u=1,3
          tstress(u,v) = tstress(u,v) + xcstress(u,v)
       end do
       end do
      end if


      call dcopy(9,tstress,1,xcstress,1)
c      p1 = (tstress(1,1)+tstress(2,2)+tstress(3,3))/3.0d0
c      write(*,*) "not finished, xc p=",p1

*     **** Core-correction Coulomb energy component : dE_core/dhuv ****
      if (psp_semicore(0)) then
        call semicore_euv_vxc(ispin,xcp,tstress)
        call daxpy(9,1.0d0,tstress,1,xcstress,1)
      end if

      call daxpy(9,1.0d0,xcstress,1,stress,1)

*     **** Ewald energy component : dE_ewald/dhuv ****
      call ewald_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
c      p1 = (tstress(1,1)+tstress(2,2)+tstress(3,3))/3.0d0
c      write(*,*) "ewald p=",p1

*     **** SIC energy component : dE_SIC/dhuv ****
      if (pspw_SIC())  then
        call electron_SIC_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
      end if

*     **** HFX energy component : dE_HFX/dhuv ****
      if (pspw_hfx())  then
        call electron_HFX_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
      end if

*     **** Dispersion energy component : dE_disp/dhuv ****
      if (ion_disp_on()) then
         call ion_disp_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
      end if


**     **** define ht ****
       do v=1,3
       do u=1,3
          ht(u,v)=lattice_unita(v,u)
       end do
       end do

 
**     **** define stress tensor ****
       call dcopy(9,0.0d0,0,sigma,1)
       scal = -1.0d0/lattice_omega()
       do v=1,3
       do u=1,3
          do s=1,3
             sigma(u,v) = sigma(u,v) + scal*stress(u,s)*ht(s,v)
          end do
       end do
       end do
       p1 = (sigma(1,1)+sigma(2,2)+sigma(3,3))/3.0d0


       scal = 1.0d0/lattice_omega()
       do ii=1,ion_nion()
         scal2 = ion_amass(ii)*scal
         do v=1,3
         do u=1,3
           sigma(u,v) = sigma(u,v) 
     >                + scal2*(ion_vion(u,ii)*ion_vion(v,ii))
         end do
         end do
       end do
       pressure = (sigma(1,1)+sigma(2,2)+sigma(3,3))/3.0d0
       p2 = pressure - p1

      return
      end


*     *******************************
*     *				    *
*     *	  cgsd_energy_gradient_md   *
*     *				    *
*     *******************************
      subroutine cgsd_energy_gradient_md(G1)
      implicit none
      real*8 G1(3,*)

      integer i,k,nion,nion1
      integer i1
      real*8  GG

*     **** external functions ****
      logical     psp_semicore,pspw_charge_found,pspw_qmmm_found
      character*4 ion_aname,pspw_charge_aname
      integer     ion_katm,ion_nion,control_version
      integer     pspw_charge_nion
      real*8      ion_rion,pspw_charge_rion
      real*8      pspw_charge_charge
      external psp_semicore,pspw_charge_found,pspw_qmmm_found
      external ion_aname,pspw_charge_aname
      external ion_katm,ion_nion,control_version
      external pspw_charge_nion
      external ion_rion,pspw_charge_rion
      external pspw_charge_charge
      logical  pspw_bqext,ion_disp_on
      external pspw_bqext,ion_disp_on
      logical  control_hess_model,pspw_Efield_found
      external control_hess_model,pspw_Efield_found

      !**** hack for implementing hessian model ****
      if (control_hess_model()) then
         call ion_hess_force(G1)
         return
      end if

      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      call dcopy(3*nion,0.0d0,0,G1,1)

      call psi_1force(G1)
      if (psp_semicore(0)) call electron_semicoreforce(G1)

      if (control_version().eq.3) call ewald_f(G1)
      if (control_version().eq.4) call ion_ion_f(G1)
      if (pspw_qmmm_found()) call pspw_qmmm_fion(G1)
      if (ion_disp_on()) call ion_disp_force(G1)

      if (pspw_charge_found()) then
        if(pspw_bqext()) then
           call pspw_charge_charge_Fion(G1)
         else
           nion1 = ion_nion()
           call pspw_charge_Fion_Fcharge(G1,G1(1,nion1+1))
           call pspw_charge_Fcharge(G1(1,nion1+1))
           call rho_1Fcharge(G1(1,nion1+1))
         end if
      end if

*     **** get forces from external Efield ****
      if (pspw_Efield_found()) then
         call pspw_Efield_Fion(G1)
      end if

*     **** remove ion forces using ion_FixIon ****
      call ion_FixIon(G1)

      return
      end
