*
* $Id$
*

***********************************************************************
*                      cgmdv5                                         *

*     This is a developing Born-Oppenheimer MD code for NWCHEM        *
*                                                                     *
***********************************************************************

      logical function cgmdv5(rtdb,flag)
      implicit none
      integer rtdb
      integer flag

#include "global.fh"
#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"

#include "nwpw_timing.fh"

*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)
      real*8 autoatm
      parameter (autoatm =290.360032539d6)

      
*     **** parallel variables ****
      integer  taskid,np,np_i,np_j
      integer  MASTER
      parameter(MASTER=0)

*     **** timing variables ****
      real*8   cpu1,cpu2,cpu3,cpu4
      real*8   t1,t2,t3,t4,av

*     **** lattice variables ****
      integer ngrid(3),nwave,nfft3d,n2ft3d
      real*8  a,b,c,alpha,beta,gamma
    

*     ***** energy variables ****
      integer ne(2),ms,ispin
      real*8  E(40),en(2)
      real*8  dipole(3),cdipole(3)
      real*8  stress(3,3),lstress(6)

*     **** gradient variables ****
      integer fion(2),fion1(2)

*     **** error variables ****
      logical value,newpsi
      integer ierr

*     **** local variables ****
      logical verlet,mulliken,SA,calc_pressure,nose,bo_cpmd,vverlet
      logical oprint,lprint,hprint,qmmm,found,notfirststep,found_bak
      real*8  gx,gy,gz,cx,cy,cz
      real*8  vgx,vgy,vgz,vcx,vcy,vcz,ekg,eki0,eki1,dt,fmass,dte
      real*8  sum,w,emotion_time_shift,wa,wr
      real*8  EV,EV0,pi,esum1,esum2,eave,evar,dr,E25,E26
      real*8  ratio,aratio
      real*8  icharge,cv
      real*8  sa_alpha(2),sa_decay(2),tt1,tt2,ssr,r,p1,p2,pressure
      real*8 Te_new,Tr_new,Te_init,Tr_init
      integer i,k,ia,nion,minimizer,mapping,icount_shift,icount
      integer nbq,mapping1d,it_in,it_out

      character*50 filename
      character*255 full_filename,full_bak


*     **** external functions ****
      real*8      psp_zv,psp_rc,ewald_rcut
      real*8      ewald_mandelung
      real*8      lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      real*8      lattice_unitg,ion_amass
      integer     ewald_ncut,ewald_nshl3d
      integer     psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      integer     psp_nprj,psp_psp_type
      character*4 ion_atom,ion_aname
      external    psp_zv,psp_rc,ewald_rcut
      external    ewald_mandelung
      external    lattice_omega,lattice_unita,lattice_ecut,lattice_wcut
      external    lattice_unitg,ion_amass
      external    ewald_ncut,ewald_nshl3d
      external    psp_lmmax,psp_lmax,psp_locp,ion_nkatm
      external    psp_nprj,psp_psp_type
      external    ion_atom,ion_aname
      
      external      psp_comment
      character*255 psp_comment,comment

      real*8   control_tole,control_tolc,control_tolr,ion_rion,ion_vion
      external control_tole,control_tolc,control_tolr,ion_rion,ion_vion
      real*8   control_time_step,control_fake_mass,control_bo_fake_mass
      external control_time_step,control_fake_mass,control_bo_fake_mass
      logical  control_read,control_move,ion_init,control_out_of_time
      external control_read,control_move,ion_init,control_out_of_time
      logical  control_translation,control_rotation,control_parallel_io
      external control_translation,control_rotation,control_parallel_io
      integer  pack_nwave_all
      integer  control_it_in,control_it_out,control_gga,control_version
      integer  control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      integer  pspw_charge_nion
      external pack_nwave_all
      external control_it_in,control_it_out,control_gga,control_version
      external control_ngrid,pack_nwave,ion_nion,ion_natm,ion_katm
      external pspw_charge_nion

      character*12 control_boundry
      external     control_boundry
      character*50 control_cell_name
      external     control_cell_name

      logical      psp_semicore
      real*8       psp_rcore,psp_ncore,psp_rlocal
      external     psp_semicore
      external     psp_rcore,psp_ncore,psp_rlocal

      logical  psi_initialize, psi_finalize
      integer  psi_ispin,psi_ne,electron_count,linesearch_count
      external psi_initialize, psi_finalize
      external psi_ispin,psi_ne,electron_count,linesearch_count
      real*8   cgsd_energy,cgsd_noit_energy
      external cgsd_energy,cgsd_noit_energy
      logical  control_Mulliken,control_DOS
      external control_Mulliken,control_DOS
      real*8   ion_TotalCharge
      external ion_TotalCharge
      logical  control_check_charge_multiplicity
      external control_check_charge_multiplicity
      logical   pspw_charge_found,ion_q_FixIon,ion_q_xyzFixIon
      external  pspw_charge_found,ion_q_FixIon,ion_q_xyzFixIon
      integer  control_minimizer,control_scf_algorithm
      external control_minimizer,control_scf_algorithm
      integer  control_ks_algorithm
      external control_ks_algorithm
      real*8   control_ks_alpha,control_kerker_g0
      external control_ks_alpha,control_kerker_g0
      logical  control_print,control_balance
      external control_print,control_balance
      integer  control_mapping,control_mapping1d,control_np_orbital
      external control_mapping,control_mapping1d,control_np_orbital
      integer  control_ks_maxit_orb,control_ks_maxit_orbs
      external control_ks_maxit_orb,control_ks_maxit_orbs
      character*14 ion_q_xyzFixIon_label
      external     ion_q_xyzFixIon_label

      logical  pspw_bqext,control_bo_cpmd
      external pspw_bqext,control_bo_cpmd

      integer  control_bo_steps_in,control_bo_steps_out
      integer  control_bo_algorithm
      integer  ion_scaling_atoms_natms,ion_scaling_atoms
      real*8   control_bo_time_step,control_rti,ion_ke,ion_Temperature
      external control_bo_steps_in,control_bo_steps_out 
      external control_bo_algorithm
      external ion_scaling_atoms_natms,ion_scaling_atoms
      external control_bo_time_step,control_rti,ion_ke,ion_Temperature
      real*8   ion_com_Temperature
      external ion_com_Temperature
      logical  control_Nose,control_SA,control_pressure,Nose_restart
      real*8   control_Nose_Te,control_Nose_Tr,control_SA_decay,Nose_ssr
      real*8   Nose_Pr,Nose_Pe,Nose_Ee0,Nose_Er0,Nose_dXr
      real*8   Nose_Tr,Nose_Te,Nose_Qr,Nose_Qe,Nose_r_energy
      integer  Nose_Mchain,Nose_Nchain
      external control_Nose,control_SA,control_pressure,Nose_restart
      external control_Nose_Te,control_Nose_Tr,control_SA_decay,Nose_ssr
      external Nose_Pr,Nose_Pe,Nose_Ee0,Nose_Er0,Nose_dXr
      external Nose_Tr,Nose_Te,Nose_Qr,Nose_Qe,Nose_r_energy
      external Nose_Mchain,Nose_Nchain
      integer  ion_nconstraints,ion_ndof
      external ion_nconstraints,ion_ndof
      logical  ion_makehmass2,control_periodic_dipole
      external ion_makehmass2,control_periodic_dipole
      logical  control_precondition,control_fractional
      external control_precondition,control_fractional


*                            |************|
*****************************|  PROLOGUE  |****************************
*                            |************|

      value = .true.
      pi = 4.0d0*datan(1.0d0)

      call nwpw_timing_init()
      call dcopy(30,0.0d0,0,E,1)


*     **** get parallel variables ****
      call Parallel_Init()
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      value = control_read(11,rtdb)
      if (.not. value) 
     > call errquit('error reading control',0, INPUT_ERR)

      call Parallel2d_Init(control_np_orbital())
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)

      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
      lprint= ((taskid.eq.MASTER).and.control_print(print_low))
      hprint= ((taskid.eq.MASTER).and.control_print(print_high))


      if (oprint) call current_second(cpu1)

*     ***** print out header ****
      if (oprint) then
         write(luout,1000)
         write(luout,1010)
         write(luout,1020)
         write(luout,1010)
         write(luout,1030)
         write(luout,1031)
         write(luout,1010)
         write(luout,1035)
         write(luout,1010)
         write(luout,1040)
         write(luout,1010)
         write(luout,1041)
         write(luout,1010)
         write(luout,1000)
         call nwpw_message(1)
         write(luout,1110)
      end if
      
      ngrid(1) = control_ngrid(1)
      ngrid(2) = control_ngrid(2)
      ngrid(3) = control_ngrid(3)
      nwave = 0
      minimizer = control_minimizer()
      mapping   = control_mapping()

*     **** initialize pspw_director ****
      call pspw_director_init(rtdb)

*     **** initialize psi_data ****
      call psi_data_init(100)


*     **** initialize D3dB data structure ****
      call D3dB_Init(1,ngrid(1),ngrid(2),ngrid(3),mapping)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     ***** Initialize double D3dB data structure ****
      if (control_version().eq.4)
     >   call D3dB_Init(2,2*ngrid(1),2*ngrid(2),2*ngrid(3),mapping)

 
*     **** initialize lattice data structure ****
      call lattice_init()
      call G_init()
      call mask_init()
      call Pack_init()
      call D3dB_pfft_init()
      call ga_sync()

*     **** read ions ****
      value = ion_init(rtdb)

*     **** initialize FixIon constraint ****
      call ion_init_FixIon(rtdb)

*     **** allocate psp data structure and read in psedupotentials into it ****
      call psp_init()
      call psp_readall()
      if (psp_semicore(0)) call semicore_check()

*     **** initialize G,mask,ke,coulomb,and ewald data structures ****
      call ke_init()
      if (control_version().eq.3) call coulomb_init()
      if (control_version().eq.4) call coulomb2_init()
      call strfac_init()
      call phafac()
      if (control_version().eq.3) then
         call ewald_init()
         call ewald_phafac()
       end if

*     **** read in wavefunctions and initialize psi ****
      if (.not.control_check_charge_multiplicity()) then
        call psi_new()
        newpsi = .true.
      else
        newpsi = .false.
      end if

*     **** Initialize 1d mappings for ne spaces ****
      call psi_get_ne(ispin,ne)
      mapping1d   = control_mapping1d()
      call Dne_init(ispin,ne,mapping1d)


*     **** read in wavefunctions and initialize psi ****
      value = psi_initialize()


      !call D3dB_n_fft_init(1,psi_ne(1)+psi_ne(2))

*     **** electron and geodesic data structures ****
      call electron_init()
      if ((minimizer.eq.4).or.(minimizer.eq.7)) then  !debug
        call geodesic2_init()
      else
         call geodesic_init()
      end if

*     **** initialize QM/MM ****
      call pspw_init_APC(rtdb)
      call pspw_qmmm_init(rtdb)
      call pspw_charge_init(rtdb)

*     **** initialize SIC and HFX ****
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)
      call pspw_init_SIC(rtdb,ne)
      call pspw_init_HFX(rtdb,psi_ispin(),ne)

*     **** initialize DFT+U ****
      call psp_U_init()

*     **** initialize Meta GGA ****
      call nwpw_meta_gga_init(control_gga())

*     **** initialize metadynamics ****
      call meta_initialize(rtdb)




*     **** initialize pressure ****
      calc_pressure = control_pressure().and.(control_version().eq.3)
      pressure      = 0.0d0
      p1            = 0.0d0
      p2            = 0.0d0
      if (calc_pressure) then
         call psp_stress_init()
         call psp_stress_readall()
      end if


*     **** initialize frac_occ data structure ****
c      call frac_occ_init(rtdb,psi_ispin(),ne)

*     **** initialize linesearching ****
      call linesearch_init()


*     ******************************
*     **** scaling ion velocity ****
*     ******************************
      call ion_init_ke(ekg,eki0,eki1)


*     **** Initialize thermostats ****
      nose = control_Nose()
      if (nose) call Nose_Init((ne(1)+ne(2)),0.01d0)

*     **** Initialize simulated annealing ****
      SA       = .false.
      Tr_init  = 0.0d0
      sa_alpha(2) = 1.0d0
      if (control_SA()) then
         if (nose) then
            SA          = .true.
            sa_decay(2) = control_SA_decay(2)
            Tr_init     = control_Nose_Tr()
         else
            dt = control_bo_time_step()
            SA          = .false.
            sa_decay(2) = control_SA_decay(2)
            sa_alpha(2) = dexp( -(dt/control_SA_decay(2)) )
         end if
      end if

      bo_cpmd = control_bo_cpmd()
      vverlet = (control_bo_algorithm().eq.1)


*                |**************************|
******************   summary of input data  **********************
*                |**************************|
      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)
      mulliken = control_Mulliken()


      if (oprint) then
         write(luout,1111) np
         write(luout,1117) np_i,np_j
         if (mapping.eq.1) write(luout,1112) 
         if (mapping.eq.2) write(luout,1113) 
         if (mapping.eq.3) write(luout,1118) 
         if (control_balance()) then
           write(luout,1114)
         else
           write(luout,1116)
         end if
         write(luout,1115)
         if (control_parallel_io()) then
           write(luout,1123)
         else
           write(luout,1124)
         end if

         write(luout,1121) control_boundry(),control_version()
         if (psi_ispin().eq.1) write(luout,1130) "restricted"
         if (psi_ispin().eq.2) write(luout,1130) "unrestricted"
         !if (qmmm) write(luout,1122)
         if (control_fractional()) write(luout,1132)


         call v_bwexc_print(luout,control_gga())

         call pspw_print_SIC(luout)
         call pspw_print_HFX(luout)
         call nwpw_meta_gga_print(luout)
         if (ion_makehmass2()) write(luout,1135)
         write(luout,1140)
         do ia = 1,ion_nkatm()
           call psp_print(ia)
c           call psp_check_print(ia)
         end do



         icharge = -(psi_ne(1)+psi_ne(psi_ispin()))
         en(1)     = psi_ne(1)
         en(psi_ispin()) = psi_ne(psi_ispin())

c         if (fractional) then
c            icharge = 0
c            do ms=1,psi_ispin()
c            en(ms) =0.0
c            do i=1,ne(ms)
c              k = fweight(1)+(i-1)+(ms-1)*ne(1)
c              icharge = icharge - (3-psi_ispin())*dbl_mb(k)
c              en(ms) = en(ms) + dbl_mb(k)
c            end do
c            end do
c         end if
         icharge = icharge + ion_TotalCharge()
         write(luout,1159) icharge

         write(luout,1160)
         write(luout,1170) (ion_atom(K),ion_natm(K),K=1,ion_nkatm())
         !if (hprint) then
         write(luout,1180)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz

         call pspw_charge_Print(luout)

         write(luout,1181)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
         write(luout,1211) ion_nconstraints(),ion_ndof()

         !end if

         !write(6,1220) psi_ne(1),psi_ne(psi_ispin()),' (Fourier space)'
c         if (fractional) then
c          write(6,1219) en(1),en(psi_ispin()),' (   fractional)'
c          write(6,1221) psi_ne(1),psi_ne(psi_ispin()),' (Fourier space)'
c         else
          write(luout,1220) psi_ne(1),psi_ne(psi_ispin()),
     >                      ' (Fourier space)'
          write(luout,1221) psi_ne(1),psi_ne(psi_ispin()),
     >                      ' (Fourier space)'
c         end if
         write(luout,1230)
         write(luout,1233) control_cell_name()
         write(luout,1241) lattice_unita(1,1),
     >                 lattice_unita(2,1),
     >                 lattice_unita(3,1)
         write(luout,1242) lattice_unita(1,2),
     >                 lattice_unita(2,2),
     >                 lattice_unita(3,2)
         write(luout,1243) lattice_unita(1,3),
     >                 lattice_unita(2,3),
     >                 lattice_unita(3,3)
         write(luout,1244) lattice_unitg(1,1),
     >                 lattice_unitg(2,1),
     >                 lattice_unitg(3,1)
         write(luout,1245) lattice_unitg(1,2),
     >                 lattice_unitg(2,2),
     >                 lattice_unitg(3,2)
         write(luout,1246) lattice_unitg(1,3),
     >                 lattice_unitg(2,3),
     >                 lattice_unitg(3,3)
         call lattice_abc_abg(a,b,c,alpha,beta,gamma)
         write(luout,1232) a,b,c,alpha,beta,gamma
         write(luout,1231) lattice_omega()
         write(luout,1250) lattice_ecut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(0),pack_nwave(0)
         write(luout,1251) lattice_wcut(),ngrid(1),ngrid(2),ngrid(3),
     >                 pack_nwave_all(1),pack_nwave(1)
         if (control_version().eq.3) then
         write(luout,1260) ewald_rcut(),ewald_ncut()
         write(luout,1261) ewald_mandelung()
         end if

         write(luout,1270)
         if (minimizer.eq.1) write(luout,1371)
         if (minimizer.eq.2) write(luout,1372)
         if (minimizer.eq.3) write(luout,1373)
         if (minimizer.eq.4) write(luout,1374)
         if (minimizer.eq.5) write(luout,1375)
         if (minimizer.eq.6) write(luout,1376)
         if (minimizer.eq.7) write(luout,1377)
         if (minimizer.eq.8) write(luout,1378)
         if (minimizer.eq.9) write(luout,1379)
         if (minimizer.eq.10) write(luout,1380)
         write(luout,1280) control_time_step(),control_fake_mass()
         write(luout,1290) control_tole(),control_tolc()
         write(luout,1281) control_it_in()*control_it_out(),
     >                 control_it_in(),control_it_out()
         if ((minimizer.eq.5).or.(minimizer.eq.8)) then
           write(luout,1291)

           if (control_ks_algorithm().eq.-1)  then
             if (control_precondition()) then
                write(luout,1292)
     >             "block preconditioned conjugate gradient"
             else
                write(luout,1292) "block conjugate gradient"
             end if
           end if
           if (control_ks_algorithm().eq.0) then
             if (control_precondition()) then
                write(luout,1292) "preconditioned conjugate gradient"
             else
                write(luout,1292) "conjugate gradient"
             end if
           end if
           if (control_ks_algorithm().eq.1) 
     >       write(luout,1292) "rmm-diis"
           write(luout,1295) control_ks_maxit_orb(),
     >                       control_ks_maxit_orbs()
           if (control_scf_algorithm().eq.0) 
     >       write(luout,1293) "simple mixing"
           if (control_scf_algorithm().eq.1) 
     >       write(luout,1293) "Anderson potential mixing"
           if (control_scf_algorithm().eq.2) 
     >       write(luout,1293) "Johnson-Pulay mixing"
           if (control_scf_algorithm().eq.3) 
     >       write(luout,1293) "Anderson density mixing"
           write(luout,1294) control_ks_alpha()
           write(luout,1301) control_kerker_g0()
         end if
         write(luout,1310)
         if (.not.control_translation()) write(luout,1271)
         if (.not.control_rotation())    write(luout,1272)
         if (.not.bo_cpmd) then
            write(luout,1320) control_bo_time_step(),
     >       control_bo_steps_in()*control_bo_steps_out(),
     >       control_bo_steps_in(),control_bo_steps_out() 
         else
            write(luout,1321) control_bo_time_step(),
     >       control_bo_fake_mass(),
     >       control_bo_steps_in()*control_bo_steps_out(),
     >       control_bo_steps_in(),control_bo_steps_out() 
        end if

         if(control_bo_algorithm().eq.0) write(6,1330) "Position Verlet"
         if(control_bo_algorithm().eq.1) write(6,1330) "Velocity Verlet"
         if(control_bo_algorithm().eq.2) write(6,1330) "Leap Frog"

         write(luout,1340) control_rti()
         call ion_scaling_atoms_print(luout)
         write(luout,1222) eki0,ekg
         write(luout,1223) eki1
         write(luout,1224) (eki1-eki0)
         if (nose) then
           write(luout,1395)
           if (Nose_restart()) then
              write(luout,*) "    thermostats resused"
           else
              write(luout,*) "    thermostats initialized"
           end if
           do i=1,Nose_Nchain()
             write(luout,1398) i,control_Nose_Tr(),Nose_Qr(i),
     >                   Nose_Pr(i),Nose_Er0(i)
           end do
         else
           write(luout,1394)
         end if
        if (calc_pressure) write(luout,1393)
         if (control_SA()) then
           write(luout,1396) sa_decay(2)
         end if
         if (mulliken) write(luout,1399)
         write(luout,1300)
         call util_flush(luout)
      end if

*                |***************************|
******************     simple MD loop        **********************
*                |***************************|
      if (taskid.eq.MASTER) call current_second(cpu2)
      if ((taskid.eq.MASTER).and.(.not.calc_pressure)) 
     >   call nwpw_message(10)
      if ((taskid.eq.MASTER).and.(calc_pressure)) 
     >   call nwpw_message(11)


*     **** write initial position to xyz data ****
      call xyz_init()          ! unit=18
      call MOTION_init(rtdb)   ! unit=19


*     *** fei io ****
      call fei_init(rtdb)

*     **** ecce print ****
      call ecce_print_module_entry('task Born-Oppenheimer')
      !call ecce_print_module_entry('driver')
      call movecs_ecce_print_off()


*     **** finalize pressure ****
      if (calc_pressure) then
         call psp_stress_end()
      end if




*     ************************************
*     **** open up other MOTION files ****
*     ************************************
       
*     **** open EMOTION file ****
      E25 = 0.0d0
      E26 = 0.0d0
      icount_shift = 0
      if (.not.btdb_cget(rtdb,'nwpw:emotion_filename',1,filename))
     >  call util_file_prefix('emotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) then

*        **** check for backup file ****
         call util_file_name_noprefix('EMOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(*,*)
            write(*,*) "EMOTION99-bak exists:"
            i=index(full_bak,' ')
            k=index(full_filename,' ')
            write(*,*) "   Copying ",full_bak(1:i),
     >                 " to ",full_filename(1:k)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if

         emotion_time_shift = 0.0d0
         icount_shift       = 0
         inquire(file=full_filename,exist=found)
         if (found) then
           open(unit=31,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
           read(31,*,end=100) emotion_time_shift,w,sum
           E25 = E25 + sum   !*** take care of running sums ***
           E26 = E26 + sum*sum
           icount_shift = icount_shift + 1
           end do
  100      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(GCC46)
           backspace 31
#endif
         else
           open(unit=31,file=full_filename,form='formatted',
     >          status='new')
         end if
      end if

*     **** open EIGMOTION file ****
      if (mulliken) then
        if (.not.btdb_cget(rtdb,'nwpw:eigmotion_filename',1,filename))
     >    call util_file_prefix('eigmotion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)
      if (taskid.eq.MASTER)
     >   open(unit=32,file=full_filename,form='formatted')
      end if

      call xyz_write()

*     ****  allocate fion ***
      nion = ion_nion()
      if (pspw_charge_found().and.
     >    (.not.pspw_bqext())) nion = nion + pspw_charge_nion()

      value = BA_push_get(mt_dbl,(3*nion),
     >                       'fion',fion(2),fion(1))
      if (vverlet) then
         value = value.and.BA_push_get(mt_dbl,(3*nion),
     >                                'fion1',fion1(2),fion1(1))
      end if
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      if (bo_cpmd) call psi_cpmd_start()

*     ***** do a simple md loop ***
      notfirststep = .false.
      dt = control_bo_time_step()
      fmass = control_bo_fake_mass()
      dte = dt*dt/fmass
      call control_reduce_print()
      EV = cgsd_energy(newpsi)
      call cgsd_energy_gradient_md(dbl_mb(fion(1)))

      r = 1.0d0
      if (nose) r =  (1.0d0-0.5d0*dt*Nose_dXr())
      call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2)*r) 
      eki1 = ion_ke()
      if (nose) call Nose_Newton_Step(0.01d0,eki1)

      call xyz_write()

      it_out = control_bo_steps_out()
      it_in  = control_bo_steps_in()
      icount = 0
      if (it_out.lt.1) goto 102


  101 continue
         icount = icount + 1

c       **** inner loop ****
        do i=1,it_in
           if (vverlet) then
              call Parallel_shared_vector_copy(.true.,3*nion,
     >                                  dbl_mb(fion(1)),
     >                                  dbl_mb(fion1(1)))
              call ion_shift21() 
           else
              call ion_shift() 
              if (nose) call Nose_shift()
           end if

           if (bo_cpmd) then
              if (notfirststep) then
                 call psi_cpmd_step(dte)
              else
                 notfirststep = .true.
                 call psi_cpmd_step(0.5d0*dte)
              end if
           end if

           EV = cgsd_energy(.false.)
           call cgsd_energy_gradient_md(dbl_mb(fion(1)))


           if (nose) then
              ssr = Nose_ssr()
              call ion_nose_step(ssr,dbl_mb(fion(1)))
              eki1 = ion_ke()
              call Nose_Verlet_Step(0.01d0,eki1)
           else if (vverlet) then
              call ion_vverlet_step(dbl_mb(fion(1)),
     >                              dbl_mb(fion1(1)))  
              call ion_vshift()  
              call ion_newton_step(dbl_mb(fion(1)),sa_alpha(2))
              eki1 = ion_ke()
           else
              call ion_verlet_step(dbl_mb(fion(1)),sa_alpha(2))
              eki1 = ion_ke()
           end if
        end do
        E25 = E25 + EV
        E26 = E26 + EV*EV


        if (calc_pressure) then
           call psi_1pressure_stress(pressure,p1,p2,stress)
        end if

        if (oprint) then
           if (calc_pressure) then
             if (nose) then
             write(luout,1350) icount*it_in,EV+eki1+Nose_r_energy(),
     >                         EV,eki1,ion_Temperature(),
     >                         pressure*autoatm
             write(31,1360) icount*it_in*dt+emotion_time_shift,
     >                      EV+eki1+Nose_r_energy(),
     >                      EV,eki1,ion_Temperature(),
     >                      pressure*autoatm
             else
             write(luout,1350) icount*it_in,EV+eki1,EV,eki1,
     >                         ion_Temperature(),
     >                         pressure*autoatm
             write(31,1360) icount*it_in*dt+emotion_time_shift,
     >                      EV+eki1,EV,eki1,ion_Temperature(),
     >                      pressure*autoatm
             end if
           else
             if (nose) then
             write(luout,1350) icount*it_in,EV+eki1+Nose_r_energy(),
     >                         EV,eki1,ion_Temperature()
             write(31,1360) icount*it_in*dt+emotion_time_shift,
     >                      EV+eki1+Nose_r_energy(),
     >                      EV,eki1,ion_Temperature()
             else
             write(luout,1350) icount*it_in,EV+eki1,EV,eki1,
     >                         ion_Temperature()
             write(31,1360) icount*it_in*dt+emotion_time_shift,
     >                      EV+eki1,EV,eki1,ion_Temperature()
             end if
           end if
           call util_flush(6)
           call util_flush(31)
        end if

        call fei_output(EV,dbl_mb(fion(1)))
        call xyz_write()
        call MOTION_write(icount*it_in*dt)

*        **** call pspw_director ****
         call pspw_director(rtdb)

*        **** update thermostats using SA decay ****
         if (SA) then
           tt2 = icount*it_in*dt/sa_decay(2)
           Tr_new = Tr_init*dexp(-tt2)
           call Nose_reset_T(Tr_new,Tr_new)
         end if

*        **** exit early ****
         if (control_out_of_time()) then
            if (oprint)
     >       write(luout,*) ' *** out of time. iteration terminated'
            go to 102
         end if


      if (icount.lt.it_out) go to 101
      if (oprint) write(luout,*) 
     > '*** arrived at the Maximum iteration.   terminated.'

*::::::::::::::::::::  end of iteration loop  :::::::::::::::::::::::::
      call control_up_print()


  102 continue
  
      if (bo_cpmd) call psi_cpmd_end()

*     **** close xyz and MOTION files ****
      call xyz_end()
      call MOTION_end()
      if (taskid.eq.MASTER) then
         close(unit=31)
         if (mulliken) close(unit=32)

*        **** remove EMOTION backup file ***
         call util_file_name_noprefix('EMOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)

      end if

*     *** close fei io ****
      call fei_end()



c*     ***** do a simple mc loop ***
c      call control_reduce_print()
c      dr = 0.08d0
c      aratio = 0.5d0
c      wa = 0.0d0
c      wr = 0.0d0
c      esum1 = 0.0d0
c      esum2 = 0.0d0
c      eave = 0.0d0
c      evar = 0.0d0
c      beta = util_random(39) !**seed the random number generator
c      beta = 1.0d0/(kb*298.15)
c      call xyz_write()
c      EV0 = cgsd_energy(newpsi)
c      do it=1,control_bo_steps()
c
c        call ion_MC_step(dr)
c        EV = cgsd_energy(.false.)
c        w = dexp(-beta*(EV-EV0))
c        alpha = util_random(0)
c
c        !*** accept the sequence ***
c        if (alpha.lt.w) then
c          EV0 = EV
c          !write(6,*) "   ----ACCEPTED----"
c          if (taskid.eq.MASTER)  then
c             write(6,1350) it,EV,EV0,w,alpha,wa,wr,eave,evar,dr
c             call util_flush(6)
c          end if
c          call xyz_write()
c          call MOTION_write(it*dt)
c          wa = wa + 1.0d0
c          esum1 = esum1 + EV
c          esum2 = esum2 + EV*EV
c          eave = esum1/wa
c          evar = esum2/wa
c          evar = evar - eave*eave
c
c        !*** reject the sequence ***
c        else
c           call ion_MC_reject_step()
c          !write(6,*) "   ----REJECTED----"
c          wr = wr + 1.0d0
c        end if
c
c        !*** adjust ***
c        if (mod(it,10).eq.0) then
c            ratio = wa/(wa+wr)
c            dr = dr*(ratio/aratio)
c            !if (ratio.lt.aratio) dr = dr*0.95d0
c            !if (ratio.gt.aratio) dr = dr*1.05d0
c        end if
c
c      end do
c
c      call control_up_print()

*                |***************************|
******************     simple MD loop        **********************
*                |***************************|


*                |***************************|
****************** report summary of results **********************
*                |***************************|

      call center_geom(cx,cy,cz)
      call center_mass(gx,gy,gz)
      call center_v_geom(vcx,vcy,vcz)
      call center_v_mass(vgx,vgy,vgz)

      if (oprint) then
         call print_elapsed_time(icount*it_in*dt)
         write(luout,1300)
         write(luout,1410)
         write(luout,1420)
         do I=1,ion_nion()
           if (ion_q_FixIon(I)) then
           write(luout,1191) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           else if (ion_q_xyzFixIon(I)) then
           write(luout,1194) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0,ion_q_xyzFixIon_label(I)
           else
           write(luout,1190) I,ion_aname(I),(ion_rion(k,i),K=1,3),
     >                   ion_amass(I)/1822.89d0
           end if
         end do
         write(luout,1200) cx,cy,cz
         write(luout,1210) gx,gy,gz

         write(luout,1421)
         write(luout,1192) (I,ion_aname(I),
     >                  (ion_vion(K,I),K=1,3),I=1,ion_nion())
         write(luout,1200) vcx,vcy,vcz
         write(luout,1210) vgx,vgy,vgz
         write(luout,1211) ion_nconstraints(),ion_ndof()
      end if

      EV = cgsd_noit_energy()

      if (oprint) then

         write(luout,1300)
         write(luout,1460)
         if (nose) then
            cx = Nose_r_energy() 
            write(luout,1461) EV+ion_ke()+cx,(EV+ion_ke()+cx)/ion_nion()
            write(luout,1462) cx,cx/ion_nion()
         else
            write(luout,1461) EV+ion_ke(),(Ev+ion_ke())/ion_nion()
         end if
         write(luout,1472) ion_ke(),ion_ke()/ion_nion()

*       **** write out Temperatures ****
         eki0 = ion_Temperature()
         write(luout,1480) eki0
         write(luout,1490) ion_com_Temperature()

         eave = E25/dble(icount+icount_shift)
         evar = E26/dble(icount+icount_shift)
         evar = evar - eave*eave
         cv = (evar)/(kb*ion_Temperature()**2)
         cv = cv/dble(ion_nion())
         write(luout,1492) eave
         write(luout,1493) evar
         write(luout,1494) cv

      end if
      call cgsd_energy_gradient(dbl_mb(fion(1)))

*     **** calculate excited state orbitals ****
      call cgsd_excited()

*     **** calculate oep eigenvalues ****
      call cgsd_oep_eigenvalues()

*     **** extra energy output for QA test ****
      if (lprint) write(6,1600) EV

*     **** calculate the spin contamination ****
      if (flag.gt.-1) call psi_spin2(dipole(1))

*     **** calculate the dipole ***
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0
      cdipole(1) = 0.0d0
      cdipole(2) = 0.0d0
      cdipole(3) = 0.0d0
      if (flag.gt.-1) then
        call rho_dipole(dipole)
        if (control_periodic_dipole()) call psi1_crystal_dipole(cdipole)
      end if


*     **** calculate the stress tensor ****
      call dcopy(9,0.0d0,0,stress,1)
      call dcopy(6,0.0d0,0,lstress,1)
      if (flag.eq.3) then
         call psp_stress_init()
         call psp_stress_readall()
         call cgsd_energy_stress(stress,lstress)
         call psp_stress_end()
      end if


*     *************************************************************
*     **** output energy, dipole, and gradient to rtdb for use ****
*     **** by task_energy and task_gradient                    ****
*     *************************************************************
      if (flag.gt.-1) then
      value = btdb_put(rtdb,'pspw:energy',mt_dbl,1,EV)
      value = value.and.
     >        btdb_put(rtdb,'pspw:dipole',mt_dbl,
     >                 3,dipole)
      value = value.and.
     >        btdb_put(rtdb,'pspw:crystal_dipole',mt_dbl,
     >                 3,cdipole)
      end if
      if (flag.gt.0) then
        value = value.and.
     >        btdb_put(rtdb,'pspw:gradient',mt_dbl,
     >                 3*nion,dbl_mb(fion(1)))
      end if
      if (flag.eq.3) then
        value = value.and.
     >        btdb_put(rtdb,'pspw:stress',mt_dbl,
     >                 9,stress)
        value = value.and.
     >        btdb_put(rtdb,'pspw:lstress',mt_dbl,
     >                 6,lstress)
      end if
      if (vverlet) value = value.and.BA_pop_stack(fion1(2))
      value = value.and.BA_pop_stack(fion(2))
      if (.not. value) call errquit('cgmdv5: error writing rtdb',0,
     &       RTDB_ERR)

      if (taskid.eq.MASTER) call current_second(cpu3)

*                |***************************|
******************         Epilogue          **********************
*                |***************************|

*     **** calculate Mulliken Populations ***
      if (control_Mulliken()) call psi_Mulliken(rtdb)




*     **** calculate Mulliken Populations ***

*     **** write wavefunctions to file and finalize psi ****
      if (flag.eq.-1) then
        value = psi_finalize(.false.)
      else
        value = psi_finalize(.true.)
      end if

*     **** write geometry to rtdb ****
      call pspw_charge_write(rtdb)
      call ion_write(rtdb)


*     **** deallocate heap memory ****
      call electron_finalize()
      if ((minimizer.eq.4).or.(minimizer.eq.7)) then
        call geodesic2_finalize()
      else
        call geodesic_finalize()
      end if
      if (control_version().eq.3) call ewald_end()
      call strfac_end()
      if (control_version().eq.3) call coulomb_end()
      if (control_version().eq.4) call coulomb2_end()
      call ke_end()
      call mask_end()
      call Pack_end()
      call G_end()
      call psp_U_end()
      call nwpw_meta_gga_end()
      call pspw_end_SIC()
      call pspw_end_HFX()
      call pspw_end_APC()
      call pspw_charge_end()
      call pspw_qmmm_end()
c      call frac_occ_end()
      call meta_finalize(rtdb)

      if (nose) call Nose_end()
      call ion_end()
      call psp_end()
      call ion_end_FixIon()
      !call D3dB_n_fft_end(1)
      call D3dB_pfft_end()
      call D3dB_end(1)
      if (control_version().eq.4) call D3dB_end(2)
      call Dne_end()
      call psi_data_end()


*     **** do anaylysis on MOTION files ****
      call cpmd_properties(rtdb)



*                |***************************|
****************** report consumed cputime   **********************
*                |***************************|
      if (lprint) then
         CALL current_second(cpu4)

         T1=CPU2-CPU1
         T2=CPU3-CPU2
         T3=CPU4-CPU3
         T4=CPU4-CPU1
         AV=T2/dble(electron_count())
         write(6,1801)
         write(6,1802)
         write(6,1803) T1
         write(6,1804) T2
         write(6,1805) T3
         write(6,1806) T4
         write(6,1807) AV,electron_count(),linesearch_count()

         call nwpw_timing_print_final(oprint,electron_count())

         write(6,*)
         CALL nwpw_MESSAGE(4)
      end if 



      call Parallel_Finalize()
      cgmdv5 = value
      return


*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::
 1000 FORMAT(10X,'****************************************************')
 1010 FORMAT(10X,'*                                                  *')
 1020 FORMAT(10X,'*               NWPW PSPW Calculation              *')
 1030 FORMAT(10X,'*     [      Born-Oppenheimer molecular       ]    *')
 1031 FORMAT(10X,'*     [          dynamics simulation          ]    *')
 1035 FORMAT(10x,'*     [  NorthWest Chemistry implementation   ]    *')
 1040 FORMAT(10X,'*            version #1.00   07/18/06              *')
 1041 FORMAT(10X,'*    This code was developed by Eric J. Bylaska.   *')
 1100 FORMAT(//)
 1110 FORMAT(10X,'================ input data ========================')
 1111 FORMAT(/' number of processors used:',I10)
 1112 FORMAT( ' parallel mapping         :      1d-slab')
 1113 FORMAT( ' parallel mapping         :   2d-hilbert')
 1114 FORMAT( ' parallel mapping         :     balanced')
 1115 FORMAT(/' options:')
 1116 FORMAT( ' parallel mapping         : not balanced')
 1117 FORMAT( ' processor grid           :',I4,' x',I4)
 1118 FORMAT( ' parallel mapping         :    2d-hcurve')
 1120 FORMAT(5X,' ionic motion         = ',A)
 1121 FORMAT(5X,' boundary conditions  = ',A,'(version', I1,')')
 1122 FORMAT(5X,' qmmm simulation')
 1123 FORMAT( ' parallel io              :        on')
 1124 FORMAT( ' parallel io              :       off')
 1130 FORMAT(5X,' electron spin        = ',A)
 1131 FORMAT(5X,' exchange-correlation = ',A)
 1132 FORMAT(5X,' using fractional occupation')
 1135 FORMAT(/' The masses of QM H atoms converted to 2.0 amu. ',
     >       /' To turn off this default',
     >        ' set nwpw:makehmass2 .false.')
 1140 FORMAT(/' elements involved in the cluster:')
 1150 FORMAT(5X,I2,': ',A4,'  core charge:',F4.1,'  lmax=',I3)
 1151 FORMAT(5X,'        cutoff =',4F8.3)
 1152 FORMAT(12X,' highest angular component      : ',i3)
 1153 FORMAT(12X,' local potential used           : ',i3)
 1154 FORMAT(12X,' number of non-local projections: ',i3)
 1155 FORMAT(12X,' semicore corrections included  : ',
     >       F6.3,' (radius) ',F6.3,' (charge)')
 1156 FORMAT(12X,' aperiodic cutoff radius        : ',F6.3)
 1157 FORMAT(12X,' comment    : ',A)
 1158 FORMAT(12X,' pseudpotential type            : ',i3)

 1159 FORMAT(/' total charge:',F8.3)
 1160 FORMAT(/' atomic composition:')
 1170 FORMAT(7(5X,A4,':',I5))
 1180 FORMAT(/' initial position of ions (au):')
 1181 FORMAT(/' initial velocity of ions after scaling (au):')

 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - atomic mass= ',F7.3,' ')
 1191 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - fixed')
 1192 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1193 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,' - z fixed')
 1194 FORMAT(5X, I4, A5, ' (',3F11.5,
     >       ' ) - atomic mass= ',F7.3,A)

 1200 FORMAT(5X,'    G.C. ',' (',3F11.5,' )')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1211 FORMAT(5X,'   number of constraints = ', I6,' ( DOF = ',I6,' )' )
 1219 FORMAT(/' number of active electrons: spin up=',F6.2,
     >        '  down=',F6.2,A)
 1220 FORMAT(/' number of active electrons: spin up=',I6,  
     >        '  down=',I6,A)
 1221 FORMAT( ' number of active orbitals : spin up=',I6,  
     >        '  down=',I6,A)

 1222 format(5x,' initial kinetic energy=',e12.5,' (ion)',2x,
     >                                     e12.5,' (c.o.m.)')
 1223 format(5x,' after scaling=         ',e12.5,' (ion)')
 1224 format(5x,' increased energy=      ',e12.5,' (ion)')
 1226 format(/' final kinetic energy= ',   e12.5,' (ion)',2x,
     >                                     e12.5,' (c.o.m.)')

 1230 FORMAT(/' supercell:')
 1231 FORMAT(5x,'             omega=',F12.1)
 1232 FORMAT(5x,' lattice:    a=    ',f8.3,' b=   ',f8.3,' c=    ',f8.3,
     >      /5x,'             alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
 1233 FORMAT(5x,' cell_name:  ',A)
 1241 FORMAT(5x,' lattice:    a1=<',3f8.3,' >')
 1242 FORMAT(5x,'             a2=<',3f8.3,' >')
 1243 FORMAT(5x,'             a3=<',3f8.3,' >')
 1244 FORMAT(5x,' reciprocal: b1=<',3f8.3,' >')
 1245 FORMAT(5x,'             b2=<',3f8.3,' >')
 1246 FORMAT(5x,'             b3=<',3f8.3,' >')

 1250 FORMAT(/5X,' density cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
 1251 FORMAT(5X,' wavefnc cutoff=',F7.3,'  fft=',I3,'x',I3,'x',I3,
     &       '( ',I8,' waves ',I8,' per task)')
  
 1260 FORMAT(5X,' Ewald summation: cut radius=',F8.2,'  and',I3)
 1261 FORMAT(5X,'                   madelung=',f14.8)


 1270 FORMAT(/' technical parameters for minimizer:')
 1271 FORMAT(5x, ' translation constrained')
 1272 FORMAT(5x, ' rotation constrained')
 1280 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1)
 1281 FORMAT(5X, ' maximum iterations =',I10, 
     >           ' ( ',I4,' inner ',I6,' outer )')
 1290 FORMAT(5X, ' tolerance=',E8.3,' (energy)',E12.3,
     &        ' (density)')
 1291 FORMAT(/' Kohn-Sham scf parameters:')
 1292 FORMAT(5X, ' Kohn-Sham algorithm  = ',A)
 1293 FORMAT(5X, ' SCF algorithm        = ',A)
 1294 FORMAT(5X, ' SCF mixing parameter =',F7.4)
 1295 FORMAT(5X, ' Kohn-Sham iterations = ',I4,
     >           ' (',I4,' outer)')
 1300 FORMAT(//)
 1301 FORMAT(5X, ' Kerker damping       =',F7.4)
 1310 FORMAT(/' molecular dynamics parameters:')
 1320 FORMAT(5X, ' time step=',F10.2,5X,' iterations=',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1321 FORMAT(5X, ' time step=',F10.2,5X,'fictitious mass=',F10.1,
     >       /5X,' iterations=',I10,
     >           ' ( ',I4,' inner ',I6,' outer )')
 1330 FORMAT(5X, ' integration algorithm=   ',A)
 1340 FORMAT(/5X, ' cooling/heatting rate= ',e12.5,' (ion)')
 1350 FORMAT(I8,2E19.10,E14.5,4F14.2,2E19.10,2F8.4)
 1360 format(100e19.10)
 1371 FORMAT(5X,' Grassmann conjugate gradient iteration')
 1372 FORMAT(5X,' Grassmann lmbfgs iteration')
 1373 FORMAT(5X,' Grassmann diis iteration')
 1374 FORMAT(5X,' Stiefel conjugate gradient iteration')
 1375 FORMAT(5X,' Kohn-Sham scf iteration (potential)')
 1376 FORMAT(5X,'Projected Grassmann lmbfgs iteration')
 1377 FORMAT(5X,' Stiefel lmbfgs iteration')
 1378 FORMAT(5X,' Kohn-Sham scf iteration (density)')
 1379 FORMAT(5X,' Grassmann cg (Stich linesearch) iteration')
 1380 FORMAT(5X,' Grassmann lmbfgs (Stich linesearch) iteration')

 1393 format(/' Pressure Output Generated         ')
 1394 format(/' Constant Energy Simulation                     ')
 1395 format(/' Nose-Hoover Simulation - Thermostat Parameters:')
 1396 format(5x, 'SA decay rates  =',e10.3,' (ion)')
 1397 format(5x, 'link = ',I3,
     > ' Te =',f8.2,' Qe =',e10.3,' 2*pi/we=',e10.3,' Ee0=',e10.3)
 1398 format(5x, 'link = ',I3,
     > ' Tr =',f8.2,' Qr =',e10.3,' 2*pi/wr=',e10.3,' Er0=',e10.3)

 1399 format(//' Mulliken Analysis Output Generated            ')
 1400 FORMAT(I3,3E18.8/3X,3E18.8)
 1410 FORMAT(10X,'=============  summary of results  =================')
 1420 FORMAT(/' final position of ions:')
 1421 FORMAT(/' final velocity of ions:')
 1460 FORMAT( ' AIMD Energies:')
 1461 FORMAT( ' Total AIMD energy       :',E19.10,' (',E15.5,'/ion)')
 1462 FORMAT( ' Nose thermostat energy  :',E19.10,' (',E15.5,'/ion)')
 1471 FORMAT(/' Kinetic energy (elc)    :',E19.10,' (',E15.5,'/elc)')
 1472 FORMAT( ' Kinetic energy (ion)    :',E19.10,' (',E15.5,'/ion)')
 1473 FORMAT( ' thermostat energy (elc) :',E19.10,' (',E15.5,'/elc)')
 1474 FORMAT( ' thermostat energy (ion) :',E19.10,' (',E15.5,'/ion)')
 1480 FORMAT( ' Temperature             :',F19.1,' K (ion)')
 1490 FORMAT( '                         :',F19.1,' K (c.o.m.)')
 1491 FORMAT( ' Temperature             :',F19.1,' K (elc)')
 1492 FORMAT( ' Eaverage                :',E19.10)
 1493 FORMAT( ' Evariance               :',E19.10)
 1494 FORMAT( ' Cv - f*kb/(2*nion)      :',E19.10)
 1600 FORMAT(/' Total PSPW energy   :',E19.10)
 1801 FORMAT(//'== Timing ==')
 1802 FORMAT(/'cputime in seconds')
 1803 FORMAT( '  prologue    : ',E14.6)
 1804 FORMAT( '  main loop   : ',E14.6)
 1805 FORMAT( '  epilogue    : ',E14.6)
 1806 FORMAT( '  total       : ',E14.6)
 1807 FORMAT( '  cputime/step: ',E14.6,
     >        '       (',I8,' evalulations,', I8,' linesearches)')
 1808 FORMAT(A,E14.6,E14.6)
 1809 FORMAT(//A,2A14)

 9010 FORMAT(//' >> job terminated due to code =',I3,' <<')

 9000 if (taskid.eq.MASTER) write(6,9010) ierr
      call Parallel_Finalize()

      cgmdv5 = value
      return
      END
